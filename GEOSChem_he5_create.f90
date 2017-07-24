PROGRAM GEOSChem_he5_create

  USE ISO_C_BINDING
  USE GEOSChem_he5_module
  IMPLICIT NONE

  ! --------------------
  ! HE5 Output variables
  ! --------------------
  CHARACTER (LEN=132)          :: he5file
  CHARACTER (LEN=5), PARAMETER :: sensor = 'TEMPO'

  ! --------------
  ! Version number
  ! --------------
  CHARACTER (LEN=4), PARAMETER :: version = 'v0p0'

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER*4 :: imon

  ! -----------------------------------------------------
  ! Set the HE5 output file name, one file for each month
  ! -----------------------------------------------------
  DO imon = 1, 12
     WRITE(he5file,'(A)') TRIM(ADJUSTL(sensor))//'_GEOS-Chem_climatology_'&
          //SwathFieldsMonth(imon)//'_'//version//'.he5'
     print*, he5file
 
     ! ------------------------------------------------------
     ! Initialize output quantity and call HE5 output routine
     ! ------------------------------------------------------
     CALL create_he5_file ( he5file, TRIM(ADJUSTL(sensor)), &
          imon)
  END DO
  STOP

END PROGRAM GEOSChem_he5_create

SUBROUTINE create_he5_file ( he5file, sensor, month)

  !------------------------------------------------------------------------------
  ! This subroutine creates an HE5 file with Trace Gases Climatologies
  !
  ! Input:
  !   he5file    - Name of HE5 output file
  !   sensor     - Name of satellite sensor
  !   month      - Month of the climatology
  !   molecules  - Molecules to be included
  !   nlon       - Number of longitudes
  !   nlat       - Number of latitudes
  !   nlev       - Number of vertical levels
  !   nhrs       - Number of hours in the day
  !------------------------------------------------------------------------------
  USE GEOSChem_he5_module
  USE OMSAO_lininterpolation_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=*),    INTENT (IN)  :: he5file, sensor
  INTEGER*4, INTENT(IN) :: month

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER :: swath_file_id  
  INTEGER :: swath_id       
  INTEGER :: he5stat, idf, ihrs
  CHARACTER (LEN=132) :: swath_name, ascii_name
  REAL*4 :: deltap, scale

  ! ------------------------------------------------
  ! Variables to read and hold climatology data from 
  ! GEOS-Chem ASCII file (needs to modify nlon and
  ! nlat to match GEOS-Chem resolution
  ! ------------------------------------------------
  INTEGER, PARAMETER :: nlon = 57, nlat = 41, nlev = 47, &
       nhrs = 24
  REAL*4, DIMENSION (nlon,nlat,nlev,nhrs) :: GC_H2O, GC_HCHO, GC_NO2, &
       GC_O3, GC_TEMP
  REAL*4, DIMENSION (nlon) :: GC_lon
  REAL*4, DIMENSION (nlat) :: GC_lat
  REAL*4, DIMENSION (nhrs) :: GC_time
  REAL*4, DIMENSION (nlon,nlat,nhrs) :: GC_PSUR

  ! --------------------------------------------
  ! Variables to define and work out output grid
  ! --------------------------------------------
  REAL*4, PARAMETER :: dlon = 0.25, dlat = 0.25
  REAL*4, PARAMETER :: flon = -160, llon = -40, flat = 15, llat = 65
  INTEGER :: noutlon, noutlat, ideg
  REAL*4, ALLOCATABLE, DIMENSION(:) :: out_lon, out_lat, out_time

  ! ---------------------------------------
  ! Variables to compute newly gridded data
  ! ---------------------------------------
  REAL*8 :: val
  INTEGER :: status, ilon, ilat, ilev
  INTEGER*4, ALLOCATABLE, DIMENSION(:,:,:) :: out_lut_data
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: out_data_3d
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:) :: out_data_4d

  ! ----------------------------------------------------
  ! LUT ozone profile index and ozone total column
  ! Different arrays for low, medium and high latitudes.
  ! ----------------------------------------------------
  INTEGER*4, PARAMETER, DIMENSION(4) :: L_idx = (/1,2,3,4/)
  INTEGER*4, PARAMETER, DIMENSION(8) :: M_idx = (/5,6,7,8,9,10,11,12/)
  INTEGER*4, PARAMETER, DIMENSION(10) :: H_idx = (/13,14,15,16,17,18,19,20,21,22/)

  REAL*4, PARAMETER, DIMENSION(4) :: L_val = (/200,250,300,350/)
  REAL*4, PARAMETER, DIMENSION(8) :: M_val = (/200,250,300,350,400,450,500,550/)
  REAL*4, PARAMETER, DIMENSION(10) :: H_val = (/100,150,200,250,300,350,400,450,500,550/)

  ! Code starts here
  he5stat = 0

  ! ----------------------------------------
  ! Initialize GEOS-Chem variables variables
  ! ----------------------------------------
  GC_lon = fill_value; GC_lat = fill_value; GC_time = fill_value
  GC_H2O = fill_value; GC_HCHO = fill_value; GC_O3 = fill_value
  GC_NO2 = fill_value; GC_PSUR = fill_value

  ! -------------------------------------------
  ! Read data to be included in the climatology
  ! One file for each hour
  ! -------------------------------------------
  READCTM: DO ihrs = 1, nhrs
     WRITE(ascii_name,"(A,I2.2,A,I2.2,A)") "./test_data/ctm_test_m"&
          ,month,"_",ihrs-1,"h.dat"
     GC_time(ihrs) = ihrs-1
     PRINT *, 'reading from file '//TRIM(ADJUSTL(ascii_name))
     CALL read_gchem_ascii_file  ( &
          TRIM(ADJUSTL(ascii_name)), nlon, nlat, nlev, &
          GC_lon(1:nlon), GC_lat(1:nlat), &
          GC_H2O(1:nlon,1:nlat,1:nlev,ihrs), &
          GC_HCHO(1:nlon,1:nlat,1:nlev,ihrs), &
          GC_NO2(1:nlon,1:nlat,1:nlev,ihrs), &
          GC_O3(1:nlon,1:nlat,1:nlev,ihrs), &
          GC_TEMP(1:nlon,1:nlat,1:nlev,ihrs), &
          GC_PSUR(1:nlon,1:nlat,ihrs))
  END DO READCTM

  ! -----------------------------------------------------------
  ! Create output grid given values for delta_lon and delta_lat
  ! -----------------------------------------------------------
  noutlon = INT(ABS(flon-llon)/dlon+1.0,KIND=2)
  noutlat = INT(ABS(flat-llat)/dlat+1.0,KIND=2)
  ALLOCATE (out_lon(noutlon)); ALLOCATE (out_lat(noutlat))
  ALLOCATE (out_time(nhrs))
  ALLOCATE (out_data_3d(noutlon,noutlat,nhrs),out_lut_data(noutlon,noutlat,nhrs))
  ALLOCATE (out_data_4d(noutlon,noutlat,nlev,nhrs))
  out_lon = [(REAL(ideg,KIND=4)*dlon, ideg = 0, noutlon-1)]+flon
  out_lat = [(REAL(ideg,KIND=4)*dlat, ideg = 0, noutlat-1)]+flat
  out_time = GC_time

  ! ---------------------------------------------------------------
  ! Open HE5 output file and check AMF_SWATH_FILE_ID ( -1 if error)
  ! ---------------------------------------------------------------
  swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(he5file)), he5f_acc_trunc )
  IF ( swath_file_id == -1 ) THEN
     WRITE (*,*) 'ERROR: HE5_SWopen failed!'; STOP 1
  END IF
     
  swath_name = TRIM(ADJUSTL(sensor))//'_Climatology'//'_'&
       //SwathFieldsMonth(month)
  
  ! --------------------------
  ! Create HE5 swath and check 
  ! --------------------------
  swath_id = HE5_SWcreate ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
  IF ( swath_id == -1 ) CALL  he5_error_stop ( 'HE5_SWcreate', '<'//   &
       TRIM(ADJUSTL(swath_name))//'>' )
  
  ! ----------------------------------
  ! Define new dimensions in HE5 swath
  ! ----------------------------------
  he5stat = HE5_SWdefdim ( swath_id, nLonDim, INT(noutlon,8) ) 
  IF ( he5stat /= 0 ) CALL he5_error_stop ( 'HE5_SWdefdim', nLonDim )
  he5stat = HE5_SWdefdim ( swath_id, nLatDim, INT(noutlat,8) )
  IF ( he5stat /= 0 ) CALL he5_error_stop ( 'HE5_SWdefdim', nLatDim )
  he5stat = HE5_SWdefdim ( swath_id, nLevDim, INT(nlev,8) )
  IF ( he5stat /= 0 ) CALL he5_error_stop ( 'HE5_SWdefdim', nLevDim )
  he5stat = HE5_SWdefdim ( swath_id, nHrsDim, INT(nhrs,8) )
  IF ( he5stat /= 0 ) CALL he5_error_stop ( 'HE5_SWdefdim', nHrsDim )

  ! --------------------------
  ! Define Geolocation Fields
  ! --------------------------
  defgeofields: DO idf = 1, ngf
     he5stat = HE5_SWsetfill( swath_id, TRIM(ADJUSTL(GeoFieldNames(idf))), &
          HE5T_NATIVE_FLOAT, fill_value)
     he5stat = HE5_SWdefgfld( swath_id, TRIM(ADJUSTL(GeoFieldNames(idf))), &
          TRIM(ADJUSTL(GeoFieldDims(idf))), " ", HE5T_Native_FLOAT, &
          he5_hdfe_nomerge )
     IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWdefgfld', &
          TRIM(ADJUSTL(GeoFieldNames(idf))) )
     CALL write_field_attributes ( 'geo', idf, swath_id )
  END DO defgeofields
  
  ! ------------------------
  ! Define data fields:
  !   1. Surface Pressure
  !   2. Temperature profile
  !   3... Gas x nmol
  !   4... Gas SD x nmol
  ! ------------------------
  fieldloop: DO idf = 1, ndf
     CALL define_data_fields(INT(nlev,8), idf, swath_id)
     CALL write_field_attributes ( 'data', idf, swath_id )
  END DO fieldloop
  
  ! -------------------------------------------------------------------------------
  ! Detach from and re-attach to created swath (recommended before adding to swath)
  ! -------------------------------------------------------------------------------
  he5stat   = HE5_SWdetach ( swath_id )
  swath_id  = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
  IF ( swath_id == -1 ) CALL  he5_error_stop ( 'HE5_SWattach', TRIM(ADJUSTL(swath_name)) )
  
  ! -----------------------
  ! Write Global Attributes
  ! -----------------------
  CALL write_swath_attributes ( swath_id )
  
  ! -----------------------------------------------
  ! Detach from HE5 swath
  ! -----------------------------------------------
  he5stat = HE5_SWdetach ( swath_id )
  IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWdetach', TRIM(ADJUSTL(swath_name)) )

  ! -----------------------------------------------------------
  ! Attach to created swath (recommended before adding to swath)
  ! -----------------------------------------------------------
  PRINT *, 'Writing to '//TRIM(ADJUSTL(he5file))
  swath_id  = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
  IF ( swath_id == -1 ) CALL  he5_error_stop ( 'HE5_SWattach', TRIM(ADJUSTL(swath_name)) )
  
  ! --------------------------
  ! Write Geolocation Fields
  ! --------------------------
  CALL write_geo_field ( INT(noutlon,8), swath_id, 'Longitude', out_lon(1:noutlon) )
  CALL write_geo_field ( INT(noutlat,8), swath_id, 'Latitude' , out_lat(1:noutlat) )
  CALL write_geo_field ( INT(nhrs,8), swath_id, 'UTC_time', out_time(1:nhrs) )
    
  ! ------------------------------
  ! Write data to climatology file
  ! ------------------------------
  ! --------------------------------------------------------
  ! For each layer we use equation 1 of Ziemke et al.,
  ! Upper tropospheric ozone derived from cloud slicing to
  ! derive layer partial column (DU) from mixing ratios
  ! Partial_column_layer = 0.79 x Delta_pressure_layer x VMR
  ! Pressure units (hPa)
  ! 0.79 constant units (DU hPa-1 ppmv-1)
  ! --------------------------------------------------------
  ! Work put Delta_presure_layer (deltap) using Hybrid grid
  ! definition: Pedge = Ap + [ Bp * Psurf ]
  ! -------------------------------------------------------  
  wrtdatafields: DO idf = 1, ndf
     SELECT CASE (TRIM(ADJUSTL(DataFieldNames(idf))))
     CASE ("SurfacePressure  ")
        print*, 'Processing SurfacePressure...'
        DO ihrs = 1, nhrs
           ! -------------
           ! Interpolation
           ! -------------
           DO ilon = 1, noutlon
              DO ilat = 1, noutlat
                 val = linInterpol(nlon,nlat,REAL(GC_lon,8),REAL(GC_lat,8), &
                      REAL(GC_PSUR(1:nlon,1:nlat,ihrs),8), REAL(out_lon(ilon),8), &
                      REAL(out_lat(ilat),8), status = status)
                 out_data_3d(ilon,ilat,ihrs) = REAL(val,4)
              END DO
           END DO
           ! -------------------------
           ! Write data for given hour
           ! -------------------------
        END DO
        out_data_3d = out_data_3d * ScaleFactors(idf)
        he5_start_3  = (/       0,       0,    0 /)
        he5_stride_3 = (/       1,       1,    1 /)
        he5_edge_3   = (/ noutlon, noutlat, nhrs /)
        he5stat = HE5_SWwrfld ( swath_id, &
             TRIM(ADJUSTL(DataFieldNames(idf))), he5_start_3, he5_stride_3, he5_edge_3, &
             out_data_3d )
     CASE ("HCHO             ")
        print*, 'Processing HCHO...'
        DO ihrs = 1, nhrs
           DO ilev = 1, nlev
              ! -------------
              ! Interpolation
              ! -------------
              DO ilon = 1, noutlon
                 DO ilat = 1, noutlat
                    deltap = (Ap(ilev) + Bp(ilev) * out_data_3d(ilon,ilat,ihrs)) - &
                          (Ap(ilev+1) + Bp(ilev+1) * out_data_3d(ilon,ilat,ihrs))
                    val = linInterpol(nlon,nlat,REAL(GC_lon,8),REAL(GC_lat,8), &
                         REAL(GC_HCHO(1:nlon,1:nlat,ilev,ihrs),8), REAL(out_lon(ilon),8), &
                         REAL(out_lat(ilat),8), status = status)
                    out_data_4d(ilon,ilat,ilev,ihrs) = REAL(val,4)*VMRtoPPM*deltap*0.79
                 END DO
              END DO
              ! -------------------------
              ! Write data for given hour
              ! -------------------------
           END DO
        END DO
        out_data_4d = out_data_4d * ScaleFactors(idf)
        he5_start_4  = (/       0,       0,    0,    0 /)
        he5_stride_4 = (/       1,       1,    1,    1 /)
        he5_edge_4   = (/ noutlon, noutlat, nlev, nhrs /)
        he5stat = HE5_SWwrfld ( swath_id, &
             TRIM(ADJUSTL(DataFieldNames(idf))), he5_start_4, he5_stride_4, he5_edge_4, &
             out_data_4d )
     CASE ("O3               ")
        print*, 'Processing O3...'
        DO ihrs = 1, nhrs
           ! -------------
           ! Interpolation
           ! -------------
           DO ilon = 1, noutlon
              DO ilat = 1, noutlat
                 DO ilev = 1, nlev
                    deltap = (Ap(ilev) + Bp(ilev) * out_data_3d(ilon,ilat,ihrs)) - &
                          (Ap(ilev+1) + Bp(ilev+1) * out_data_3d(ilon,ilat,ihrs))
                    val = linInterpol(nlon,nlat,REAL(GC_lon,8),REAL(GC_lat,8), &
                         REAL(GC_O3(1:nlon,1:nlat,ilev,ihrs),8), REAL(out_lon(ilon),8), &
                         REAL(out_lat(ilat),8), status = status)
                    out_data_4d(ilon,ilat,ilev,ihrs) = REAL(val,4)*VMRtoPPM*deltap*0.79
                    ! Work out corresponding LUT ozone profile
                 END DO              
                 !Fill up the LUT Ozone Profile output array based on Xiong's climatology
                 !Since the dummy climatology I'm using now doesn't have information on
                 !stratospheric ozone for now I create random columns. Based on the scale
                 !ozone colum (in DU) and the specifics of Xiong's ozone climatology I assing
                 !at each scale value (depending on latitude) a given LUT ozone index. The
                 !correspondence of LUT indices is:
                 !1-->L200, 2-->L250, 3-->L300, 4-->L350
                 !5-->M200, 6-->M250, 7-->M300, 8-->M350, 9-->M400, 10-->M450, 11-->M500, 12-->M550
                 !13-->H100, 14-->H150, 15-->H200, 16-->H250, 17-->H300, 18-->H350, 19-->H400
                 !20-->H450, 21-->H500, 22-->H550
                 !For a given latitudinal band just pick the ozone profile with the closest total
                 !column to scale
                 CALL RANDOM_NUMBER(scale)
                 scale = SUM(out_data_4d(1,1,1:nlev,ihrs))+scale*(350.0)+180.0
                 IF ( ABS(out_lat(ilat)) .LE. 30.0 ) THEN
                    out_lut_data(ilon,ilat,ihrs) = L_idx(MINLOC(ABS(L_val-scale),1))
                 ELSE IF ( ABS(out_lat(ilat)) .GT. 30.0 .AND. ABS(out_lat(ilat)) .LE. 55.0 ) THEN
                    out_lut_data(ilon,ilat,ihrs) = M_idx(MINLOC(ABS(M_val-scale),1))
                 ELSE
                    out_lut_data(ilon,ilat,ihrs) = H_idx(MINLOC(ABS(H_val-scale),1))
                 ENDIF                    
              END DO 
           END DO
        END DO
        out_data_4d = out_data_4d * ScaleFactors(idf)
        he5_start_4  = (/       0,       0,    0,    0 /)
        he5_stride_4 = (/       1,       1,    1,    1 /)
        he5_edge_4   = (/ noutlon, noutlat, nlev, nhrs /)
        he5stat = HE5_SWwrfld ( swath_id, &
             TRIM(ADJUSTL(DataFieldNames(idf))), he5_start_4, he5_stride_4, he5_edge_4, &
             out_data_4d )
     CASE ("H2O              ")
        print*, 'Processing H2O...'
        DO ihrs = 1, nhrs
           DO ilev = 1, nlev
              ! -------------
              ! Interpolation
              ! -------------
              DO ilon = 1, noutlon
                 DO ilat = 1, noutlat
                    deltap = (Ap(ilev) + Bp(ilev) * out_data_3d(ilon,ilat,ihrs)) - &
                          (Ap(ilev+1) + Bp(ilev+1) * out_data_3d(ilon,ilat,ihrs))
                    val = linInterpol(nlon,nlat,REAL(GC_lon,8),REAL(GC_lat,8), &
                         REAL(GC_H2O(1:nlon,1:nlat,ilev,ihrs),8), REAL(out_lon(ilon),8), &
                         REAL(out_lat(ilat),8), status = status)
                    out_data_4d(ilon,ilat,ilev,ihrs) = REAL(val,4)*VMRtoPPM*deltap*0.79
                 END DO
              END DO
              ! -------------------------
              ! Write data for given hour
              ! -------------------------
           END DO 
        END DO 
        out_data_4d = out_data_4d * ScaleFactors(idf)
        he5_start_4  = (/       0,       0,    0,    0 /)
        he5_stride_4 = (/       1,       1,    1,    1 /)
        he5_edge_4   = (/ noutlon, noutlat, nlev, nhrs /)
        he5stat = HE5_SWwrfld ( swath_id, &
             TRIM(ADJUSTL(DataFieldNames(idf))), he5_start_4, he5_stride_4, he5_edge_4, &
             out_data_4d )
     CASE ("NO2             ")
        print*, 'Processing NO2...'
        DO ihrs = 1, nhrs
           DO ilev = 1, nlev
              ! -------------
              ! Interpolation
              ! -------------
              DO ilon = 1, noutlon
                 DO ilat = 1, noutlat
                    deltap = (Ap(ilev) + Bp(ilev) * out_data_3d(ilon,ilat,ihrs)) - &
                          (Ap(ilev+1) + Bp(ilev+1) * out_data_3d(ilon,ilat,ihrs))
                    val = linInterpol(nlon,nlat,REAL(GC_lon,8),REAL(GC_lat,8), &
                         REAL(GC_NO2(1:nlon,1:nlat,ilev,ihrs),8), REAL(out_lon(ilon),8), &
                         REAL(out_lat(ilat),8), status = status)
                    out_data_4d(ilon,ilat,ilev,ihrs) = REAL(val,4)*VMRtoPPM*deltap*0.79
                 END DO
              END DO
              ! -------------------------
              ! Write data for given hour
              ! -------------------------
           END DO
        END DO
        out_data_4d = out_data_4d * ScaleFactors(idf)
        he5_start_4  = (/       0,       0,    0,    0 /)
        he5_stride_4 = (/       1,       1,    1,    1 /)
        he5_edge_4   = (/ noutlon, noutlat, nlev, nhrs /)
        he5stat = HE5_SWwrfld ( swath_id, &
             TRIM(ADJUSTL(DataFieldNames(idf))), he5_start_4, he5_stride_4, he5_edge_4, &
             out_data_4d )
        CASE ("LUT Ozone Profile")
           print*, 'Writing LUT ozone profile...'
           he5_start_3 = (/ 0, 0, 0 /)
           he5_stride_3 = (/ 1, 1, 1 /)
           he5_edge_3 = (/ noutlon, noutlat, nhrs /)
           he5stat = HE5_SWwrfld ( swath_id, &
                TRIM(ADJUSTL(DataFieldNames(idf))), he5_start_3, he5_stride_3, he5_edge_3, &
                out_lut_data )
     CASE DEFAULT
        WRITE(*,*) 'Warning!! Field '//TRIM(ADJUSTL(DataFieldNames(idf)))//' not found.' 
     END SELECT
  END DO wrtdatafields

  DEALLOCATE (out_data_3d,out_lut_data,out_data_4d)
  ! ----------------------
  ! Detach from this swath
  ! ----------------------
  he5stat   = HE5_SWdetach ( swath_id )
     
  ! -----------------------------------------------
  ! Close HE5 output file
  ! -----------------------------------------------
  he5stat = HE5_SWclose  ( swath_file_id )
  IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWclose', TRIM(ADJUSTL(he5file)) )

  RETURN
END SUBROUTINE create_he5_file

SUBROUTINE  read_gchem_ascii_file  ( ascii_name, nlon, nlat, nlev, &
                                     longitudes, latitudes,        &
                                     H2O, HCHO, NO2, O3,           &
                                     TEMP, PSUR )

  USE GEOSChem_he5_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER, INTENT (IN) :: nlon, nlat, nlev
  CHARACTER (LEN=*),    INTENT (IN) :: ascii_name
  ! ---------------
  ! Output variables
  ! ---------------
  REAL (KIND=r4), INTENT (OUT), DIMENSION (nlon)           :: longitudes
  REAL (KIND=r4), INTENT (OUT), DIMENSION (nlat)           :: latitudes
  REAL (KIND=r4), INTENT (OUT), DIMENSION (nlon,nlat)      :: PSUR
  REAL (KIND=r4), INTENT (OUT), DIMENSION (nlon,nlat,nlev) :: H2O, &
       HCHO, NO2, O3, TEMP
  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER, PARAMETER :: aunit = 99
  INTEGER            :: ilon, ilat, ios

  OPEN (UNIT=aunit, FILE=TRIM(ADJUSTL(ascii_name)), STATUS='OLD', ACTION='READ', IOSTAT=ios)
  IF ( ios /= 0 ) THEN
     WRITE (*, '(A)') 'ERROR reading from ASCII file '//TRIM(ADJUSTL(ascii_name))
     STOP 1
  END IF

  DO ilon = 1, nlon
     DO ilat = 1, nlat
        READ (UNIT=aunit, FMT=*) longitudes(ilon), latitudes(ilat), PSUR(ilon,ilat)
        READ (UNIT=aunit, FMT=*) H2O(ilon,ilat,1:nlev)
        READ (UNIT=aunit, FMT=*) O3(ilon,ilat,1:nlev)
        READ (UNIT=aunit, FMT=*) HCHO(ilon,ilat,1:nlev)
        READ (UNIT=aunit, FMT=*) NO2(ilon,ilat,1:nlev)
        READ (UNIT=aunit, FMT=*) TEMP(ilon,ilat,1:nlev)
     END DO
  END DO
  CLOSE (aunit)

  RETURN
END SUBROUTINE read_gchem_ascii_file
