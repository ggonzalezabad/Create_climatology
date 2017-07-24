MODULE GEOSChem_he5_module

  USE ISO_C_BINDING

  IMPLICIT NONE
  INCLUDE 'hdfeos5.inc'


  ! --------------
  ! Precision KIND
  ! --------------
  INTEGER, PARAMETER :: r8 = KIND(1.0D0)
  INTEGER, PARAMETER :: r4 = KIND(1.0  )

  ! ---------------------
  ! Data field fill value
  ! ---------------------
  REAL(KIND=r4), PARAMETER :: fill_value = -1.0E30
  INTEGER, PARAMETER :: fill_value_i4 = -32767

  ! --------------------------
  ! Volume mixing ratio to PPM
  ! --------------------------
  REAL*4, PARAMETER :: VMRtoPPM = 1E6

  ! ----------------------
  ! Hybrid grid definition
  ! ----------------------
  ! ---------------------------------
  ! GMAO GEOS-5 hybrid grid Ap and Bp
  ! ---------------------------------
  REAL(KIND=r4), DIMENSION(48), PARAMETER :: Ap=(/0.000000E+00, 4.804826E-02, 6.593752E+00, 1.313480E+01, &
       1.961311E+01, 2.609201E+01, 3.257081E+01, 3.898201E+01, &
       4.533901E+01, 5.169611E+01, 5.805321E+01, 6.436264E+01, &
       7.062198E+01, 7.883422E+01, 8.909992E+01, 9.936521E+01, &
       1.091817E+02, 1.189586E+02, 1.286959E+02, 1.429100E+02, &
       1.562600E+02, 1.696090E+02, 1.816190E+02, 1.930970E+02, &
       2.032590E+02, 2.121500E+02, 2.187760E+02, 2.238980E+02, &
       2.243630E+02, 2.168650E+02, 2.011920E+02, 1.769300E+02, &
       1.503930E+02, 1.278370E+02, 1.086630E+02, 9.236572E+01, &
       7.851231E+01, 5.638791E+01, 4.017541E+01, 2.836781E+01, &
       1.979160E+01, 9.292942E+00, 4.076571E+00, 1.650790E+00, &
       6.167791E-01, 2.113490E-01, 6.600001E-02, 1.000000E-02/)
  REAL(KIND=r4), DIMENSION(48), PARAMETER :: Bp=(/1.000000E+00, 9.849520E-01, 9.634060E-01, 9.418650E-01, &
       9.203870E-01, 8.989080E-01, 8.774290E-01, 8.560180E-01, &
       8.346609E-01, 8.133039E-01, 7.919469E-01, 7.706375E-01, &
       7.493782E-01, 7.211660E-01, 6.858999E-01, 6.506349E-01, &
       6.158184E-01, 5.810415E-01, 5.463042E-01, 4.945902E-01, &
       4.437402E-01, 3.928911E-01, 3.433811E-01, 2.944031E-01, &
       2.467411E-01, 2.003501E-01, 1.562241E-01, 1.136021E-01, &
       6.372006E-02, 2.801004E-02, 6.960025E-03, 8.175413E-09, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
       0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00/)

  ! ------------
  ! HE5 Routines
  ! ------------
  INTEGER (KIND = 4), EXTERNAL ::                                               &
       he5_swattach,   he5_swclose,  he5_swcreate, he5_swdefdfld, he5_swdefgfld, &
       he5_swdefdim,   he5_swdetach, he5_swopen, he5_swwrfld,   he5_swwrlattr, &
       he5_swdefcomch, he5_swwrattr, he5_swsetfill

  ! ----------------------------------
  ! Start, Stride, and Endge Variables
  ! ----------------------------------
  INTEGER(C_LONG)               :: he5_start_1, he5_stride_1, he5_edge_1
  INTEGER(C_LONG), DIMENSION(2) :: he5_start_2, he5_stride_2, he5_edge_2
  INTEGER(C_LONG), DIMENSION(3) :: he5_start_3, he5_stride_3, he5_edge_3
  INTEGER(C_LONG), DIMENSION(4) :: he5_start_4, he5_stride_4, he5_edge_4
  
  ! -----------------------------------------------------------------------
  ! Parameters and variables associated with field compression and chunking
  ! -----------------------------------------------------------------------
  INTEGER(C_LONG), DIMENSION (4) :: cchunk_dim
  INTEGER                        :: comp_par, comp_type, ncchunkdim

  ! ----------------------------------------
  ! General Swath Fields, one for each month
  ! ----------------------------------------
  INTEGER, PARAMETER :: nmonth = 12

  CHARACTER (LEN=3), DIMENSION (nmonth), PARAMETER :: SwathFieldsMonth = (/    &
       "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" /)

  ! ------------------------------
  ! Names of Dimensions in Swaths
  ! ------------------------------
  CHARACTER (LEN= 4), PARAMETER :: nLonDim = 'nLon'
  CHARACTER (LEN= 4), PARAMETER :: nLatDim = 'nLat'
  CHARACTER (LEN= 4), PARAMETER :: nLevDim = 'nLev'
  CHARACTER (LEN= 4), PARAMETER :: nHrsDim = 'nHrs'
  CHARACTER (LEN= 9), PARAMETER :: n2Dim   = 'nLon,nLat'
  CHARACTER (LEN=14), PARAMETER :: n3Dim   = 'nLon,nLat,nHrs'
  CHARACTER (LEN=19), PARAMETER :: n4Dim   = 'nLon,nLat,nLev,nHrs'

  ! ----------------------------------------------
  ! Number of data and geolocation fields in swath
  ! ----------------------------------------------
  INTEGER, PARAMETER :: ndf = 6, ngf = 3

  ! ------------------------------
  ! Names of Data Fields in Swaths
  ! ------------------------------
  CHARACTER (LEN=13), DIMENSION (ngf), PARAMETER :: GeoFieldNames = (/ "Longitude", "Latitude ", "UTC_time "/)
  CHARACTER (LEN=4), DIMENSION (ngf), PARAMETER :: GeoFieldDims = (/ nLonDim, nLatDim, nHrsDim /)

  ! ---------------------------------------
  ! Titles of Geolocation Fields in Swaths
  ! --------------------------------------
  CHARACTER (LEN=40), DIMENSION (ngf), PARAMETER :: GeoFieldTitles = (/ &
       "Geodetic Longitude at Center of Grid Box", &
       "Geodetic Latitude at Center of Grid Box ", &
       "Local time                              " /)

  ! -------------------------------------
  ! Units of Geolocation Fields in Swaths
  ! -------------------------------------
  CHARACTER (LEN=8), DIMENSION (ngf), PARAMETER :: GeoFieldUnits = (/ &
       "deg     ",  &
       "deg     ",  &
       "hours   " /) 

  ! ------------------------------
  ! Names of Data Fields in Swaths
  ! ------------------------------
  CHARACTER (LEN=17), DIMENSION (ndf), PARAMETER :: DataFieldNames =(/ &
       "SurfacePressure  ", &
       "O3               ", &
       "HCHO             ", &
       "H2O              ", &
       "NO2              ", &
       "LUT Ozone Profile" /)
       

  ! ------------------------------
  ! Units of Data Fields in Swaths
  ! ------------------------------
  CHARACTER (LEN=8), DIMENSION (ndf), PARAMETER :: DataFieldUnits = (/ &
       "hPa     ", &
       "DU      ", &
       "DU      ", &
       "DU      ", &
       "DU      ", &
       "Unitless" /) 
 
  ! --------------------------------------
  ! Scale factors of Data Fields in Swaths
  ! --------------------------------------
  REAL (KIND=r4), DIMENSION (ndf), PARAMETER :: ScaleFactors = (/ &
       1.0E-00_r4,  &
       1.0E-00_r4,  &
       1.0E-00_r4,  &
       1.0E-00_r4,  &
       1.0E-00_r4,  &
       1.0E-00_r4 /)
 
  ! ------------------------------
  ! Titles of Data Fields in Swaths
  ! ------------------------------
  CHARACTER (LEN=25), DIMENSION (ndf), PARAMETER :: DataFieldTitles = (/ &
       "Surface Pressure         ", &
       "O3 Layer Partial Column  ", &
       "HCHO Layer Partial Column", &
       "H2O Layer Partial Column ", &
       "NO2 Layer Partial Column ", &
       "LUT Ozone Profile        "  /)

  ! ------------------------------
  ! Names of Data Provider
  ! ------------------------------
  CHARACTER (LEN=3), PARAMETER :: DataProvider = "TBD"

  ! -------------------------------
  ! Version of GEOS-Chem Simulation
  ! -------------------------------
  CHARACTER (LEN=3), PARAMETER :: DataDescription = "TBD"

CONTAINS

  SUBROUTINE define_data_fields ( nlev, ifield, swath_id )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,               INTENT (IN) :: ifield, swath_id
    INTEGER (KIND=C_LONG), INTENT (IN) :: nlev

    ! --------------
    ! Local variable
    ! --------------
    INTEGER :: he5stat

    ! --------------------------------------------
    ! Define array value data fields in HE5 swath.
    ! SurfacePressure is the only 3D variable.
    ! Everything else is 4D. 
    ! For now I'm not using compression.
    ! --------------------------------------------
    IF ( TRIM(ADJUSTL(DataFieldNames(ifield))) == 'SurfacePressure' ) THEN
       he5stat = HE5_SWsetfill( swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), &
            HE5T_NATIVE_FLOAT, fill_value)
       he5stat = HE5_SWdefdfld ( swath_id, &
            TRIM(ADJUSTL(DataFieldNames(ifield))),  n3Dim, " ", HE5T_Native_FLOAT, he5_hdfe_nomerge )
    ELSE IF ( TRIM(ADJUSTL(DataFieldNames(ifield))) == 'LUT Ozone Profile' ) THEN 
       he5stat = HE5_SWsetfill( swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), &
            HE5T_NATIVE_INT, fill_value_i4)
       he5stat = HE5_SWdefdfld ( swath_id, &
            TRIM(ADJUSTL(DataFieldNames(ifield))),  n3Dim, " ", HE5T_Native_INT, he5_hdfe_nomerge )
    ELSE
       he5stat = HE5_SWsetfill( swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), &
            HE5T_NATIVE_FLOAT, fill_value)
       comp_par   = 9 ; comp_type = HE5_HDFE_COMP_SHUF_DEFLATE
       ncchunkdim = 4 ; cchunk_dim(1:ncchunkdim) = (/ INT(6,8), INT(6,8), INT(nlev,8), INT(6,8)/)
       he5stat = HE5_SWdefcomch ( &
            swath_id, comp_type, comp_par, ncchunkdim, cchunk_dim(1:ncchunkdim) )
       IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWdefcomch', TRIM(ADJUSTL(DataFieldNames(ifield))) )
       he5stat = HE5_SWdefdfld ( swath_id, &
            TRIM(ADJUSTL(DataFieldNames(ifield))),  n4Dim, " ", HE5T_Native_FLOAT, he5_hdfe_nomerge )
    END IF
    IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWdefdfld', TRIM(ADJUSTL(DataFieldNames(ifield))) )

    RETURN
  END SUBROUTINE define_data_fields

  SUBROUTINE write_data_fields ( nlon, nlat, nlev, nhrs, ifield, swath_id, gcdata)

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,                                            INTENT (IN) :: ifield, swath_id
    INTEGER (KIND=C_LONG),                              INTENT (IN) :: nlon, nlat, nlev, nhrs
    REAL    (KIND=r4), DIMENSION (nlon,nlat,nlev,nhrs), INTENT (IN) :: gcdata

    ! --------------
    ! Local variable
    ! --------------
    INTEGER :: he5stat

    ! ---------------------------------------------------------------------
    ! Write array value data fields in HE5 swath. SurfacePressure is the
    ! only 3D variable, everything else is 4D/
    ! ---------------------------------------------------------------------
    IF ( TRIM(ADJUSTL(DataFieldNames(ifield))) == 'SurfacePressure' ) THEN
       he5_start_3  = (/    0,    0,    0 /)
       he5_stride_3 = (/    1,    1,    1 /)
       he5_edge_3   = (/ nlon, nlat, nhrs /)

       he5stat = HE5_SWwrfld ( swath_id, &
            TRIM(ADJUSTL(DataFieldNames(ifield))), he5_start_3, he5_stride_3, he5_edge_3, &
            gcdata(1:nlon,1:nlat,1,1:nhrs) )
    ELSE

       he5_start_4  = (/    0,    0,    0,    0 /)
       he5_stride_4 = (/    1,    1,    1,    1 /)
       he5_edge_4   = (/ nlon, nlat, nlev, nhrs /)

       he5stat = HE5_SWwrfld ( swath_id, &
            TRIM(ADJUSTL(DataFieldNames(ifield))), he5_start_4, he5_stride_4, he5_edge_4, &
            gcdata(1:nlon,1:nlat,1:nlev,1:nhrs) )
    END IF
    IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWwrfld', TRIM(ADJUSTL(DataFieldNames(ifield))) )


    RETURN
  END SUBROUTINE write_data_fields

  SUBROUTINE write_geo_field ( ndim, swath_id, geofield, gfdata )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,                             INTENT (IN) :: swath_id
    INTEGER (KIND=C_LONG),               INTENT (IN) :: ndim
    REAL (KIND=r4), DIMENSION (ndim), INTENT (IN) :: gfdata
    CHARACTER (LEN=*),                   INTENT (IN) :: geofield

    ! --------------
    ! Local variable
    ! --------------
    INTEGER :: he5stat

    he5_start_1  = 0 ; he5_stride_1 = 1 ; he5_edge_1   = ndim

    he5stat = HE5_SWwrfld ( swath_id, &
            TRIM(ADJUSTL(geofield)), he5_start_1, he5_stride_1, he5_edge_1, gfdata(1:ndim) )
    IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWwrfld', TRIM(ADJUSTL(geofield)) )

    RETURN
  END SUBROUTINE write_geo_field

  SUBROUTINE write_field_attributes ( geodat, ifield, swath_id )
  
    IMPLICIT NONE
  
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,           INTENT (IN) :: ifield, swath_id
    CHARACTER (LEN=*), INTENT (IN) :: geodat
    ! --------------
    ! Local variable
    ! --------------
    INTEGER(C_LONG) :: n_attr, n_one = 1
    INTEGER         :: he5stat
  
  
    ! -------------------------------------------------
    ! Attributes for data fields but geolocation fields
    ! -------------------------------------------------
    IF ( INDEX ( TRIM(ADJUSTL(geodat)), 'geo') > 0 ) THEN

       n_attr = LEN_TRIM(ADJUSTL(GeoFieldTitles(ifield)))
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(GeoFieldNames(ifield))), "Title", &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(GeoFieldTitles(ifield))) )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'Title')

       he5stat = HE5_SWwrlattr (                                            &
            swath_id, TRIM(ADJUSTL(GeoFieldNames(ifield))), "ScaleFactor",  &
            HE5T_NATIVE_FLOAT, n_one, 1.0E+00                                 )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'ScaleFactor')
  
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(GeoFieldNames(ifield))), "MissingValue", &
            HE5T_NATIVE_FLOAT, n_one, fill_value                              )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'MissingValue')
  
       n_attr = LEN_TRIM(ADJUSTL(GeoFieldUnits(ifield)))
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(GeoFieldNames(ifield))), "Units",         &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(GeoFieldUnits(ifield))) )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'Units')

    ELSE
       n_attr = LEN_TRIM(ADJUSTL(DataFieldTitles(ifield)))
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "Title", &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(DataFieldTitles(ifield))) )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'Title')

       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "ScaleFactor",  &
            HE5T_NATIVE_FLOAT, n_one, ScaleFactors(ifield)                    )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'ScaleFactor')

       n_attr = LEN_TRIM(ADJUSTL(DataFieldUnits(ifield)))
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "Units", &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(DataFieldUnits(ifield))) )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'Units')
  
    END IF  

    RETURN
  END SUBROUTINE write_field_attributes


  SUBROUTINE write_swath_attributes ( swath_id )
  
    IMPLICIT NONE
  
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER, INTENT (IN) :: swath_id

    ! --------------
    ! Local variable
    ! --------------
    INTEGER(C_LONG) :: n_attr
    INTEGER         :: he5stat
    CHARACTER(LEN=8):: date
  
  
    n_attr = LEN_TRIM(ADJUSTL("Gonzalo Gonzalez Abad"))
    he5stat = HE5_SWwrattr (                                             &
         swath_id, "Author", &
         HE5T_NATIVE_CHAR, n_attr, "Gonzalo Gonzalez Abad" )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrgattr', 'Author')

    n_attr = LEN_TRIM(ADJUSTL("Harvard-Smithsonian Center for Astrophysics"))
    he5stat = HE5_SWwrattr (                                             &
         swath_id, "AuthorAffiliation", &
         HE5T_NATIVE_CHAR, n_attr, "Harvard-Smithsonian Center for Astrophysics" )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrgattr', 'AuthorAffiliation')
  
    n_attr = LEN_TRIM(ADJUSTL("ggonzale@cfa.harvard.edu"))
    he5stat = HE5_SWwrattr (                                             &
         swath_id, "AuthorContact", &
         HE5T_NATIVE_CHAR, n_attr, "ggonzale@cfa.harvard.edu" )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrgattr', 'AuthorContact')
  
    CALL date_and_time(DATE=date)
    n_attr = LEN_TRIM(ADJUSTL(date))
    he5stat = HE5_SWwrattr (                                             &
         swath_id, "DateOfProduction", &
         HE5T_NATIVE_CHAR, n_attr, date )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrgattr', 'DateOfProduction')

    n_attr = LEN_TRIM(ADJUSTL(DataDescription))
    he5stat = HE5_SWwrattr (                                                &
         swath_id, "DataDescription", &
         HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(DataDescription))   )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrgattr', 'DataDescription')

    n_attr = LEN_TRIM(ADJUSTL(DataProvider))
    he5stat = HE5_SWwrattr (                                             &
         swath_id, "DataProvider", &
         HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(DataProvider)) )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrgattr', 'DataProvider')

    RETURN
  END SUBROUTINE write_swath_attributes


  SUBROUTINE he5_error_stop ( he5routine, fieldname )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT (IN) :: he5routine, fieldname

    WRITE (*,'(A)') 'ERROR: '//TRIM(ADJUSTL(he5routine))//' failed for '//TRIM(ADJUSTL(fieldname))
    STOP 1

    RETURN
  END SUBROUTINE he5_error_stop

END MODULE GEOSChem_he5_module
