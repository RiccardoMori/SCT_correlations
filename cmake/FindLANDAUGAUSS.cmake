# print system information
#
# if you are building in-source, this is the same as CMAKE_SOURCE_DIR, otherwise 
# this is the top level directory of your build tree 
MESSAGE( STATUS "Searching for LandauGauss:" )
find_path(LANDAUGAUS_INCLUDE_DIR_dummy landgausFit.h  
          HINTS ${PROJECT_SOURCE_DIR}/extern/LandauGaus/include/ 
          "${LANDAUGAUS}/main/include" 
          "$ENV{LANDAUGAUS}/main/include"  
          )


  find_library(LANDAUGAUS_LIBRARY NAMES LandauGaussFit
    HINTS 
    ${PROJECT_SOURCE_DIR}/extern/LandauGaus/release/lib
    "${LANDAUGAUS}/lib" "$ENV{LANDAUGAUS}/lib")
    
    
set(LANDAUGAUS_LIBRARY ${LANDAUGAUS_LIBRARY} )
set(LANDAUGAUS_INCLUDE_DIR ${LANDAUGAUS_INCLUDE_DIR_dummy} )
set(LANDAUGAUS_DEFINITIONS "-DUSE_LANDAUGAUS" )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(LANDAUGAUSS  DEFAULT_MSG
                                  LANDAUGAUS_LIBRARY LANDAUGAUS_INCLUDE_DIR)
                                  
                                  