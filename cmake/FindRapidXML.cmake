# print system information
#
# if you are building in-source, this is the same as CMAKE_SOURCE_DIR, otherwise 
# this is the top level directory of your build tree 
MESSAGE( STATUS "Searching for Rapid XML:" )
find_path(RapidXML_INCLUDE_DIR rapidxml.hpp  HINTS ${PROJECT_SOURCE_DIR}/extern/rapidxml-1.13/)





set(RapidXML_LIBRARIES "" )
set(RapidXML_INCLUDE_DIR ${RapidXML_INCLUDE_DIR} )
set(RapidXML_DEFINITIONS "-DUSE_RapidXML" )



include(FindPackageHandleStandardArgs)


find_package_handle_standard_args(RapidXML DEFAULT_MSG
                                   RapidXML_INCLUDE_DIR)
