FIND_PACKAGE( ROOT REQUIRED)
find_package(TCLAP REQUIRED)
find_package(RapidXML REQUIRED)


INCLUDE_DIRECTORIES(  ${ROOT_INCLUDE_DIR}  ./inc/ ${TCLAP_INCLUDE_DIR} ${SCT_PLOTS_INCLUDES}  ${RapidXML_INCLUDE_DIR})

add_executable(apix_tel_correlations    src/apix_tel_correlations.cxx)
set_property(TARGET apix_tel_correlations   PROPERTY PROJECT_LABEL "bin_apix_tel_correlations") 


target_link_libraries(apix_tel_correlations    ${PROJECT_NAME})

INSTALL(TARGETS apix_tel_correlations  
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
