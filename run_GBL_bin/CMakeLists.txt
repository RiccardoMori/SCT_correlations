FIND_PACKAGE( ROOT REQUIRED)
find_package(TCLAP REQUIRED)
find_package(RapidXML REQUIRED)


INCLUDE_DIRECTORIES(  ${ROOT_INCLUDE_DIR}  ./inc/ ${TCLAP_INCLUDE_DIR} ${SCT_PLOTS_INCLUDES} ${RapidXML_INCLUDE_DIR})

add_executable(run_gbl    src/run_GBL.cxx )
set_property(TARGET run_gbl   PROPERTY PROJECT_LABEL "bin_run_gbl") 


target_link_libraries(run_gbl    ${PROJECT_NAME} )

INSTALL(TARGETS run_gbl  
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
