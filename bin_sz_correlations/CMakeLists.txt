FIND_PACKAGE( ROOT REQUIRED)
find_package(TCLAP REQUIRED)
find_package(RapidXML REQUIRED)


INCLUDE_DIRECTORIES(  ${ROOT_INCLUDE_DIR}  ./inc/ ${TCLAP_INCLUDE_DIR} ${SCT_PLOTS_INCLUDES}  ${RapidXML_INCLUDE_DIR})

add_executable(sz_correlations    src/run_sz_correlations.cxx)
set_property(TARGET sz_correlations   PROPERTY PROJECT_LABEL "bin_sz_correlations") 


target_link_libraries(sz_correlations    ${PROJECT_NAME})

INSTALL(TARGETS sz_correlations  
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
