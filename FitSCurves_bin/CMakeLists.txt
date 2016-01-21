FIND_PACKAGE( ROOT REQUIRED)
find_package(TCLAP REQUIRED)
find_package(LANDAUGAUSS REQUIRED)
find_package(RapidXML REQUIRED)


INCLUDE_DIRECTORIES(  ${ROOT_INCLUDE_DIR}  ./inc/ ${TCLAP_INCLUDE_DIR} ${SCT_PLOTS_INCLUDES} ${LANDAUGAUS_INCLUDE_DIR}  ${RapidXML_INCLUDE_DIR})

add_executable(FitSCurves    src/FitSCurves.cxx src/getScurvesFromTTree.cc src/scurve_fit_collection.cc)
set_property(TARGET FitSCurves   PROPERTY PROJECT_LABEL "bin_FitSCurves") 


target_link_libraries(FitSCurves    ${PROJECT_NAME} ${LANDAUGAUS_LIBRARY})

INSTALL(TARGETS FitSCurves  
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
