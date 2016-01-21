FIND_PACKAGE( ROOT REQUIRED)
find_package(TCLAP REQUIRED)
find_package(RapidXML REQUIRED)


INCLUDE_DIRECTORIES(  ${ROOT_INCLUDE_DIR}  ./inc/ ${TCLAP_INCLUDE_DIR} ${SCT_PLOTS_INCLUDES} ${RapidXML_INCLUDE_DIR})

add_executable(ProcessCollection    src/ProcessCollection.cxx )
set_property(TARGET ProcessCollection   PROPERTY PROJECT_LABEL "bin_ProcessCollection") 


target_link_libraries(ProcessCollection    ${PROJECT_NAME} )

INSTALL(TARGETS ProcessCollection  
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
