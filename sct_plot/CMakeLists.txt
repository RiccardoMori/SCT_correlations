FIND_PACKAGE( ROOT REQUIRED)
find_package(RapidXML REQUIRED)




AUX_SOURCE_DIRECTORY( src library_sources )
AUX_SOURCE_DIRECTORY( src/internal library_sources_internals )
AUX_SOURCE_DIRECTORY( src/sct_events library_sources_sct_events )
AUX_SOURCE_DIRECTORY( src/processors library_sources_processors)
SOURCE_GROUP(internals FILES ${library_sources_internals})
SOURCE_GROUP(event FILES ${library_sources_processors})
SOURCE_GROUP(processors FILES ${library_sources_sct_events})
SET(ProjectDict_name SCT_plotsDict.cxx)


INCLUDE_DIRECTORIES(  ${ROOT_INCLUDE_DIR} ${RapidXML_INCLUDE_DIR} )





ADD_LIBRARY( ${PROJECT_NAME} SHARED ${library_sources} ${library_sources_internals} ${library_sources_sct_events} ${library_sources_processors} )



SOURCE_GROUP(Headers  ./inc/*.h)

target_link_libraries( ${PROJECT_NAME} ${ROOT_LIBRARIES} ${ADDITIONAL_LIBRARIES})


INSTALL(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
