FIND_PACKAGE( ROOT REQUIRED)
find_package(TCLAP REQUIRED)


INCLUDE_DIRECTORIES(  ${ROOT_INCLUDE_DIR} ${TCLAP_INCLUDE_DIR} )

add_executable(main    main.cxx  )
set_property(TARGET main   PROPERTY PROJECT_LABEL "bin_main") 



# ${ADDITIONAL_LIBRARIES} is only set if e.g. the native reader processor is built (EUTelescope/LCIO)
target_link_libraries(main    ${PROJECT_NAME} )

INSTALL(TARGETS main  
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
