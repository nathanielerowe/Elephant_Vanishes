add_executable(PROfit PROfit.cxx)

target_compile_features(PROfit PRIVATE cxx_std_14)
#target_link_libraries(PROfit PRIVATE PROfitLib tinyxml2 eigen lbfgspp ${FMT_LIBRARY_PATH}/libfmt.a)
target_link_libraries(PROfit PRIVATE PROfitLib tinyxml2 eigen lbfgspp )





file(GLOB SRHeaders ${SBNANAOBJ_INCLUDE_DIR}/sbnanaobj/StandardRecord/*.h)
execute_process(COMMAND rootcling -f ${CMAKE_CURRENT_BINARY_DIR}/PROfit_dict.cxx -c ${SRHeaders} ${CMAKE_SOURCE_DIR}/dict/LinkDef.h -I${SBNANAOBJ_INCLUDE_DIR})
add_library(PROfit_dict ${CMAKE_CURRENT_BINARY_DIR}/PROfit_dict.cxx)
target_include_directories(PROfit_dict PRIVATE ${SBNANAOBJ_INCLUDE_DIR})
target_link_directories(PROfit_dict PRIVATE ${CMAKE_BINARY_DIR})
target_link_libraries(PROfit_dict PRIVATE ${ROOT_LIBRARIES})


add_executable(PROreader PROreader.cxx ${CMAKE_CURRENT_BINARY_DIR}/PROfit_dict.cxx )
target_include_directories(PROreader PRIVATE ${SBNANAOBJ_INCLUDE_DIR})
target_compile_features(PROreader PRIVATE cxx_std_14)
target_link_directories(PROfit_dict PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(PROreader PRIVATE PROfitLib tinyxml2 eigen lbfgspp PROfit_dict)


