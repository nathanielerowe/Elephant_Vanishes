add_executable(PROfit PROfit.cxx)

target_compile_features(PROfit PRIVATE cxx_std_14)
#target_link_libraries(PROfit PRIVATE PROfitLib tinyxml2 eigen lbfgspp ${FMT_LIBRARY_PATH}/libfmt.a)
target_link_libraries(PROfit PRIVATE PROfitLib tinyxml2 eigen lbfgspp )


if(MPI_ENABLE)
    add_executable(PROmfa4u PROmfa4u.cxx)
    target_compile_features(PROmfa4u PRIVATE cxx_std_14)
    target_link_libraries(PROmfa4u PRIVATE PROfitLib tinyxml2 eigen lbfgspp )
endif(MPI_ENABLE)




add_executable(PROreader PROreader.cxx ${CMAKE_CURRENT_BINARY_DIR}/PROfit_dict.cxx )
target_include_directories(PROreader PRIVATE ${SBNANAOBJ_INCLUDE_DIR})
target_compile_features(PROreader PRIVATE cxx_std_14)
target_link_directories(PROfit_dict PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(PROreader PRIVATE PROfitLib tinyxml2 eigen lbfgspp PROfit_dict)


