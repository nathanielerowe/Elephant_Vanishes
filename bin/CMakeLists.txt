add_executable(PROfit PROfit.cxx)

target_compile_features(PROfit PRIVATE cxx_std_14)
#target_link_libraries(PROfit PRIVATE PROfitLib tinyxml2 eigen lbfgspp ${FMT_LIBRARY_PATH}/libfmt.a)
target_link_libraries(PROfit PRIVATE PROfitLib tinyxml2 eigen lbfgspp )


