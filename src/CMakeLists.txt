file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/inc/*.h")
add_library(PROfitLib PROconfig.cxx PROspec.cxx PROmfa.cxx ${HEADER_LIST})
#add_library(PROfitLib PROconfig.cxx PROspec.cxx PROcovariancegen.cxx PROmfa.cxx ${HEADER_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(PROfitLib PUBLIC ${PROJECT_SOURCE_DIR}/inc ${install_dir}/src/mfa/include/ ${MPI_INCLUDE_PATH} ${DIY_INCLUDE_DIR} ${FMT_INCLUDE_DIR} )
target_link_libraries(PROfitLib PRIVATE eigen tinyxml2 lbfgspp ${ROOT_LIBRARIES} ${libraries})
target_compile_features(PROfitLib PUBLIC cxx_std_14)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/inc"
    PREFIX "Header Files"
    FILES ${HEADER_LIST})


