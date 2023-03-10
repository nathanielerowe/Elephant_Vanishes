file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/inc/*.h")
add_library(PROfitLib PROconfig.cxx PROspec.cxx ${HEADER_LIST})

# We need this directory, and users of our library will need it too
target_include_directories(PROfitLib PUBLIC ${PROJECT_SOURCE_DIR}/inc)
target_link_libraries(PROfitLib PRIVATE eigen tinyxml2 lbfgspp ${ROOT_LIBRARIES})
target_compile_features(PROfitLib PUBLIC cxx_std_14)

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/inc"
    PREFIX "Header Files"
    FILES ${HEADER_LIST})


