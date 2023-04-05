file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/inc/*.h")

    add_library(PROfitLib PROconfig.cxx PROspec.cxx ${HEADER_LIST})
    #add_library(PROfitLib PROconfig.cxx PROspec.cxx  PROcess_weights.cxx ${HEADER_LIST})
if(MPI_ENABLE)
    target_include_directories(PROfitLib PUBLIC ${PROJECT_SOURCE_DIR}/inc ${install_dir}/src/mfa/include/ ${MPI_INCLUDE_PATH} ${DIY_INCLUDE_DIR} ${FMT_INCLUDE_DIR} )
    target_link_libraries(PROfitLib PRIVATE eigen tinyxml2 lbfgspp ${ROOT_LIBRARIES} ${libraries})
    target_compile_features(PROfitLib PUBLIC cxx_std_14)
else ()
    target_include_directories(PROfitLib PUBLIC ${PROJECT_SOURCE_DIR}/inc  )
    target_link_libraries(PROfitLib PRIVATE eigen tinyxml2 lbfgspp ${ROOT_LIBRARIES} ${libraries})
    target_compile_features(PROfitLib PUBLIC cxx_std_14)
endif(MPI_ENABLE) 

# IDEs should put the headers in a nice place
source_group(
    TREE "${PROJECT_SOURCE_DIR}/inc"
    PREFIX "Header Files"
    FILES ${HEADER_LIST})


