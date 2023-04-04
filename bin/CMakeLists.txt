add_executable(PROfit PROfit.cxx)
target_compile_features(PROfit PRIVATE cxx_std_14)
target_link_libraries(PROfit PRIVATE PROfitLib tinyxml2 eigen lbfgspp)

add_executable(PROreader PROreader.cxx)
target_include_directories(PROreader PRIVATE ${PROJECT_SOURCE_DIR}/sbnanaobj)
target_link_directories(PROreader PRIVATE ${CMAKE_SOURCE_DIR})
target_compile_features(PROreader PRIVATE cxx_std_14)
target_link_libraries(PROreader PRIVATE PROfitLib tinyxml2 eigen lbfgspp PROfit_dict)

#ROOT_GENERATE_DICTIONARY(PROreader_dict
#     ${PROJECT_SOURCE_DIR}/sbnanaobj/sbnanaobj/StandardRecord/*.h
#     OPTIONS -I${PROJECT_SOURCE_DIR}/sbnanaobj
#     MODULE PROreader
#     LINKDEF LinkDef.h)

