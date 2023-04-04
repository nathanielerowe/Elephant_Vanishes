#Lets aim to be version 3.19 upwards
cmake_minimum_required (VERSION 3.19)

project (
    PROfit
    VERSION 0.1
    DESCRIPTION "A PROfessional PROfitable Fitter"
    LANGUAGES CXX)

include(FetchContent)
include(ExternalProject)


#Allow for nice support folders in IDEs like VScode
set_property(GLOBAL PROPERTY USE_FOLDERS ON)


find_package(Doxygen)
if(Doxygen_FOUND)
    add_subdirectory(docs)
else()
    message(STATUS "Doxygen not found, not building docs")
endif()


##### Some Dependancies ######

FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git 
    GIT_TAG        3.4
    )
FetchContent_MakeAvailable(eigen)

FetchContent_Declare(
    lbfgspp
    GIT_REPOSITORY https://github.com/yixuan/LBFGSpp.git 
    GIT_TAG        v0.2.0
    )
FetchContent_MakeAvailable(lbfgspp)

FetchContent_Declare(
    tinyxml2
    GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git 
    GIT_TAG        9.0.0
    )
FetchContent_MakeAvailable(tinyxml2)

#./configure --prefix=$MPICH_DIR/mpich-install --disable-f08 --disable-collalgo-tests
#ExternalProject_Add(mpich
#    DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
#    URL https://www.mpich.org/static/downloads/4.0.3/mpich-4.0.3.tar.gz
#    UPDATE_COMMAND ""
#    BUILD_IN_SOURCE 1
#    CONFIGURE_COMMAND ./configure  --disable-f08 --disable-collalgo-tests
#    INSTALL_COMMAND ""
#    )

ExternalProject_Add(mfa
    GIT_REPOSITORY https://markrosslonergan@github.com/markrosslonergan/mfa4profit.git
    GIT_TAG 0e0cc41e027603668a9ca67b9a9195b3d6be50a3
    USES_TERMINAL_DOWNLOAD ON 
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    CMAKE_FLAGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/mfa_build/
    INSTALL_COMMAND ""
    )
ExternalProject_Get_Property(mfa install_dir)

message("MFA installed to: ${install_dir}")

set(MFA_INCLUDE_DIR ${install_dir}/src/mfa/include)
set(DIY_INCLUDE_DIR ${MFA_INCLUDE_DIR}/diy/include)
set(FMT_LIBRARY_PATH ${MFA_INCLUDE_DIR}/fmt)
set(FMT_INCLUDE_DIR ${FMT_LIBRARY_PATH}/include)

message("MFA include dir: ${MFA_INCLUDE_DIR}")
message("DIY include dir: ${DIY_INCLUDE_DIR}")

if                          (MFA_INCLUDE_DIR)
    include_directories     (SYSTEM ${MFA_INCLUDE_DIR} ${DIY_INCLUDE_DIR} ${FMT_INCLUDE_DIR})
    message                 ("MFA_INCLUDE_DIR = " ${MFA_INCLUDE_DIR})
    add_definitions         (-DMFA_NO_WEIGHTS)
    add_definitions         (-DMFA_SERIAL)
    message                 ("Weights are disabled for MFA")
endif  ()


#FetchContent_Declare( got thorugh MFA
#    diy
#    GIT_REPOSITORY https://github.com/diatomic/diy.git 
#    GIT_TAG        ef460f828c36e8f533ccec1f551024f4da896165
#    )
#FetchContent_MakeAvailable(diy)


find_package(Boost REQUIRED)

###### ROOT ROOT ##############33

FIND_PACKAGE(ROOT COMPONENTS MathCore MathMore)

if(ROOT_FOUND)
    message("Found ROOT using CMAKE FIND_PACKAGE")
    include(${ROOT_USE_FILE})
ELSE()
    message("ROOT not compiled using CMAKE, use root-config (BAD!)")
    FIND_PROGRAM(ROOT_CONFIG root-config)
    EXEC_PROGRAM(${ROOT_CONFIG} ARGS --cflags --glibs OUTPUT_VARIABLE ROOT_CXX_FLAGS)
    message("ROOT cxx flags set to: ${ROOT_CXX_FLAGS}")
ENDIF()

find_package(MPI REQUIRED)
set(libraries ${libraries} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
message("MPI libraries: ${libraries}")
message("MPI include path: ${MPI_INCLUDE_PATH}")



# The compiled library code lives inside here
add_subdirectory(src)

# The executable code will be placed here
add_subdirectory(bin)


