cmake_minimum_required(VERSION 2.8)


project(SimulationBase)

cmake_minimum_required(VERSION 2.8)

# Requires C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
# Setup of ROOT ANALYSIS : optional.
option(WITH_ANALYSIS_USE "Build example with analysis objects" OFF)

if(WITH_ANALYSIS_USE)
    add_definitions(-DANALYSIS_USE)
    # You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
    #   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
    #   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

    #---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
    find_package(ROOT REQUIRED COMPONENTS RIO)

    #---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
    include(${ROOT_USE_FILE})

    include_directories(${ROOT_INCLUDE_DIRS})
else()
endif()


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                     ${Geant4_INCLUDE_DIR} )

if(WITH_ANALYSIS_USE)
        EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
endif(WITH_ANALYSIS_USE)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#------------------------------------------------------------
add_executable(${PROJECT_NAME} "main.cc" ${sources} ${headers})

if(WITH_ANALYSIS_USE)
        EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LD_FLAGS})
endif(WITH_ANALYSIS_USE)

target_link_libraries(${PROJECT_NAME} ${Geant4_LIBRARIES})
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TestEm14. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
set(${PROJECT_SOURCE_DIR}_SCRIPTS
#    compton.mac
#    conv.mac
    Output.mac
    Geant4Reader.py
    spectraReader.py
    macroWriter.py
    scoring.mac
    spectraScoring.mac
    init.mac
    FlexiSourceMacro.mac
    iridium_source_primary.mac
    vis.mac
    geometry.gdml
  )

foreach(_script ${${PROJECT_SOURCE_DIR}_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

