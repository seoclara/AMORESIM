#----------------------------------------------------------------------------
# Project Name: AmoreSim
# Description: Simulation library & executable for AmoreSim package
#----------------------------------------------------------------------------
# Setup the project and require CMake version above 3.5
#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(AmoreSim)

#----------------------------------------------------------------------------
# Find ROOT and setup it
#----------------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Find Geant4 and setup it
#----------------------------------------------------------------------------
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
    add_definitions(-DG4VIS_USE)
else()
    find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})
#For backward compatibility
set(CMAKE_CXX_STANDARD 17)


#----------------------------------------------------------------------------
# Add variables for run scripts and macros
#----------------------------------------------------------------------------
set(MUON_PATH "/proj/internal_group/cup_mc/data/spectrum/Muon/muGen/atHandeok_BIGNUM")
set(NEUTRON_PATH "/proj/internal_group/cup_mc/data/spectrum/Neutron/nGen/dist3500")

set(AMORESIM_JOB_NAME "scratch" CACHE STRING "Representative name for this compilation")
if("${AMORESIM_JOB_NAME}" STREQUAL "")
    set(AMORESIM_JOB_NAME "scratch" FORCE)
endif()


set(AMORESIM_DECAY_SOURCE_NAME "U238 Th232 K40 U235 Pb210 Na22 I125 I126 Te121 Te121m Te123m Te125m Te127m H3" CACHE STRING "Target isotopes list for decay simulation (Seperated by space)")
set(AMORESIM_DECAY_RATE "4.9E-18 1.6E-18 1.8E-17 3.1E-17 1E-9 8.5E-9 1.4E-7 6.2E-7 4.2E-7 4.9E-8 6.7E-8 1.4E-7 7.6E-8 1.8E-9" CACHE STRING "Half life list for target isotopes (Seperated by space, Unit in Hz)")
set(AMORESIM_NCORE "100" CACHE STRING "Total number of jobs for simulation of each isotope in the list")
set(AMORESIM_SOURCE_CATEGORY "internal" CACHE STRING "Simulation category for decay simulation (internal or external-pmt)")

#----------------------------------------------------------------------------
# Set advanced properties to variables
#----------------------------------------------------------------------------
mark_as_advanced(AMORESIM_DECAY_SOURCE_NAME)
mark_as_advanced(AMORESIM_DECAY_RATE)
mark_as_advanced(AMORESIM_NCORE)
mark_as_advanced(AMORESIM_SOURCE_CATEGORY)

#----------------------------------------------------------------------------
# Set variables for run scripts and macros
#----------------------------------------------------------------------------
set(SIMOUT_PATH ${PROJECT_BINARY_DIR}/output)

set(AMORESIM_WORK_DIR ${CMAKE_BINARY_DIR})
set(RUN_SCRIPTS
	olaf_jobsubmit.sh
    run_II_muon.sh
    run_II_neut.sh
    run_II_decay.sh
    run_I_neut.sh
    run_I_decay.sh
    run_Pilot_decay.sh
    )
set(SIM_MACROS
    II_muonbckg.mac
    II_neutbckg.mac
    II_dc_internal.mac
    I_neutbckg.mac
    I_dc_internal.mac
    Pilot_dc_external-pmt.mac
    Pilot_dc_internal.mac
    )
set(SESSION_MACROS
    init_vis.mac
    vis.mac
    gui.mac
    )
set(ENVIRONMENT_SCRIPTS environment.sh)

#----------------------------------------------------------------------------
# Configure the scripts and macros
#----------------------------------------------------------------------------
foreach(_SCRIPT ${RUN_SCRIPTS})
    configure_file (
        ${PROJECT_SOURCE_DIR}/run/${_SCRIPT}.in
        ${PROJECT_BINARY_DIR}/run/${_SCRIPT}
        @ONLY
        )
endforeach()
foreach(_SCRIPT ${ENVIRONMENT_SCRIPTS})
    configure_file (
        ${PROJECT_SOURCE_DIR}/${_SCRIPT}.in
        ${PROJECT_BINARY_DIR}/${_SCRIPT}
        @ONLY
        )
endforeach()

foreach(_MACRO ${SIM_MACROS})
    configure_file (
        ${PROJECT_SOURCE_DIR}/mac/${_MACRO}.in
        ${PROJECT_BINARY_DIR}/mac/${_MACRO}
        @ONLY
        )
endforeach()
foreach(_MACRO ${SESSION_MACROS})
    configure_file (
        ${PROJECT_SOURCE_DIR}/${_MACRO}
        ${PROJECT_BINARY_DIR}/${_MACRO}
        @ONLY
        )
endforeach()

#----------------------------------------------------------------------------
# Copy the data files to build directory
#----------------------------------------------------------------------------
file(COPY ${PROJECT_SOURCE_DIR}/mac ${PROJECT_SOURCE_DIR}/run ${PROJECT_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR}
    PATTERN *.in EXCLUDE)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output/log)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output/mac)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output/root)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#----------------------------------------------------------------------------
file(GLOB AmoreSim_EXEC_SOURCE ${PROJECT_SOURCE_DIR}/test/amoresim.cc)
file(GLOB AmoreSim_GitRev_SOURCE ${PROJECT_SOURCE_DIR}/src/AmoreSimGitRevision.cc)
file(GLOB AmoreSim_LIB_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB AmoreSim_LIB_HEADERS ${PROJECT_SOURCE_DIR}/AmoreSim/*.hh)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Add definitions for git commit hash
#----------------------------------------------------------------------------
set_source_files_properties(${AmoreSim_EXEC_SOURCE} PROPERTIES
    COMPILE_DEFINITIONS "AmoreSim_GIT_BRANCH=${AmoreSim_GIT_BRANCH};AmoreSim_GIT_COMMIT_HASH=${AmoreSim_GIT_COMMIT_HASH}" )
set_source_files_properties(${AmoreSim_GitRev_SOURCE} PROPERTIES
    COMPILE_DEFINITIONS "AmoreSim_GIT_BRANCH=${AmoreSim_GIT_BRANCH};AmoreSim_GIT_COMMIT_HASH=${AmoreSim_GIT_COMMIT_HASH}" )

#----------------------------------------------------------------------------
# Add libraries and executables, and link it to the ROOT and Geant4 framework library
#----------------------------------------------------------------------------
add_library(AmoreSimL SHARED ${AmoreSim_LIB_SOURCES})
target_link_libraries(AmoreSimL CupSimL ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
add_executable(amoresim ${AmoreSim_EXEC_SOURCE})
target_link_libraries(amoresim AmoreSimL CupSimL ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Check dependencies for this project and set include directories and libraries
#----------------------------------------------------------------------------
if(DEFINED CupSim)
    include_directories(${CupSim_INCLUDE_DIR})
    target_include_directories(amoresim PUBLIC ${CupSim_INCLUDE_DIR})
    target_include_directories(AmoreSimL PUBLIC ${CupSim_INCLUDE_DIR})
    target_link_libraries(amoresim CupSimL)
else()
    message(SEND_ERROR "This project cannot be built without CupSim. Please check about this.")
endif()

if(DEFINED MCObjs)
    include_directories(${MCObjs_INCLUDE_DIR})
    target_include_directories(amoresim PUBLIC ${MCObjs_INCLUDE_DIR})
    target_include_directories(AmoreSimL PUBLIC ${MCObjs_INCLUDE_DIR})
    target_link_libraries(amoresim MCObjs)
else()
    message(SEND_ERROR "This project cannot be built without MCObjs. Please check about this.")
endif()

#----------------------------------------------------------------------------
# Setup global and local include directories for this project
#----------------------------------------------------------------------------
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${Geant4_INCLUDE_DIRS})
include_directories(${ROOT_INCLUDE_DIRS})
target_include_directories(amoresim PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(amoresim PUBLIC ${Geant4_INCLUDE_DIRS})
target_include_directories(amoresim PUBLIC ${ROOT_INCLUDE_DIRS})
target_include_directories(AmoreSimL PUBLIC ${PROJECT_SOURCE_DIR})
target_include_directories(AmoreSimL PUBLIC ${Geant4_INCLUDE_DIRS})
target_include_directories(AmoreSimL PUBLIC ${ROOT_INCLUDE_DIRS})
if(CMAKE_COMPILER_IS_GNUCXX)
    target_compile_options(AmoreSimL PUBLIC -fdiagnostics-color=always)
    target_compile_options(amoresim PUBLIC -fdiagnostics-color=always)
endif()

#----------------------------------------------------------------------------
# Expose this public includes and library to other subprojects through cache
# variable.
#----------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIR ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)
set(${PROJECT_NAME}_LIBRARY ${PROJECT_BINARY_DIR}/lib${PROJECT_NAME}.so
    CACHE INTERNAL "${PROJECT_NAME}: Path to library" FORCE)
set(${PROJECT_NAME} "Exist!"
    CACHE INTERNAL "${PROJECT_NAME}: Variable for checking dependencies" FORCE)
