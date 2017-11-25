cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(GammaJetAnalyzer)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
list(APPEND CMAKE_PREFIX_PATH $ENV{CLHEP_PATH})
find_package(CLHEP REQUIRED)

#manually set headers path for cmssw, using the submodule
#procedure will be changed completely when cmssw will be built as a library
#and GammaJetAnalyzer will be built properly as a cmssw module :)
set (CMSSW_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/cmssw)

include_directories(${ROOT_INCLUDE_DIRS} ${CMSSW_INCLUDE_DIRS})

#add_executable(GammaJetAnalyzer 
add_library(GammaJetAnalyzer 
	${CMAKE_CURRENT_SOURCE_DIR}/src/GammaJetAnalyzer.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/GammaJetAnalyzer.h
)

#target_link_libraries(GammaJetAnalyzer ${ROOT_LIBRARIES} CLHEP::CLHEP )
