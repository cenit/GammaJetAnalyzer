cmake_minimum_required(VERSION 2.8)
project(GammaJetAnalyzer)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED)

include_directories(${ROOT_INCLUDE_DIRS})


add_executable(GammaJetAnalyzer 
	${CMAKE_CURRENT_SOURCE_DIR}/src/GammaJetAnalyzer.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/GammaJetAnalyzer.h
)

target_link_libraries(GammaJetAnalyzer ${ROOT_LIBRARIES})
