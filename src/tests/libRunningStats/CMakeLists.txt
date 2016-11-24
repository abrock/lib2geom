cmake_minimum_required(VERSION 2.8)
project(RunningStats)

if (NOT TARGET RunningStats)
  add_library(RunningStats runningstats.cpp)
endif()
