set(CMAKE_C_COMPILER "/home/alexander/promotion/afl-2.13b/afl-gcc")
set(CMAKE_CXX_COMPILER "/home/alexander/promotion/afl-2.13b/afl-g++")

ADD_SUBDIRECTORY (2geom)
ADD_SUBDIRECTORY (gtest)
ADD_SUBDIRECTORY (tests)
ADD_SUBDIRECTORY (performance-tests)
ADD_SUBDIRECTORY (toys)

ADD_SUBDIRECTORY (cython)
ADD_SUBDIRECTORY (py2geom)
