# Use this variable for tests which provide their own main().
SET(2GEOM_TESTS_SRC
#bezier-utils-test
#lin_alg_test
sbasis-text-test
root-find-test
implicitization-test
#timing-test
#rtree-performance-test
)

find_package(Boost COMPONENTS system filesystem REQUIRED)
add_subdirectory(libRunningStats)

# Use this variable for GTest tests which should have a default main().
SET(2GEOM_GTESTS_SRC
affine-test
angle-test
bezier-test
choose-test
circle-test
convex-hull-test
coord-test
ellipse-test
elliptical-arc-test
intersection-graph-test
line-test
nl-vector-test
path-test
point-test
polynomial-test
rect-test
sbasis-test
half-outline-test
)

file(COPY "half-outline-difficult-paths" DESTINATION ".")
file(COPY "half-outline-simple-paths" DESTINATION ".")
file(COPY "half-outline-known-issues" DESTINATION ".")
file(COPY "half-outline-embroidery-paths" DESTINATION ".")

FOREACH(source ${2GEOM_GTESTS_SRC})
	ADD_EXECUTABLE(${source} ${source}.cpp)
	TARGET_LINK_LIBRARIES(${source} 2geom gtest gtest_main)
	ADD_TEST(NAME ${source} COMMAND ${source})
ENDFOREACH()

FOREACH(source ${2GEOM_TESTS_SRC})
    ADD_EXECUTABLE(${source} ${source}.cpp)
    TARGET_LINK_LIBRARIES(${source} 2geom ${LINK_GSL} ${GTK2_LINK_FLAGS} gtest)
    ADD_TEST(NAME ${source} COMMAND ${source})
ENDFOREACH(source)

target_link_libraries(half-outline-test RunningStats)
target_link_libraries(half-outline-test ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(half-outline-test geom-pathstroke 2geom)


