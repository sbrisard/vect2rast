cmake_minimum_required(VERSION 3.15)
project(vect2rast C)

set(CMAKE_C_STANDARD 11)

add_library(vect2rast SHARED src/vect2rast.c src/vect2rast.h src/v2r_ndsphere.c src/v2r_spheroid.c)
target_link_libraries(vect2rast m)

find_package(PkgConfig REQUIRED)
pkg_check_modules(deps REQUIRED IMPORTED_TARGET glib-2.0)

add_executable(test_vect2rast src/test_vect2rast.c src/v2r_test_utils.c src/v2r_test_utils.h src/v2r_test_ndsphere.c src/v2r_test_spheroid.c)
target_link_libraries(test_vect2rast PRIVATE PkgConfig::deps)
target_link_libraries(test_vect2rast PRIVATE vect2rast)
