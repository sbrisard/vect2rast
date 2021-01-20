# add_executable(test_vect2rast
#   test_vect2rast_main.cpp
#   test_vect2rast.cpp
# )

add_executable(test_vect2rast
  test_vect2rast.c
  v2r_test_utils.c
  v2r_test_ndsphere.c
  v2r_test_spheroid.c)

target_include_directories(test_vect2rast PUBLIC ../include)

if (MATH_LIBRARY)
    target_link_libraries(test_vect2rast ${MATH_LIBRARY})
endif ()

#find_package(Catch2 REQUIRED)
target_link_libraries(test_vect2rast vect2rast)

add_test(NAME test_vect2rast COMMAND test_vect2rast)
