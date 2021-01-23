# add_executable(test_vect2rast
#   test_vect2rast_main.cpp
#   test_vect2rast.cpp
# )

add_executable(test_vect2rast
        test_vect2rast.cpp
        v2r_test_utils.cpp
        test_hypersphere.cpp
        v2r_test_spheroid.cpp)

target_include_directories(test_vect2rast PUBLIC ../include)

if (MATH_LIBRARY)
    target_link_libraries(test_vect2rast ${MATH_LIBRARY})
endif ()
add_compile_definitions(_USE_MATH_DEFINES)

#find_package(Catch2 REQUIRED)
target_link_libraries(test_vect2rast vect2rast)

add_test(NAME test_vect2rast COMMAND test_vect2rast)
