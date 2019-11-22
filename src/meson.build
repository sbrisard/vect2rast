project('vect2rast', 'c')
math = meson.get_compiler('c').find_library('m', required : false)
glib = dependency('glib-2.0')

vect2rast_lib = library('vect2rast',
			['vect2rast.c', 'v2r_ndsphere.c'],
			install: true)
vect2rast_test = executable('test_vect2rast', 'test_vect2rast_sphere.c',
			    dependencies: [math, glib],
			    link_with: vect2rast_lib)
