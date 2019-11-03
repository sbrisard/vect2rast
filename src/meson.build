project('vect2rast', 'c')
math = meson.get_compiler('c').find_library('m', required : false)
glib = dependency('glib-2.0')

v2r_lib = library('v2r', ['v2r.c', 'v2r.h'], dependencies: glib, install: true)
v2r_test = executable('test_v2r', 'test_v2r.c', dependencies: glib, link_with: v2r_lib)

vect2rast_lib = library('vect2rast',
			['vect2rast.c', 'v2r_ndsphere.c'],
			install: true)
vect2rast_test = executable('test_vect2rast', 'test_vect2rast_sphere.c',
			    dependencies: [math, glib],
			    link_with: vect2rast_lib)

test('tests of vect2rast (legacy)', v2r_test)
