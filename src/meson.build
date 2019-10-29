project('v2r', 'c')
glib = dependency('glib-2.0')

v2r_lib = library('v2r', ['v2r.c', 'v2r.h'], dependencies: glib, install: true)
v2r_test = executable('test_v2r', 'test_v2r.c', dependencies: glib, link_with: v2r_lib)

v2r2_lib = library('v2r2', ['v2r2.c', 'v2r2.h'], install: true)
v2r2_exe = executable('v2r2', ['v2r2.c', 'v2r2.h'])
v2r2_test = executable('test_v2r_sphere', 'test_v2r_sphere.c', dependencies: glib, link_with: v2r2_lib)

test('tests of v2r', v2r_test)
