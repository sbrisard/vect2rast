project('voxelize', 'c')
glib = dependency('glib-2.0')

voxelize_lib = library('voxelize', ['voxelize.c', 'voxelize.h'], dependencies: glib, install: true)
voxelize_test = executable('test_voxelize', 'test_voxelize.c', dependencies: glib, link_with: voxelize_lib)

test('tests of voxelize', voxelize_test)
