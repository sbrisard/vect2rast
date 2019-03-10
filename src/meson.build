project('Voxelize', 'c')
glib = dependency('glib-2.0')
voxlib = library('voxelize', ['voxelize.c', 'voxelize.h'], dependencies: glib)
executable('test_voxelize', 'test_voxelize.c', dependencies: glib, link_with: voxlib)