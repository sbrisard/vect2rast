`voxelize` is a library for the pixelization/voxelization of 2d/3d continuous
objects such as spheres, spheroids, …

`voxelize` is released under a BSD 3-Clause License.

## Installation

### Installation of the C library

First install

- the [Meson build system](https://mesonbuild.com/)
- the [GLib](https://developer.gnome.org/glib/)

Clone the repository of the project

```
git clone https://github.com/sbrisard/voxelize.git
```

and `cd` into the `src` directory

```
cd voxelize/src
meson build
cd build
ninja install
```

(note that the `build/` subdirectory is automatically created; do *not* create
it). To install the library in a custom location, replace the second line with

```
meson -Dprefix=voxelize/installation/directory build
```

To run the tests, stay in the `build/` subdirectory and run the following
command

```
meson test
```

or (more verbose output)

```
./test_voxelize
```

### Installation of the Python wrapper

<!-- Local Variables: -->
<!-- fill-column: 80 -->
<!-- End: -->
