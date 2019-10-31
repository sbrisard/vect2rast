`vect2rast` is a library for the vector → raster conversion of simple
 geometrical objects such as spheres, spheroids, …

`vect2rast` is released under a BSD 3-Clause License.

## Installation

### Installation of the C library

First install

- the [Meson build system](https://mesonbuild.com/)
- the [GLib](https://developer.gnome.org/glib/)

Clone the repository of the project

```
git clone https://github.com/sbrisard/v2r.git
```

and `cd` into the `src` directory

```
cd vect2rast/src
meson build
cd build
ninja install
```

(note that the `build/` subdirectory is automatically created; do *not* create
it). To install the library in a custom location, replace the second line with

```
meson -Dprefix=vect2rast/installation/directory build
```

To run the tests, stay in the `build/` subdirectory and run the following
command

```
meson test
```

or (more verbose output)

```
./test_vect2rast
```

### Installation of the Python wrapper

<!-- Local Variables: -->
<!-- fill-column: 80 -->
<!-- End: -->
