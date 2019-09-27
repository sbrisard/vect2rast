import ctypes
import h5py
import numpy as np

from ctypes import c_byte, c_double, c_size_t

c_byte_p = ctypes.POINTER(c_byte)
c_double_p = ctypes.POINTER(c_double)
c_size_t_p = ctypes.POINTER(c_size_t)

discr = ctypes.cdll.LoadLibrary('libdiscretize.so')
discr.discretize3d.argtypes = [c_double_p,
                               c_double,
                               c_double_p,
                               c_size_t_p,
                               c_byte_p]
discr.discretize3d.restype = None

discr.discretize2d.argtypes = [c_double_p,
                               c_double,
                               c_double_p,
                               c_size_t_p,
                               c_byte_p]
discr.discretize2d.restype = None


def discretize(center, radius, dim, grid):
    size = np.array(grid.shape, dtype=np.intp)
    func = discr.discretize3d if center.shape[0] == 3 else discr.discretize2d
    return func(center.ctypes.data_as(c_double_p),
                radius,
                dim.ctypes.data_as(c_double_p),
                size.ctypes.data_as(c_size_t_p),
                grid.ctypes.data_as(c_byte_p))


if __name__ == '__main__':
    with h5py.File('spheres/spheres-40p100-000.h5', 'r') as f:
        dim = np.array(f['dim'])
        center = np.array(f['center'])
        radius = np.array(f['radius'])

    grid = np.zeros((512, 512, 512), dtype=np.uint8)
    for x_k, r_k in zip(center, radius):
        discretize(x_k, r_k, dim, grid)

    with h5py.File('grid.h5', 'w') as f:
        f['grid'] = grid
