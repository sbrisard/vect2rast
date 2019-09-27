import h5py
import numpy as np
import matplotlib.pyplot as plt

from discretize import discretize

if __name__ == '__main__':
    d = 2 # Number of space dimensions
    nincl = 100 # Number of inclusions
    if nincl > 256:
        raise ValueError('Number of inclusions should be smaller than 256')
    f = 0.5 # Volume fraction of inclusions
    nvox = 256 # Number of voxels in each dimension

    ncfg = 1000 # Number of configurations

    # volume of unit-sphere
    v0 = np.pi *(1 if d == 2 else 4/3)

    # Radius, volume and number density of inclusions
    a = 0.5
    v = v0*a**d
    rho = -np.log(f)/v

    # Dimensions of simulation box
    V = nincl/rho
    L = V**(1/d)

    # Size of voxel
    h = L/nvox

    print('Dimensions of simulation box')
    print('----------------------------')
    print('L = {}'.format(L))
    print('L / (2a) = {}'.format(L))


    print('Voxel size')
    print('----------')
    print('h = {}'.format(h))
    print('a/h = {}'.format(a/h))

    radius = np.full((nincl, 1), fill_value=a, dtype=np.float64)
    dim = np.array(d*(L,), dtype=np.float64)


    template = 'bool2d/bool2d-vfrac_{:0.2f}-nincl_{:03d}-{:03d}.h5'
    group_name = 'x'.join(d*(str(nvox),))
    for icfg in range(ncfg):
        print('Generating configuration {}/{}'.format(icfg+1, ncfg))
        center = L*np.random.rand(nincl, d)
        phase = np.zeros(d*(nvox,), dtype=np.uint8)
        for c, r in zip(center, radius):
            discretize(c, r, dim, phase)
        phase = phase.astype(np.bool8).astype(np.int8)

        name = template.format(f, nincl, icfg)

        with h5py.File(name, 'w') as fd:
            group = fd.create_group(group_name)
            fd['dim'] = dim
            fd['center'] = center
            fd['radius'] = radius
            group['phase'] = phase
