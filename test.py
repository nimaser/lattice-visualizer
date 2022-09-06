# @author nimaser
# @version 0.1.0
# @date Sep 9, 2022

import app
import numpy as np

def main():
    dim = 2
    u, v = (1, 2), (3, 4)
    lattice = app.Lattice(dim, [u, v])

    assert lattice.dimension == 2
    assert np.array_equal(lattice.basis, [u, v])

    assert np.array_equal(lattice.trans_lat_to_std([(1, 0), (0, 1)]), [u, v])
    assert np.array_equal(lattice.trans_std_to_lat([u, v]), [(1, 0), (0, 1)])

    assert np.array_equal(lattice.valid_lattice_points([(2, 1), (1, 2)]), [False, True])

if __name__ == "__main__":
    main()