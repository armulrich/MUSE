### Below is the original mgrid_301.py file. Note running this in this environment will not work as the packages are outdated and fixes need to be implemented
### Runable version of this is mgrid_301.py (which includes patches to adress those issues)

# from re import M
import numpy as np
import matplotlib.pyplot as plt
from magtense import magtense
from coilpy import *

def updateTiles2FICUS(tiles,fname):
    nmag = tiles.n
    # MUSEdata = np.zeros((nmag, 42))
    FICUSdata = np.zeros((nmag,15))
    # header0 = "centroid(3), dimension-HLW(3), n1-H(3), n2-L(3), n3-W(3), magnetic-orientation(3), 8-corner-points(24)"
    header1 = "X [m], Y[m], Z[m], n1x, n1y, n1z, n2x, n2y, n2z, H [m], L [m], M [A m], mx, my, mz"

    rMatrix = np.zeros((nmag,3,3))
    for i in range(nmag):
        rMatrix[i] = rotation_matrix(tiles.rot[i,0],tiles.rot[i,1],tiles.rot[i,2], xyz=True).T
        
    # MUSEdata[:,0:3] = tiles.offset
    # MUSEdata[:,3:6] = tiles.size
    # MUSEdata[:,6:9] = rMatrix[:,0,:] #n1
    # MUSEdata[:,9:12] = rMatrix[:,1,:] #n2
    # MUSEdata[:,12:15] = rMatrix[:,2,:] #n3
    # # MUSEdata[:,15:18] = tiles.M * tiles.size.prod(1)
    # MUSEdata[:,15:18] = tiles.M/np.linalg.norm(tiles.M, axis=1).reshape((-1,1))
    # MUSEdata[:,18:] = 0
    
    FICUSdata[:,:3] = tiles.offset
    FICUSdata[:,3:6] = rMatrix[:,0,:] 
    FICUSdata[:,6:9] = rMatrix[:,1,:] 
    FICUSdata[:,9:11] = tiles.size[:,:2]
    FICUSdata[:,11] = np.linalg.norm(tiles.M, axis=1)
    FICUSdata[:,12:] = tiles.M/np.linalg.norm(tiles.M, axis=1).reshape((-1,1))
    
    # np.savetxt("magTense"+fname,MUSEdata,delimiter=',',header=header0)
    np.savetxt(fname,FICUSdata,delimiter=',',header=header1) 


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--muea", type=float, default=1.00)
    parser.add_argument("--muoa", type=float, default=1.00)
    args = parser.parse_args()

    rmin    =  0.200000000000000E+000
    rmax    =  0.400000000000000E+000
    zmin    = -0.100000000000000E+000
    zmax    =  0.100000000000000E+000
    pmin    =  0.000000000000000E+000
    pmax    =   6.28318530717959     
    nr      =         301
    nz      =         301
    nphi    =         72
    nfp = 2

    # generate coordinates
    rr = np.linspace(rmin, rmax, nr)
    zz = np.linspace(zmin, zmax, nz)
    phi = np.linspace(0, 2 * np.pi / nfp, nphi, endpoint=False)
    gridr, gridp, gridz = np.meshgrid(rr, phi, zz, indexing="ij")
    gridx = gridr * np.cos(gridp)
    gridy = gridr * np.sin(gridp)
    targets = np.transpose([gridx.ravel(), gridy.ravel(), gridz.ravel()])
    np.save("./Intermediate/evaluation_points_301", targets)

    #Finite mu calculation
    mu0 = 4 * np.pi * 1e-7
    prism = muse2magntense("./Input/magtense_zot80_3d.csv", magnetization=1.1658E6, mu=[1.05, 1.05])
    # mu_arr = [[1.00, 1.00], [1.05, 1.05], [1.05, 1.15]]
    mu = [args.muea, args.muoa]
    print("mu = ", mu)
    prism.set_mu_r_ea(mu[0])
    prism.set_mu_r_oa(mu[1])
    (updated_tiles, H_mu) = magtense.run_simulation(prism, targets, plot=False)
    
    #Magnetization updates & H output
    np.save("./Intermediate/bfield_grid_301_{0:.2f}_{1:.2f}".format(mu[0],mu[1]), H_mu*mu0)
    np.save("./Intermediate/Tiles_{0:.2f}_{1:.2f}".format(mu[0],mu[1]),updated_tiles)
    updateTiles2FICUS(updated_tiles,"./PM_FICUS/FICUS_zot80_3d_{0:.2f}_{1:.2f}.csv".format(mu[0],mu[1]))
    print("Write FICUS PM Input:", "./PM_FICUS/FICUS_zot80_3d_{0:.2f}_{1:.2f}.csv".format(mu[0],mu[1]))


if __name__=='__main__':
    main()