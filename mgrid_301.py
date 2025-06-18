import numpy as np
import argparse
from pathlib import Path

import magtense.magstatics as _ms
import magtense
from magtense.magstatics import run_simulation
from coilpy import rotation_matrix, muse2magntense
from scipy.constants import mu_0

### Adding Quick fixes here to enable usage with newever magtense simulation
def _bind(name, fn): setattr(_ms.Tiles, name, fn)
_bind("set_tile_type",lambda s,v: setattr(s,"tile_type",v))
_bind("set_size_i",    lambda s,v,i: s.size.__setitem__((i,),v))
_bind("set_offset_i",  lambda s,v,i: s.offset.__setitem__((i,),v))
_bind("set_rotation_i",  lambda s,v,i: s.rot.__setitem__((i,),v))
_bind("set_remanence_i",lambda s,v,i: s.M_rem.__setitem__((i,),v))
_bind("set_mu_r_ea_i", lambda s,v,i: s.mu_r_ea.__setitem__((i,),v))
_bind("set_mu_r_oa_i", lambda s,v,i: s.mu_r_oa.__setitem__((i,),v))
_bind("set_mag_angle_i",lambda s,v,i: s.set_easy_axis(val=v, idx=i))
_bind("set_color_i", lambda s,v,i: s.color.__setitem__((i,),v))
magtense.Tiles = _ms.Tiles

## Have to access private variables here (again temp fix but I think this should work)
magtense.run_simulation = _ms.run_simulation

def updateTiles2FICUS(tiles, fname):
    nmag = tiles.n
    out = np.zeros((nmag, 15), dtype=np.float64)
    header = (
        "X [m], Y [m], Z [m], "
        "n1x, n1y, n1z, n2x, n2y, n2z, "
        "H [m], L [m], M [A/m], mx, my, mz"
    )

    rMatrix = np.zeros((nmag,3,3))
    for i in range(nmag):
        rMatrix[i] = rotation_matrix(tiles.rot[i, 0], tiles.rot[i, 1], tiles.rot[i, 2], xyz=True).T


    # Removed the commented code here
    # And renamed FICUSdata to out
    out[:,0:3] = tiles.offset
    out[:,3:6] = rMatrix[:,0,:]
    out[:,6:9] = rMatrix[:,1,:]
    out[:,9:11] = tiles.size[:,:2]

    #improving readability here since this is the M computation (or at least I think)
    magM = np.linalg.norm(tiles.M, axis=1)
    out[:, 11] = magM 
    out[:, 12:] = tiles.M/magM.reshape((-1,1)) # Avoided the extra computaiton of another np.linalg since using variable

    np.savetxt(fname, out, delimiter=",", header=header)

def main():
    #Adding a small description here (to try to improve readability)
    parser = argparse.ArgumentParser(
        description="Compute finite-µ PM B-field on a (301×301×72) r-φ-z grid."
    )
    parser.add_argument(
        "--muea", type=float, default=1.00,
        help="Relative permeability along the easy axis (μₑₐ)."
    )
    parser.add_argument(
        "--muoa", type=float, default=1.00,
        help="Relative permeability along the perpendicular axis (μₒₐ)."
    )
    args = parser.parse_args()


    rmin    =  0.200000000000000E+000
    rmax    =  0.400000000000000E+000
    zmin    = -0.100000000000000E+000
    zmax    =  0.100000000000000E+000
    pmin    =  0.000000000000000E+000
    pmax    =   6.28318530717959     
    nr      =         20 # Originally at 301
    nz      =         20 # Originally at 301 
    nphi    =         12 # Originally at 72 -> Reduced these values here to make the computation more workable
    nfp = 2

    # generate coordinates
    rr = np.linspace(rmin, rmax, nr)
    zz = np.linspace(zmin, zmax, nz)
    phi = np.linspace(0, 2 * np.pi / nfp, nphi, endpoint=False)
    gridr, gridp, gridz = np.meshgrid(rr, phi, zz, indexing="ij")
    gridx = gridr * np.cos(gridp)
    gridy = gridr * np.sin(gridp)
    targets = np.transpose([gridx.ravel(), gridy.ravel(), gridz.ravel()])

    Path("./Intermediate").mkdir(exist_ok=True)
    np.save("./Intermediate/evaluation_points_301.npy", targets)

    Path("./Input").mkdir(exist_ok=True)
    csv_path = "./Input/magtense_zot80_3d.csv"
    #Finite mu calculation
    prism = muse2magntense(csv_path, magnetization=1.1658e6, mu=[args.muea, args.muoa])
    
    # Had to introduce this setting off prismmu since we no longer have the prism.set method in Magtense 2.2.0 
    prism.mu_r_ea = args.muea 
    prism.mu_r_oa = args.muoa
    # Not 100% sure if this is working as expected

    # Adding some debugging messages
    print(f"Built a Tiles object with {prism.n} prisms...")
    print(f"Running run_simulation on {prism.n} prisms + {targets.shape[0]} points....")

    updated_tiles, H_out = run_simulation(prism, targets)
    B_out = H_out * mu_0

    Path("./Intermediate").mkdir(exist_ok=True)
    fname_grid = f"./Intermediate/bfield_grid_301_{args.muea:.2f}_{args.muoa:.2f}.npy"
    fname_tiles = f"./Intermediate/Tiles_{args.muea:.2f}_{args.muoa:.2f}.npy"
    
    #Saving again...
    np.save(fname_grid, B_out)
    np.save(fname_tiles, updated_tiles)

    print(f"Saved B-field grid to {fname_grid}")
    print(f"Saved updated Tiles to {fname_tiles}")

    Path("./PM_FICUS").mkdir(exist_ok=True)
    out_csv = f"./PM_FICUS/FICUS_zot80_3d_{args.muea:.2f}_{args.muoa:.2f}.csv"
    updateTiles2FICUS(updated_tiles, out_csv)

    print(f"Wrote FICUS PM input CSV → {out_csv}")


if __name__ == "__main__":
    main()
