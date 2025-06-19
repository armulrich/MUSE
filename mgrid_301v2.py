"""
Micromagnetic LaBonte–SOR sweep for the 9 736 Nd–Fe–B blocks in
    ./Input/magtense_zot80_3d.csv

Field terms at every tile centre
  H_demag   – MagTense (all dipolar interactions,incl. self)
  H_anis   – 2K/(μ0 Ms²) (m·u) u (uniaxial)
  H_coil    – coil field (default,0)
"""
from pathlib import Path
import argparse, numpy as np
import time
import magtense, magtense.magstatics as _ms
from magtense.magstatics import get_demag_tensor, get_H_field
from coilpy import rotation_matrix, muse2magntense
from scipy.constants import mu_0
from scipy.optimize import minimize

# For loaing coil fields
from simsopt.util.permanent_magnet_helper_functions import read_focus_coils
from simsopt.field import Coil, BiotSavart

## This small section is devoted to fixes I incorporated to make the magtense -> coilpy inteaction work correctly (I already tried to reduce this)
### Context: this is because coilpy was written with an older version of magtense in mind ###
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


## coil field reader helper
def load_coil_BiotSavart(path: Path) -> BiotSavart:
    """
    Read a *FOCUS* ‘tf-coils.focus’ file with SIMSOPT and return a
    `simsopt.field.BiotSavart` object (already wired with the correct
    currents)  
    """
    curves, currents, ncoils = read_focus_coils(str(path))
    # build Coil objects
    coils = [Coil(curves[i], currents[i]) for i in range(ncoils)]
    bs = BiotSavart(coils)
    return bs

## field helper funcs
def anisotropy_field(M, u_ea, K, Ms):
    """
    Local uniaxial anisotropy
        H_anis_i = 2 K/(μ0 Ms²)·(m_i·u_i)·u_i
    where u_i = tiles.u_ea[i].
    """
    proj = np.einsum('ij,ij->i', M, u_ea)           # (m·u) per tile
    return (2*K/(mu_0*Ms**2)) * proj[:, None] * u_ea  # (N,3)

def coil_field_at(pts, const_H, use_coils, bs_coil: BiotSavart):
    """
    Return H (A m⁻¹).  If TF-coils are requested we evaluate the Biot–Savart
    field with SIMSOPT; otherwise we fall back to a user-supplied uniform H.
    """
    if use_coils and bs_coil is not None:
        bs_coil.set_points(np.asarray(pts))
        return bs_coil.B() / mu_0
    # fallback uniform (This is just returning 0 or some arbitrary constants. Allows us to do the GUT check)
    N = pts.shape[0]
    H = np.broadcast_to(const_H, (N,3)).copy()
    return H

## Labonte: accelerated iteration (note that I am recommputing the H terms at each iterate)
# Adopted this from "A Numerical Study of LaBonte’s Iteration: An Approach to Acceleration" (Kosavisutte, Hayaslii, 1996)
# Link: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=539347
def labonte_sweep(tiles, centres, demag_tensor, Ms, K, omega, max_it, tol, Hcoil_func, log_every, demag_only=False):
    
    chi = tiles.mu_r_ea[0] - 1.0 # Linear scaling since iterate magnetization does not have the mu values we need (just an approixmation)
    stall_counter = 0
    last_res      = np.inf

    for it in range(1, max_it+1):

        H_demag = get_H_field(tiles, centres, demag_tensor)
        H_anis  = anisotropy_field(tiles.M, tiles.u_ea, K, Ms) if not demag_only else np.zeros_like(H_demag)
        H_coil  = Hcoil_func(centres)                           if not demag_only else np.zeros_like(H_demag)

        # chi‐correction on the internal field H_int = H_demag + H_anis + H_coil
        if chi > 0:
            m_hat   = tiles.M / np.linalg.norm(tiles.M, axis=1, keepdims=True)
            H_int   = H_demag + H_anis + H_coil
            tiles.M = Ms * m_hat + chi * H_int
            # must re‐compute H_demag after changing M
            H_demag = get_H_field(tiles, centres, demag_tensor)

        H_eff = H_demag + H_anis + H_coil

        norm = np.linalg.norm(H_eff, axis=1, keepdims=True)
        norm[norm < 1e-12] = 1.0
        h_hat = H_eff / norm

        m_old  = tiles.M.copy() # Deep copy here
        tiles.M = (1-omega)*m_old + omega*Ms*h_hat # Labonte 

        # compute field magnitudes
        h_abs = np.linalg.norm(H_eff, axis=1)
        
        # get a safe epsilon: machine epsilon scaled by the largest field
        eps = np.finfo(h_abs.dtype).eps * max(h_abs.max(), 1.0)
        
        # clip so nothing is below eps
        h_safe = np.clip(h_abs, eps, None)
        
        # Using cross product rule....
        res = np.max(
            np.linalg.norm(np.cross(tiles.M, H_eff), axis=1)
            / (Ms * h_safe)
        )

        if it==1 or it%log_every==0:
            dot  = np.einsum('ij,ij->i', m_old, tiles.M) / (np.linalg.norm(m_old,axis=1) * np.linalg.norm(tiles.M,axis=1))
            ## Adding some iter logging
            max_tilt = np.degrees(np.arccos(np.clip(dot,-1,1))).max()
            print(
                f"  iter {it:3d}: residual={res:.2e} | max tilt={max_tilt:6.2f}° "
                f"| max|H_demag|={np.linalg.norm(H_demag,axis=1).max():.2e} "
                f"| max|H_anis|={np.linalg.norm(H_anis,axis=1).max():.2e} "
                f"| max|H_coil|={np.linalg.norm(H_coil,axis=1).max():.2e} "
                f"| max|H_eff|={np.linalg.norm(H_eff,axis=1).max():.2e}")

                    
        # simple stagnation guard (DOUBLE GUARD)
        if res < tol:
            print(f"  converged in {it} iterations (residual={res:.2e})")
            break
        if abs(res-last_res)/res < 1e-4:   # So a ... <0.01 % improvement
            stall_counter += 1
            if stall_counter >= 5:
                print(f"  stalled after {it} iterations (residual≈{res:.2e})")
                break
        else:
            stall_counter = 0
        last_res = res
    return tiles

def export_ficus_csv(tiles, fname):
    n = tiles.n
    out = np.zeros((n,15))
    R   = np.zeros((n,3,3))
    for i in range(n):
        R[i] = rotation_matrix(*tiles.rot[i], xyz=True).T
    out[:, :3] = tiles.offset
    out[:, 3:6] = R[:,0,:]
    out[:, 6:9] = R[:,1,:]
    out[:, 9:11] = tiles.size[:,:2]
    out[:,11] = np.linalg.norm(tiles.M,axis=1)
    out[:,12:] = tiles.M/np.linalg.norm(tiles.M,axis=1)[:,None]
    np.savetxt(fname,out,delimiter=',',header="X[m],Y[m],Z[m],n1x,n1y,n1z,n2x,n2y,n2z,H[m],L[m],M[A m],mx,my,mz")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--muea",  type=float, default=1.00)
    ap.add_argument("--muoa",  type=float, default=1.00)
    ap.add_argument("--omega", type=float, default=0.3, help="SOR over-relaxation factor w (≈0.5–1.9)")
    ap.add_argument("--max-it",type=int,   default=50)
    ap.add_argument("--tol",   type=float, default=1e-4)
    ap.add_argument("--Hcoil", nargs=3, type=float, default=[0,0,0], metavar=("Hx","Hy","Hz"), help="Uniform coil field components A/m") 
    ap.add_argument("--use-coils",action="store_true",help="add TF-coil field with FICUS analytic model")
    ap.add_argument("--demag-only", action="store_true", help="only include H_demag (zero out H_anis and H_coil)")
    ap.add_argument("--coil-file",type=Path,help="numpy .npy or csv with columns x y z nx ny nz I a " "(ignored unless --use-coils)")
    args = ap.parse_args()

    tag = f"{args.muea:.2f}_{args.muoa:.2f}"
    Path("Intermediate").mkdir(exist_ok=True)
    Path("PM_FICUS").mkdir(exist_ok=True)

    # Step 1. Loadding magnets and coils.... (from completed ficus optimization)
    # Extra note used magnetization as defined in def muse2magntense(muse_file, mu=(1.05, 1.05), magnetization=1.16e6, **kwargs): source code
    # Link: https://zhucaoxiang.github.io/CoilPy/_modules/coilpy/magtense_interface.html#muse2magntense
    tiles = muse2magntense("./Input/magtense_zot80_3d.csv", magnetization=1.1658e6, mu=[args.muea, args.muoa])
    
    #tiles.magnet_type[:] = 3  #3 = soft + constant μ_r
    #  Note: I think this has no effect and will internally be treated as 1->hard since we do not invoke iterate magnizaton. 
    
    tiles.mu_r_ea, tiles.mu_r_oa = args.muea, args.muoa

    centres = tiles.offset.copy()

    # Using values from one of the magtense papers (Rasmus, Roberto 2023) 
    # Link: https://backend.orbit.dtu.dk/ws/portalfiles/portal/313199195/1_s2.0_S0304885323001592_main_1_.pdf
    Ms = 1.281e6 # A/m  (st. μ0 Ms=1.61T)
    K  = 4.3e6  # J/m³

    #Precompute tensor
    demag_tensor = get_demag_tensor(tiles, centres)


    #Loading coils...
    #intended to be used with: https://github.com/tmqian/Muse-Design-Paper/blob/main/FAMUS-files/tf-coils.focus
    bs_coil = None
    if args.use_coils:
        if args.coil_file is None:
            ap.error("--coil-file must be given with --use-coils")
        bs_coil = load_coil_BiotSavart(args.coil_file)
        print(f"TF-coil geometry loaded from {args.coil_file}")
    
    const_H = np.asarray(args.Hcoil,dtype=float)
    Hcoil_func = lambda pts: coil_field_at(pts, const_H, args.use_coils, bs_coil)
    
    # Step 2. Iteration until Browns 1st Eqn. Is satisfied 
    t0 = time.perf_counter()
    tiles = labonte_sweep(tiles,centres, demag_tensor,Ms,K,args.omega,args.max_it,args.tol,Hcoil_func,log_every=5, demag_only=args.demag_only)
    t1 = time.perf_counter()
    print(f"\nTotal LaBonte-Accelerated iteration time: {t1-t0:.2f} s\n")

    np.save(f"Intermediate/Tiles_{tag}.npy", tiles)
    print(f"→  Intermediate/Tiles_{tag}.npy written")
    export_ficus_csv(tiles, f"PM_FICUS/FICUS_{tag}.csv")
    print(f"→  PM_FICUS/FICUS_{tag}.csv written")

if __name__ == "__main__":
    main()
