from pathlib import Path
import argparse, numpy as np
import time
import copy
import magtense, magtense.magstatics as _ms
from magtense.magstatics import get_demag_tensor, get_H_field
from coilpy import rotation_matrix, muse2magntense
from scipy.constants import mu_0
from scipy.optimize import minimize
from scipy.linalg import cholesky, solve_triangular, LinAlgError, svd, norm, solve
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import svds, minres, gmres, LinearOperator

# For loaing coil fields
from simsopt.util.permanent_magnet_helper_functions import read_focus_coils
from simsopt.field import Coil, BiotSavart

# LA - Directsolve Test
import scipy.linalg as _la

## Creating new float type
_F64 = dict(dtype=np.float64)

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
# todo: remove this since macromagnetic model + device scale magnets don't need anisostropic term
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
    Return H (A m⁻¹).  If TF-coils are requested we evaluate the Biot Savart
    field with SIMSOPT; otherwise we fall back to a user-supplied uniform H.
    """
    if use_coils and bs_coil is not None:
        bs_coil.set_points(np.asarray(pts))
        return bs_coil.B() / mu_0
    # fallback uniform (This is just returning 0 or some arbitrary constants. Allows us to do the GUT check)
    N = pts.shape[0]
    H = np.broadcast_to(const_H, (N,3)).copy()
    return H

# Inspired by MUSE magstatics run simulation -> Iterate Magnet Solution (hard case) fortran code 
# LINK: https://github.com/cmt-dtu-energy/MagTense/blob/0dbf1442170fd734e73496e3af458bd0e3b2d8ef/source/DemagField/DemagField/IterateMagnetSolution.f90#L32
def sor_sweep(tiles, centres, demag_tensor, Ms, K, omega, max_it, tol, Hcoil_func, log_every, demag_only=False):
    lambda_ = float(1.0)
    H_prev  = np.zeros((tiles.n, 3), **_F64)
    M_prev = None
    maxDiff = [0.0, 0.0, 0.0, 0.0]

    for it in range(1, max_it + 1):
        # update magnetization from relaxed field
        u = tiles.u_ea
        H_par  = (H_prev * u).sum(1, **_F64)[:, None] * u
        H_perp = H_prev - H_par
        chi_par  = (tiles.mu_r_ea - 1)[:, None]
        chi_perp = (tiles.mu_r_oa - 1)[:, None]
        M_rem  = tiles.M_rem[:, None] * u
        tiles.M = M_rem + (chi_par * H_par) + (chi_perp * H_perp)

        # compute new field terms
        H_demag = get_H_field(tiles, centres, demag_tensor).astype(np.float64) 
        H_anis  = np.zeros_like(H_demag).astype(np.float64)
        # H_anis = anisotropy_field(tiles.M, u, K, Ms) if not demag_only else np.zeros_like(H_demag) 
        H_coil  = Hcoil_func(centres).astype(np.float64)                if not demag_only else np.zeros_like(H_demag).astype(np.float64)
        H_eff   = H_demag + H_anis + H_coil

        # successive over-relaxation
        H_prev = H_prev + lambda_ * (H_eff - H_prev)

        # convergence on |M|
        if M_prev is None:
            err = float('inf')
        else:
            delta_m = np.abs(np.linalg.norm(tiles.M, axis=1) - np.linalg.norm(M_prev, axis=1))
            rel = delta_m / np.maximum(np.linalg.norm(M_prev, axis=1), 1e-30)
            err = rel.max()

        # adapt relaxation parameter
        maxDiff = [err] + maxDiff[:3]
        if it > 3:
            c1 = maxDiff[0] > maxDiff[1] < maxDiff[2] > maxDiff[3]
            c2 = maxDiff[3] > maxDiff[2] > maxDiff[1] > maxDiff[0]
            c3 = maxDiff[3] > maxDiff[2] > maxDiff[1] < maxDiff[0]
            if c1 or c2 or c3:
                lambda_ *= 0.5

        # compute angle between M and field
        m_hat = tiles.M / np.linalg.norm(tiles.M, axis=1, keepdims=True)
        h_hat = H_eff / np.maximum(np.linalg.norm(H_eff, axis=1, keepdims=True), 1e-30)
        max_angle_MH = np.degrees(np.arccos(np.clip((m_hat * h_hat).sum(1), -1, 1))).max()

        # logging
        if it == 1 or it % log_every == 0:
            print(
                f"  iter {it:3d}: lambda={lambda_:.3f} | max‖ΔM‖/‖M‖={err:.2e} max∠MH={max_angle_MH:6.2f}° "
                f"| max|H_demag|={np.linalg.norm(H_demag, axis=1).max():.2e} "
                f"| max|H_anis|={np.linalg.norm(H_anis, axis=1).max():.2e} "
                f"| max|H_coil|={np.linalg.norm(H_coil, axis=1).max():.2e} "
                f"| max|H_eff|={np.linalg.norm(H_eff, axis=1).max():.2e}"
            )

        if err < tol * lambda_ and it > 3:
            print(f"  converged in {it} iterations (err={err:.2e})")
            break

        M_prev = tiles.M.copy()

    return tiles


# Used this for more comprehensive direct solve which included logging of different methods will remove this but wanted to add this git history
def magstatics_direct_solve_v1(tiles, centres, demag_tensor, ms, k, Hcoil_func, demag_only=False, krylov_tol=1e-8, krylov_it=200):
    n = tiles.n
    u = tiles.u_ea
    chi_pa = tiles.mu_r_ea - 1
    chi_pe = tiles.mu_r_oa - 1
    m_rem  = tiles.M_rem
    chi = np.empty((n, 3, 3), **_F64)
    for i in range(n):
        uuT = np.outer(u[i], u[i])
        chi[i] = chi_pa[i]*uuT + chi_pe[i]*(np.eye(3) - uuT)

    H_a = Hcoil_func(centres).astype(np.float64) if not demag_only else np.zeros((n, 3))
    b   = np.vstack([m_rem[i]*u[i] + chi[i] @ H_a[i] for i in range(n)]).ravel()

    N3 = 3*n
    A  = np.zeros((N3, N3), **_F64)
    I3 = np.eye(3)
    for i in range(n):
        r = slice(3*i, 3*i+3)
        # assembling dense system matrix...
        # demag_tensor holds –NM, so I – chi_i @ demag_tensor = I + chi_i @ NM (which is what we want)
        A[r, r] = I3 - chi[i] @ demag_tensor[i, i]
        for j in range(n):
            if i == j:
                continue
            c = slice(3*j, 3*j+3)
            A[r, c] = -chi[i] @ demag_tensor[i, j]
    nnz = (A !=0).sum()
    density = nnz / (N3*N3)
    sym = np.allclose(A, A.T, atol=1e-12, rtol=1e-10)

    #  iterative sparse SVD only if it's really sparse and large... 
    if density <= 0.10:
        s_min = svds(csc_matrix(A), k=1, which='SM', return_singular_vectors=False)[0]
        s_max = svds(csc_matrix(A), k=1, which='LM', return_singular_vectors=False)[0]
    else:
        s = svd(A, compute_uv=False)
        s_min, s_max = s[-1], s[0]
    print(f"density {density:.2%}, condition numbe = {s_max/s_min:.2e}, sym: {sym}")

    x = None
    spd = False
    if sym:
        try:
            t0 = time.perf_counter()
            L = cholesky(A, lower=True, check_finite=False)
            spd = True
            y = solve_triangular(L, b, lower=True, check_finite=False)
            x = solve_triangular(L.T, y,lower=False,check_finite=False)
            err_fact = norm(A - L @ L.T) / norm(A)
            print(f"Chol {time.perf_counter()-t0:.2f}s, res {norm(A@x-b)/norm(b):.2e}, LLT err {err_fact:.2e}")
        except LinAlgError:
            pass
    if x is None:
        t0 = time.perf_counter()
        x  = solve(A, b, assume_a='gen')
        print(f"LU {time.perf_counter()-t0:.2f}s, res {norm(A@x-b)/norm(b):.2e}")
    linop = LinearOperator(A.shape, matvec=lambda v: A@v, rmatvec=lambda v: A.T@v)
    if sym:
        t0 = time.perf_counter()
        xk, info = minres(linop, b, tol=krylov_tol, maxiter=krylov_it)
        print(f"MINRES {time.perf_counter()-t0:.2f}s, res {norm(A@xk-b)/norm(b):.2e}, info {info}")
    else:
        t0 = time.perf_counter()
        xk, info = gmres(linop, b, maxiter=krylov_it)
        print(f"GMRES  {time.perf_counter()-t0:.2f}s, res {norm(A@xk-b)/norm(b):.2e}, info {info}")

    tiles.M = x.reshape(n, 3)
    return tiles


def construct_A(n, chi, demag_tensor):
    N3 = 3*n
    A  = np.zeros((N3, N3), **_F64)
    I3 = np.eye(3)
    for i in range(n):
 
        mats = -np.einsum('pr,jrq->jpq', chi[i], demag_tensor[i], optimize=True)
        mats[i] += I3
        r = slice(3*i, 3*i+3)
        for j in range(n):
            c = slice(3*j, 3*j+3)
            A[r, c] = mats[j]
    return A

def construct_b(n, m_rem, u, chi, H_a):
    return np.vstack([m_rem[i]*u[i] + chi[i] @ H_a[i] for i in range(n)]).ravel()

def magstatics_direct_solve(tiles, centres, demag_tensor, ms, k, Hcoil_func, demag_only=False,krylov_tol=1e-3, krylov_it=200):
    n      = tiles.n
    u      = tiles.u_ea
    chi_pa = tiles.mu_r_ea - 1
    chi_pe = tiles.mu_r_oa - 1
    m_rem  = tiles.M_rem

    chi = np.empty((n, 3, 3), **_F64)
    for i in range(n):
        uuT    = np.outer(u[i], u[i])
        chi[i] = chi_pa[i]*uuT + chi_pe[i]*(np.eye(3) - uuT)

    H_a = Hcoil_func(centres).astype(np.float64) if not demag_only else np.zeros((n, 3))
    b   = construct_b(n, m_rem, u, chi, H_a)
    A  = construct_A(n, chi, demag_tensor)
    
    sym = np.allclose(A, A.T, atol=1e-12, rtol=1e-10)
    print(f"symmetry check: {sym}")

    linop = LinearOperator(A.shape, matvec=lambda v: A @ v,
                                      rmatvec=lambda v: A.T @ v)

    if sym:
        x, info = minres(linop, b, rtol=krylov_tol, maxiter=krylov_it)
        print(f"MINRES   res  {norm(A@x-b)/norm(b):.2e}, info {info}")
    else:
        x, info = gmres( linop, b, rtol=krylov_tol, maxiter=krylov_it )
        print(f"GMRES    res  {norm(A@x-b)/norm(b):.2e}, info {info}")


    tiles.M = x.reshape(n, 3)
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
    ap.add_argument("--max-it",type=int,   default=100)
    ap.add_argument("--tol",   type=float, default=1e-10)
    ap.add_argument("--Hcoil", nargs=3, type=float, default=[0,0,0], metavar=("Hx","Hy","Hz"), help="Uniform coil field components A/m") 
    ap.add_argument("--use-coils",action="store_true",help="add TF-coil field with FICUS analytic model")
    ap.add_argument("--demag-only", action="store_true", help="only include H_demag (zero out H_anis and H_coil)")
    ap.add_argument("--coil-file",type=Path,help="numpy .npy or csv with columns x y z nx ny nz I a " "(ignored unless --use-coils)")
    ap.add_argument("--solver", choices=["direct","sor","compare"], default="direct", help="Which magnetostatic solver to run (default: direct)")

    args = ap.parse_args()

    tag = f"{args.muea:.2f}_{args.muoa:.2f}"
    Path("Intermediate").mkdir(exist_ok=True)
    Path("PM_FICUS").mkdir(exist_ok=True)

    # Step 1. Loadding magnets and coils.... (from completed ficus optimization)
    # Extra note used magnetization as defined in def muse2magntense(muse_file, mu=(1.05, 1.05), magnetization=1.16e6, **kwargs): source code
    # Link: https://zhucaoxiang.github.io/CoilPy/_modules/coilpy/magtense_interface.html#muse2magntense
    tiles = muse2magntense("./Input/magtense_zot80_3d.csv", magnetization=1.1658e6, mu=[args.muea, args.muoa])    
    tiles.mu_r_ea, tiles.mu_r_oa = args.muea, args.muoa
    centres = tiles.offset.copy()

    # Using values from one of the magtense papers (Rasmus, Roberto 2023) 
    # Link: https://backend.orbit.dtu.dk/ws/portalfiles/portal/313199195/1_s2.0_S0304885323001592_main_1_.pdf
    Ms = np.float64(1.281e6)  # A/m  (st. μ0 Ms=1.61T)
    K  = np.float64(4.3e6) # J/m³ 

    #Precompute tensor
    demag_tensor = get_demag_tensor(tiles, centres).astype(np.float64)

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
    if args.solver == "sor":
        tiles = sor_sweep(
            tiles, centres, demag_tensor,
            Ms, K,
            omega=args.omega,
            max_it=args.max_it,
            tol=args.tol,
            Hcoil_func=Hcoil_func,
            log_every=5,
            demag_only=args.demag_only
        )
        solver_name = "SOR"
    elif args.solver == "direct":
        tiles = magstatics_direct_solve(
            tiles, centres, demag_tensor,
            Ms, K,
            Hcoil_func=Hcoil_func,
            demag_only=args.demag_only
        )
        solver_name = "Direct"
    elif args.solver == "compare": 
        n = tiles.n
        u      = tiles.u_ea
        chi_pa = tiles.mu_r_ea - 1
        chi_pe = tiles.mu_r_oa - 1
        m_rem  = tiles.M_rem
        chi = np.empty((n, 3, 3), **_F64)
        for i in range(n):
            uuT    = np.outer(u[i], u[i])
            chi[i] = chi_pa[i]*uuT + chi_pe[i]*(np.eye(3) - uuT)
        H_a = Hcoil_func(centres).astype(np.float64) if not args.demag_only else np.zeros((n, 3))
    
        b = construct_b(n, m_rem, u, chi, H_a)
        A = construct_A(n, chi, demag_tensor)

        tiles_sor = copy.deepcopy(tiles)
        tiles_sor = sor_sweep(tiles_sor, centres, demag_tensor, Ms, K,
                              omega=args.omega, max_it=args.max_it,
                              tol=args.tol, Hcoil_func=Hcoil_func,
                              log_every=5, demag_only=args.demag_only)
        res_sor = norm(A @ tiles_sor.M.ravel() - b)
        tiles_dir = copy.deepcopy(tiles)
        tiles_dir = magstatics_direct_solve(tiles_dir, centres, demag_tensor,
                                            Ms, K, Hcoil_func,
                                            demag_only=args.demag_only)
        res_dir = norm(A @ tiles_dir.M.ravel() - b) 

        print(f"SOR residual    ||A·M–b||/||b|| = {res_sor:.2e}")
        print(f"Direct residual ||A·M–b||/||b|| = {res_dir:.2e}")

        tiles = tiles_dir
        solver_name = "COMPARE"

    t1 = time.perf_counter()
    print(f"\nTotal {solver_name} solve time: {t1-t0:.2f} s\n")

    np.save(f"Intermediate/Tiles_{tag}.npy", tiles)
    print(f"→  Intermediate/Tiles_{tag}.npy written")
    export_ficus_csv(tiles, f"PM_FICUS/FICUS_{tag}.csv")
    print(f"→  PM_FICUS/FICUS_{tag}.csv written")

if __name__ == "__main__":
    main()
