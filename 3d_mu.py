from pathlib import Path
import argparse
import numpy as np
from simsopt.field import DipoleField
from simsopt.geo import SurfaceRZFourier, Surface
from coilpy import muse2magntense
import magtense.magstatics as _ms
import magtense  

def _bind(name, fn):
    setattr(_ms.Tiles, name, fn)

_bind("set_tile_type",   lambda s, v    : setattr(s, "tile_type", v))
_bind("set_size_i",      lambda s, v, i : s.size.__setitem__((i,), v))
_bind("set_offset_i",    lambda s, v, i : s.offset.__setitem__((i,), v))
_bind("set_rotation_i",  lambda s, v, i : s.rot.__setitem__((i,), v))
_bind("set_remanence_i", lambda s, v, i : s.M_rem.__setitem__((i,), v))
_bind("set_mu_r_ea_i",   lambda s, v, i : s.mu_r_ea.__setitem__((i,), v))
_bind("set_mu_r_oa_i",   lambda s, v, i : s.mu_r_oa.__setitem__((i,), v))
_bind("set_mag_angle_i", lambda s, v, i : s.set_easy_axis(val=v, idx=i))
_bind("set_color_i",     lambda s, v, i : s.color.__setitem__((i,), v))
magtense.Tiles = _ms.Tiles

def csv_brick_volumes(csv_path: Path) -> np.ndarray:
    # Not sure about htis part
    # read half-lengths (cols 3–5) from MagTense CSV to compute volumes (I think)
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    hlw = data[:, 3:6]
    return 8.0 * np.prod(hlw, axis=1)

def magnetisation_to_moment(M: np.ndarray, hlw: np.ndarray) -> np.ndarray:
    vol = 8.0 * np.prod(hlw, axis=1)
    return M * vol[:, None]

def write_vertex_vtp(pts, vals, name, fname):
    n = len(pts)
    with open(fname, "w") as f:
        f.write(
            '<?xml version="1.0"?>\n'
            '<VTKFile type="PolyData" byte_order="LittleEndian">\n'
            f' <PolyData>\n  <Piece NumberOfPoints="{n}" NumberOfVerts="{n}">\n'
            '   <Verts>\n'
            '    <DataArray type="Int32" Name="connectivity" format="ascii">\n'
            + " ".join(map(str, range(n))) +
            "\n    </DataArray>\n"
            '    <DataArray type="Int32" Name="offsets" format="ascii">\n'
            + " ".join(str(i+1) for i in range(n)) +
            "\n    </DataArray>\n"
            '   </Verts>\n'
            '   <Points>\n'
            '    <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for x, y, z in pts:
            f.write(f"{x:.10e} {y:.10e} {z:.10e}\n")
        f.write(
            "</DataArray>\n   </Points>\n"
            f'<PointData Scalars="{name}">\n'
            '<DataArray type="Float64" Name="' + name +
            '" format="ascii">\n')
        for v in vals:
            f.write(f"{v:.10e}\n")
        f.write(
            "</DataArray>\n   </PointData>\n"
            "</Piece>\n </PolyData>\n</VTKFile>\n")

def write_surface_deltaBn(surf, B1, B2, ntheta, nphi, fname, label):
    normals = surf.unitnormal().reshape(-1, 3)
    delta = np.einsum("ij,ij->i", B1 - B2, normals).reshape((ntheta, nphi)).T[:, :, None]
    surf.to_vtk(fname, extra_data={label: delta})
    print(f"Wrote {fname} | max‖{label}‖ = {np.abs(delta).max():.3e} T")
    print("max|B_ficus| =", np.linalg.norm(B2,axis=1).max())
    print("max|B_SOR|   =", np.linalg.norm(B1, axis=1).max())
    rel = np.abs(delta).max() / np.linalg.norm(B2, axis=1).max()
    print(f"max|ΔBn| / max|B_ficus| = {rel:.2%}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--muea",          type=float, required=True)
    p.add_argument("--muoa",          type=float, required=True)
    p.add_argument("--ref-muea",      type=float)
    p.add_argument("--ref-muoa",      type=float)
    p.add_argument("--surface-file",  default="./scripts/launch/input.5pga19")
    p.add_argument("--ntheta",        type=int, default=28)
    p.add_argument("--nphi",          type=int, default=49)
    p.add_argument("--add-ficus-comparison", action="store_true")
    p.add_argument(
        "--csv-file",
        default="./Input/magtense_zot80_3d.csv",
        help="path to original MagTense CSV for brick geometry"
    )
    args = p.parse_args()

    single = args.ref_muea is None or args.ref_muoa is None
    tag = f"{args.muea:.2f}_{args.muoa:.2f}"
    ref_tag = None if single else f"{args.ref_muea:.2f}_{args.ref_muoa:.2f}"

    out = Path("Output")
    out.mkdir(exist_ok=True)

    vols = csv_brick_volumes(Path(args.csv_file))
    tiles_new = np.load(Path("Intermediate") / f"Tiles_{tag}.npy", allow_pickle=True).item()
    SLICE_T = 1.5875e-3
    xyz = tiles_new.offset.copy()
    hlw_tiles = tiles_new.size / 2.0
    vol_slice = (2 * hlw_tiles[:, 0]) * (2 * hlw_tiles[:, 1]) * SLICE_T
    m_new = tiles_new.M * vol_slice[:, None]

    b_new = DipoleField(xyz, m_new.ravel(), nfp=2, coordinate_flag="cartesian")
    surf = SurfaceRZFourier.from_vmec_input(
        Path(args.surface_file),
        range=Surface.RANGE_FULL_TORUS,
        ntheta=args.ntheta,
        nphi=args.nphi
    )
    pts = surf.gamma().reshape(-1, 3)

    b_new.set_points(pts)
    b_new._toVTK(out / f"dipoles_mu_{tag}")
    print(f"Wrote dipoles_mu_{tag}")

    if not single:
        tiles_ref = np.load(Path("Intermediate") / f"Tiles_{ref_tag}.npy", allow_pickle=True).item()
        m_ref = magnetisation_to_moment(tiles_ref.M, hlw_tiles)
        b_ref = DipoleField(xyz, m_ref.ravel(), nfp=2, coordinate_flag="cartesian")
        b_ref.set_points(pts)
        b_ref._toVTK(out / f"dipoles_mu_{ref_tag}")

        Bn, Br = b_new.B(), b_ref.B()
        dB = (np.linalg.norm(Bn, axis=1) - np.linalg.norm(Br, axis=1))
        dB = dB.reshape((args.ntheta, args.nphi)).T[:, :, None]
        surf.to_vtk(out / f"surface_deltaB_mu_{tag}", extra_data={"Δ|B|": dB})
        print(f"Wrote surface_deltaB_mu_{tag}")

        dot = np.einsum("ij,ij->i", tiles_new.M, tiles_ref.M)
        # Compute norms per tile (axis=1) instead of global flatten
        norms_new = np.linalg.norm(tiles_new.M, axis=1) # per‐row 2‐norm
        norms_ref = np.linalg.norm(tiles_ref.M, axis=1) # per‐row 2‐norm
        norm_prod = norms_new * norms_ref # elementwise product

        norm_prod[norm_prod < 1e-12] = 1.0
        dtheta = np.degrees(np.arccos(np.clip(dot / norm_prod, -1, 1)))
        write_vertex_vtp(xyz, dtheta, "delta_theta_deg", out / f"magnet_delta_theta_mu_{tag}.vtp")
        print(f"Wrote magnet_delta_theta_mu_{tag}")

    if args.add_ficus_comparison:
        tiles_fic = muse2magntense(
            "./Input/magtense_zot80_3d.csv",
            magnetization=1.1658e6,
            mu=[1.00, 1.00]
        )
        
        m_fic = tiles_fic.M * vol_slice[:, None]
        
        b_fic = DipoleField(
            xyz,
            m_fic.ravel(),
            nfp=2, coordinate_flag="cartesian"
        )
        b_fic.set_points(pts)
        write_surface_deltaBn(
            surf,
            b_new.B(),
            b_fic.B(),
            args.ntheta,
            args.nphi,
            out / (f"surface_deltaBn_SORminusFICUS_{tag}"),
            "ΔBn"
        )

    print("Finished; outputs in", out)

if __name__ == "__main__":
    main()
