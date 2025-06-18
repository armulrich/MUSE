import argparse
from pathlib import Path
import numpy as np
from simsopt.field import DipoleField
from simsopt.geo import SurfaceRZFourier, Surface
from coilpy import muse2magntense
import magtense, magtense.magstatics as _ms
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

####
# Generate 3D VTK visualizations of dipole fields and "optionally" compare (using flags) two different μ‐values.
#
# Modes:
# - Single‐μ mode (default):  Produce a .vtu glyph file of the dipole field for the given μ.
# - Comparison mode (if --ref-muea or--ref-muoa are supplied):  
# In addition to the single‐μ mode output also:
#    – export the reference‐μ glyph,
#    – compute and write Δ|B| on the plasma surface (.vtu),
#    – compute per‐magnet Δθ of M vectors and write as a .vtp heat map.
####

def write_vertex_vtp(pts, values, name, out_path):
    n = len(pts)
    with open(out_path, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n')
        f.write(f' <PolyData>\n  <Piece NumberOfPoints="{n}" NumberOfVerts="{n}">\n')
        f.write('   <Verts>\n    <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        f.write(" ".join(map(str, range(n))) + "\n    </DataArray>\n")
        f.write('    <DataArray type="Int32" Name="offsets" format="ascii">\n')
        f.write(" ".join(str(i+1) for i in range(n)) + "\n    </DataArray>\n   </Verts>\n")
        f.write('   <Points>\n    <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for x,y,z in pts:
            f.write(f"     {x:.8e} {y:.8e} {z:.8e}\n")
        f.write("    </DataArray>\n   </Points>\n")
        f.write(f'   <PointData Scalars="{name}">\n<DataArray type="Float64" Name="{name}" format="ascii">\n')
        for v in values:
            f.write(f"     {v:.8e}\n")
        f.write("    </DataArray>\n   </PointData>\n")
        f.write("  </Piece>\n </PolyData>\n</VTKFile>\n")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--muea",  type=float, required=True)
    p.add_argument("--muoa",  type=float, required=True)
    p.add_argument("--ref-muea", type=float) # -> If comparison mode
    p.add_argument("--ref-muoa", type=float)
    p.add_argument("--surface-file", default="./scripts/launch/input.5pga19")
    p.add_argument("--ntheta",type=int, default=28) # Using these default values from previous MUSE code
    p.add_argument("--nphi", type=int, default=49)
    args = p.parse_args()

    single = (args.ref_muea is None or args.ref_muoa is None)
    tag = f"{args.muea:.2f}_{args.muoa:.2f}"
    ref_tag = f"{args.ref_muea:.2f}_{args.ref_muoa:.2f}" if not single else None

    intermediate = Path("Intermediate")
    out_dir = Path("Output")
    out_dir.mkdir(exist_ok=True)

    tiles_new = np.load(intermediate / f"Tiles_{tag}.npy", allow_pickle=True).item()
    xyz, m_new = tiles_new.offset, tiles_new.M
    b_new = DipoleField(xyz, m_new.ravel(), nfp=2, coordinate_flag="cartesian") 
    #Not 100% sure if cartesian is the right flag here but ouput seems correct

    surface = SurfaceRZFourier.from_vmec_input(
        Path(args.surface_file),
        range=Surface.RANGE_FULL_TORUS,
        nphi=args.nphi, ntheta=args.ntheta
    )
    pts = surface.gamma().reshape(-1,3)

    # SIDE STEP: also dump the original (pre-iteration) dipoles from the FICUS input
    tiles_orig = muse2magntense("./Input/magtense_zot80_3d.csv",magnetization=1.1658e6,mu=[args.ref_muea, args.ref_muoa])
    xyz0, m0 = tiles_orig.offset, tiles_orig.M
    b_orig    = DipoleField(xyz0, m0.ravel(), nfp=2, coordinate_flag="cartesian")
    b_orig.set_points(pts)
    out0 = out_dir / "dipoles_original.vtu"
    b_orig._toVTK(out0)
    print("Wrote", out0)

    # Contiun
    out1 = out_dir / f"dipoles_mu_{tag}.vtu"
    b_new.set_points(pts)
    b_new._toVTK(out1)
    print("Wrote",out1)

    if not single:

        tiles_ref = np.load(intermediate / f"Tiles_{ref_tag}.npy", allow_pickle=True).item()
        m_ref     = tiles_ref.M
        b_ref     = DipoleField(xyz, m_ref.ravel(), nfp=2, coordinate_flag="cartesian")
        out2 = out_dir / f"dipoles_mu_{ref_tag}.vtu"
        b_ref.set_points(pts)
        b_ref._toVTK(out2)
        print("Wrote",out2)

        B_new = b_new.B(); B_ref = b_ref.B()
        mag_new = np.linalg.norm(B_new,axis=1)
        mag_ref = np.linalg.norm(B_ref,axis=1)
        deltaB = (mag_new - mag_ref).reshape((args.ntheta,args.nphi)).T[:,:,None]
        out3 = out_dir / f"surface_deltaB_mu_{tag}.vtu"
        surface.to_vtk(out3, extra_data={"ΔB": deltaB})
        print("Wrote",out3)

        dot = np.einsum('ij,ij->i', m_new, m_ref) # inner product
        norms = np.linalg.norm(m_new,axis=1)*np.linalg.norm(m_ref,axis=1)
        norms = np.where(norms>1e-12, norms,1.0)
        deltaθ = np.degrees(np.arccos(np.clip(dot/norms,-1,1)))
        out4 = out_dir / f"magnet_delta_theta_mu_{tag}.vtp"
        write_vertex_vtp(xyz, deltaθ, "delta_theta_deg", out4)
        print("Wrote", out4)

    mode = "single-mu" if single else "comparison"
    print(f"Finished ({mode} mode). Outputs in {out_dir}")

if __name__=="__main__":
    main()
