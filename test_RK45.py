import argparse
import numpy as np
from magtense.micromag import MicromagProblem
## Just a testing script to check how Magtense RK45 scheme works

def rk45_norm_drift_tol(tol, thres):
    prob = MicromagProblem(
        res=[2,2,1],
        grid_L=(1e-6,1e-6,1e-6),
        alpha=0.1,
        gamma=2.211e5,
        cuda=False,
        cvode=True,
        tol=tol,
        thres=thres,
    )
    prob.m0 = None  # random unit initial state
    h_ext = lambda t: np.zeros((len(t), 3))
    t, M_out, *_ = prob.run_simulation(
        t_end=1e-9, nt=10,
        fct_h_ext=h_ext, nt_h_ext=2
    )

    M = M_out[..., 0, :]
    norms = np.linalg.norm(M, axis=2)
    return np.max(np.abs(norms - 1.0))

def rk45_norm_drift_scale(scale, tol = 1e-4, thres = 1e-6):
    prob = MicromagProblem(
        res=[2,2,1],
        grid_L=(1e-6,1e-6,1e-6),
        alpha=0.1,
        gamma=2.211e5,
        cuda=False,
        cvode=False,
        tol=tol,
        thres=thres,
    )
    prob.m0 = None        # random unit
    prob.m0 *= scale      # blow it up
    h_ext = lambda t: np.zeros((len(t), 3))
    t, M_out, *_ = prob.run_simulation(
        t_end=1e-9, nt=10,
        fct_h_ext=h_ext, nt_h_ext=2
    )
    M = M_out[..., 0, :]
    norms = np.linalg.norm(M, axis=2)
    return np.max(np.abs(norms - 1.0))

def main():
    parser = argparse.ArgumentParser(
        description="test RK45 norm‐conservation under varying tol/threshold or Scale"
    )
    parser.add_argument(
        "mode", choices=["tolerance","scale"],
        help="Which test to run"
    )
    args = parser.parse_args()

    if args.mode == "tolerance":
        for tol, thres in [(1e-4,1e-6), (1e-5,1e-7), (1e-6,1e-8)]:
            drift = rk45_norm_drift_tol(tol, thres)
            print(f"tol={tol:.1e}, thres={thres:.1e} → max‖m‖−1 = {drift:.3e}")

    else: 
        for scale in (1.0, 1e3, 1e100):
            try:
                drift = rk45_norm_drift_scale(scale)
                print(f"scale={scale:>4.1e} → max‖m‖−1 = {drift:.3e}")
            except Exception as e:
                print(f"scale={scale:>4.1e} → integration failed ({e})")

if __name__ == "__main__":
    main()