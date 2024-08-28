import sys
import os
import logging
from typing import Tuple, Optional

sys.path.append(os.getenv('DELPHI_PATH', '/home/dockeruser/delphi/'))
from DELPHI import *
from SOLEILII_parameters.SOLEILII_TDR_parameters import *
from scipy.constants import pi
import pandas as pd
import h5py as hp
from machine_data.soleil import v2366, v2366_v2, v2366_v3
from tqdm import tqdm
from utils import get_parser_for_delphi
from scipy.interpolate import interp1d

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def read_impedance(filename: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read impedance data from a file."""
    try:
        f = pd.read_csv(filename, sep="\t", header=None)
        valid_data = f[0] > 0
        freqs = f[0][valid_data]
        Zre = f[1][valid_data]
        Zim = f[2][valid_data]
        Z = np.vstack((Zre, Zim)).T
        return freqs, Z
    except Exception as e:
        logger.error(f"Failed to read impedance file: {filename}, Error: {e}")
        return np.array([]), np.array([])

def get_and_interpolate_sigmas(file_path: str, current: float, is_sigma_z: bool = False) -> Tuple[Optional[float], Optional[float]]:
    """Interpolate sigma values from the input file based on the current."""
    try:
        df_input = np.loadtxt(file_path, delimiter='\t', usecols=(1, 2, 3), dtype=np.float64, skiprows=1)
        (I, sigmas, sigmadp) = df_input.T
        sigma = interp1d(I * 1e-3, sigmas, kind="linear")(current)
        sigma_dp = interp1d(I * 1e-3, sigmadp, kind="linear")(current)
        return sigma, sigma_dp
    except FileNotFoundError:
        logger.warning(f"No valid filename for {'bunch length' if is_sigma_z else 'impedance'}.")
        return None, None

def process_eigenmodes(
    nx: int,
    Ib: np.ndarray,
    taub: float,
    ring,
    plane: str,
    Z: np.ndarray,
    freqs: np.ndarray,
    Qp: float,
    damper_gain: float,
    impedance_multiplier: int,
    M: int,
    Q_s: float,
    g0,
    a,
    b,
    coefdamper,
    coefZ
) -> pd.DataFrame:
    """Process eigenmodes and return results as a DataFrame."""
    results = pd.DataFrame(columns=[
        "BunchCurrent", "eigvals_re", "eigvals_im", "Qp", "BunchLength", "nx",
        "id_state", "damper_gain", "damper_phase", "impedance_multiplier"
    ])
    eigvecs_list = []

    for i, bunch_current in enumerate(tqdm(Ib)):
        Nb = bunch_current / e / ring.f0
        eigvals, eigvecs, *_ = eigenmodesDELPHI_converged(
            nx=nx, M=M, omegaksi=Qp / ring.ac * ring.omega0, omega0=ring.omega0,
            tunefrac=(ring.tune[1 if plane == "vertical" else 0] - np.floor(ring.tune[1 if plane == "vertical" else 0])),
            a=a, b=b, taub=taub, g=g0, Z=Z, kmax=5, freqZ=freqs, coefdamper=coefdamper,
            coefZ=coefZ, omegas=Q_s * ring.omega0, flageigenvect=False, optimized_convergence=True,
            lmaxold=-1, nmaxold=-1, crit=5e-2, abseps=1e-3
        )
        for k in range(len(eigvals)):
            result = pd.DataFrame({
                "BunchCurrent": bunch_current, "eigvals_re": np.real(eigvals[k]) / ring.omega0,
                "eigvals_im": np.imag(eigvals[k]) / ring.omega0, "Qp": Qp, "BunchLength": taub,
                "nx": nx, "id_state": ring.IDs, "damper_gain": damper_gain, "damper_phase": 0,
                "impedance_multiplier": impedance_multiplier
            }, index=[k])
            results = pd.concat([results, result], ignore_index=True)
        eigvecs_list.append(eigvecs)

    return results

def run_bunch_current_scan(
    impedance_filename="Zydip.dat", id_state="open", sigma_z_filename="sigmaz.npy",
    sigma_z=9e-12, plane="vertical", Qp=0.0, min_value=0.1, max_value=5.0,
    n_scan_points=50, M=1, Q_s=0.0021, n_max=0, n_min=0, damper_gain=0,
    impedance_multiplier=2
) -> pd.DataFrame:
    """Run a single bunch current scan and save results to CSV."""
    freqs, Z = read_impedance(os.path.join("/home/dockeruser/delphi/input/", impedance_filename))
    Z *= impedance_multiplier
    ring = v2366_v3(IDs=id_state)
    tune = ring.tune[0] if plane == "horizontal" else ring.tune[1]
    Ib = 1e-3 * np.linspace(min_value, max_value, n_scan_points)
    sigma, sigma_dp = get_and_interpolate_sigmas(sigma_z_filename, Ib[0], is_sigma_z=True)
    taub = 4 * (sigma if sigma else sigma_z)
    g0, a, b = longdistribution_decomp(taub, typelong="Gaussian")

    logger.info("Running bunch current scan...")
    results = process_eigenmodes(
        nx=n_min, Ib=Ib, taub=taub, ring=ring, plane=plane, Z=Z, freqs=freqs,
        Qp=Qp, damper_gain=damper_gain, impedance_multiplier=impedance_multiplier,
        M=M, Q_s=Q_s, g0=g0, a=a, b=b, coefdamper=0, coefZ=0
    )

    results.to_csv(
        f"/home/dockeruser/delphi/data/delphi(sigma_z={taub/4:.1e},ID={id_state},plane={plane},Qp={Qp},M={M},Q_s={Q_s:.1e},n_max={n_max},damper_gain={damper_gain:.1e},ximpedance={impedance_multiplier:.1f}).csv",
        sep="\t"
    )
    return results

def run_chroma_scan(
    impedance_filename="Zydip.dat", sigma_z_filename="sigmaz.npy", id_state="open",
    sigma_z=9e-12, plane="vertical", Ib=1.2e-3, min_value=0.1, max_value=5.0,
    n_scan_points=50, M=1, Q_s=0.0021, n_max=0, n_min=0, damper_gain=0,
    impedance_multiplier=2
) -> pd.DataFrame:
    """Run a chromaticity scan and save results to CSV."""
    freqs, Z = read_impedance(os.path.join("/home/dockeruser/delphi/input/", impedance_filename))
    Z *= impedance_multiplier
    ring = v2366_v2(IDs=id_state)
    tune = ring.tune[0] if plane == "horizontal" else ring.tune[1]
    Qp_range = np.linspace(min_value, max_value, n_scan_points)
    sigma, sigma_dp = get_and_interpolate_sigmas(sigma_z_filename, Ib, is_sigma_z=True)
    taub = 4 * (sigma if sigma else sigma_z)
    g0, a, b = longdistribution_decomp(taub, typelong="Gaussian")

    logger.info("Running chromaticity scan...")
    results = pd.DataFrame(columns=[
        "Qp", "BunchCurrent", "eigvals_re", "eigvals_im", "BunchLength", "nx",
        "id_state", "damper_gain", "damper_phase", "impedance_multiplier"
    ])

    for Qp in tqdm(Qp_range):
        result = process_eigenmodes(
            nx=n_min, Ib=[Ib], taub=taub, ring=ring, plane=plane, Z=Z, freqs=freqs,
            Qp=Qp, damper_gain=damper_gain, impedance_multiplier=impedance_multiplier,
            M=M, Q_s=Q_s, g0=g0, a=a, b=b, coefdamper=0, coefZ=0
        )
        results = pd.concat([results, result], ignore_index=True)

    results.to_csv(
        f"/home/dockeruser/delphi/data/delphi(sigma_z={taub/4:.1e},ID={id_state},plane={plane},Ib={Ib:.1e},M={M},Q_s={Q_s:.1e},n_max={n_max},damper_gain={damper_gain:.1e},ximpedance={impedance_multiplier:.1f}).csv",
        sep="\t"
    )
    return results

if __name__ == "__main__":
    parser = get_parser_for_delphi()
    args = parser.parse_args()
    run_bunch_current_scan(
        impedance_filename=args.impedance_filename,
        id_state=args.id_state,
        sigma_z_filename=args.sigma_z_filename,
        sigma_z=args.sigma_z,
        plane=args.plane,
        Qp=args.Qp,
        min_value=args.min_value,
        max_value=args.max_value,
        n_scan_points=args.n_scan_points,
        M=args.M,
        Q_s=args.Q_s,
        n_max=args.n_max,
        n_min=args.n_min,
        damper_gain=args.damper_gain,
        impedance_multiplier=args.impedance_multiplier
    )
