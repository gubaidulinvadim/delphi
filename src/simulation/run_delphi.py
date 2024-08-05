import sys, os

sys.path.append('/home/dockeruser/delphi/')
from DELPHI import *
from SOLEILII_parameters.SOLEILII_TDR_parameters import *
from scipy.constants import pi
import pandas as pd
import h5py as hp
from machine_data.soleil import v2366, v2366_v2
from tqdm import tqdm
from utils import get_parser_for_delphi
from scipy.interpolate import interp1d


def read_impedance(filename):
    f = pd.read_csv(filename, sep="\t", header=None)
    freqs = f[0][f[0] > 0]
    Zre = f[1][f[0] > 0]
    Zim = f[2][f[0] > 0]
    Z = np.vstack((Zre, Zim)).T
    return freqs, Z

def get_and_interpolate_sigmas(file_path, current, is_sigma_z=False):
    try:
        df_input = np.loadtxt(file_path, delimiter='\t', usecols=(1, 2, 3), dtype=np.float64, skiprows=1)
        (I, sigmas, sigmadp) = df_input.T
        sigma = interp1d(I*1e-3, sigmas, kind="linear")(current)
        sigma_dp = interp1d(I*1e-3, sigmadp, kind="linear")(current)
        return sigma, sigma_dp
    except FileNotFoundError:
        print(f"No valid filename for {'bunch length' if is_sigma_z else 'impedance'}.")
        return None, None

def run_bunch_current_scan(
    impedance_filename="Zydip.dat",
    id_state="open",
    sigma_z_filename="sigmaz.npy",
    sigma_z=9e-12,
    plane="vertical",
    Qp=0.0,
    min_value=0.1,
    max_value=5.0,
    n_scan_points=50,
    M=1,
    Q_s=0.0021,
    n_max=0,
    n_min=0,
    damper_gain=0,
    impedance_multiplier=2,
):
    freqs, Z = read_impedance("/home/dockeruser/delphi/input/" + impedance_filename)
    Z *= impedance_multiplier
    ring = v2366_v2(IDs=id_state)
    tune = ring.tune[0] if plane == "horizontal" else ring.tune[1]

    taub = 4 * sigma_z
    Ib = 1e-3 * np.linspace(min_value, max_value, n_scan_points)
    omega_ksi = Qp / ring.ac * ring.omega0
    kmax = 5
    results = pd.DataFrame(columns=[
        "BunchCurrent",
        "eigvals_re",
        "eigvals_im",
        "Qp",
        "BunchLength",
        "nx",
        "id_state",
        "damper_gain",
        "damper_phase",
        "impedance_multiplier",
            ])
    eigvecs_list = []
    for nx in range(n_min, n_max + 1):
        print(f'Coupled-bunch mode {nx}')
        for i, bunch_current in enumerate(tqdm(Ib)):
            Nb = bunch_current / e / ring.f0
            # try:
            #     df_input = np.loadtxt('/home/dockeruser/delphi/input/'+sigma_z_filename, delimiter='\t', usecols=(1, 2, 3), dtype=np.float64, skiprows=1)
            #     (I, sigmas, sigmadp) = df_input.T
            #     sigma = interp1d(I*1e-3, sigmas, kind="linear")(bunch_current)
            #     sigma_dp = interp1d(I*1e-3, sigmadp, kind="linear")(bunch_current)
            #     taub = 4 * sigma
            #     Q_s = ring.ac*sigma_dp/sigma/ring.omega0
            #     print('BUNCH LENGTH IS ', taub)
            # except:
            #     print(
            #         "No valid filename for bunch length. Using default bunch length {:.1e} ps"
            #         .format(sigma_z / 1e-12))
            #     taub = 4 * sigma_z
            sigma, sigma_dp = get_and_interpolate_sigmas(sigma_z_filename, bunch_current, is_sigma_z=True)
            if sigma is None:
                sigma_ = sigma_z
            taub = 4 * sigma
            g0, a, b = longdistribution_decomp(taub, typelong="Gaussian")
            print(g0, a, b)
            damper_phase = 0#-np.pi/2
            coefdamper, coefZ = computes_coef(
                f0=ring.f0,
                dmax=damper_gain,
                b=b,
                g0=g0[0],
                dnormfactor=1,
                taub=taub,
                dphase=damper_phase,
                M=M,
                Nb=Nb,
                gamma=ring.gamma,
                Q=tune,
                particle="electron",
            )
            beta_mbtrack2 = (
                ring.optics.local_beta[1]
                if plane == "vertical"
                else ring.optics.local_beta[0]
            )
            beta_delphi = ring.L / tune / 2 / np.pi
            coefZ *= beta_mbtrack2 / beta_delphi
            (
                eigvals,
                eigvecs,
                lmaxold,
                nmaxold,
                max_freq_old,
                matdamperold,
                matZold,
            ) = eigenmodesDELPHI_converged(
                nx=nx,
                M=M,
                omegaksi=omega_ksi,
                omega0=ring.omega0,
                tunefrac=(tune - np.floor(tune)),
                a=a,
                b=b,
                taub=taub,
                g=g0,
                Z=Z,
                kmax=kmax,
                freqZ=freqs,
                coefdamper=coefdamper,
                # flag_trapz=True,
                coefZ=coefZ,
                omegas=Q_s * ring.omega0,
                flageigenvect=False,
                optimized_convergence=True,
                lmaxold=-1,
                nmaxold=-1,
                crit=5e-2,
                abseps=1e-3,
            )
            for k in range(len(eigvals)):
                result = pd.DataFrame(
                    {
                        "BunchCurrent": bunch_current,
                        "eigvals_re": np.real(eigvals[k]) / ring.omega0,
                        "eigvals_im": np.imag(eigvals[k]) / ring.omega0,
                        "Qp": Qp,
                        "BunchLength": taub,
                        "nx": nx,
                        "id_state": id_state,
                        "damper_gain": damper_gain,
                        "damper_phase": damper_phase,
                        "impedance_multiplier": impedance_multiplier
                    },
                    index=[k],
                )
                results = pd.concat([results, result], ignore_index=True)
            eigvecs_list.append(eigvecs)
    results.to_csv(
        path_or_buf=
        "/home/dockeruser/delphi/data/delphi(sigma_z={:.1e},ID={:},plane={:},Qp={:},M={:},Q_s={:.1e},n_max={:},damper_gain={:.1e},ximpedance={:.1f}).csv"
        .format(taub/4, id_state, plane, Qp, M, Q_s, n_max, damper_gain, impedance_multiplier),
        sep="\t",
    )
    # np.save('delphi_eigvecs(sigma_z={:.1e},plane={:},Qp={:}).npy'.format(
    # sigma_z, plane, Qp), np.array(eigvecs_list, dtype=object))
    return results


def run_chroma_scan(
    impedance_filename="Zydip.dat",
    sigma_z_filename="sigmaz.npy",
    id_state="open",
    sigma_z=9e-12,
    plane="vertical",
    Ib=1.2e-3,
    min_value=0.1,
    max_value=5.0,
    n_scan_points=50,
    M=1,
    Q_s=0.0021,
    n_max=0,
    n_min=0,
    damper_gain=0,
    impedance_multiplier=1,
):
    freqs, Z = read_impedance("/home/dockeruser/delphi/input/" + impedance_filename)
    Z *= impedance_multiplier
    ring = v2366_v2(IDs=id_state)
    tune = ring.tune[0] if plane == "horizontal" else ring.tune[1]
    taub = 4 * sigma_z
    kmax = 5
    results = pd.DataFrame(columns=[
        "BunchCurrent",
        "eigvals_re",
        "eigvals_im",
        "Qp",
        "BunchLength",
        "nx",
        "id_state",
    ])

    chroma = np.linspace(min_value, max_value, n_scan_points)
    Nb = Ib / e / ring.f0
    eigvecs_list = []
    for nx in range(n_min, n_max + 1):
        print(f'Coupled-bunch mode {nx}')
        for i, Qp in enumerate(tqdm(chroma)):
            try:
                df_input = np.loadtxt('/home/dockeruser/delphi/input/'+sigma_z_filename, delimiter='\t', usecols=(1, 2, 3), dtype=np.float64, skiprows=1)
                (I, sigmas, sigmadp) = df_input.T
                sigma = interp1d(I*1e-3, sigmas, kind="linear")(Ib)
                sigma_dp = interp1d(I*1e-3, sigmadp, kind="linear")(Ib)
                Q_s = ring.ac*sigma_dp/sigma/ring.omega0
                taub = 4 * sigma
            except:
                print(
                    "No valid filename for bunch length. Using default bunch length {:.1e} ps"
                    .format(sigma_z / 1e-12))
                taub = 4 * sigma_z
            g0, a, b = longdistribution_decomp(taub, typelong="Gaussian")
            omega_ksi = Qp / ring.ac * ring.omega0
            damper_phase = 0#-np.pi/2
            coefdamper, coefZ = computes_coef(
                f0=ring.f0,
                dmax=damper_gain,
                b=b,
                g0=g0,
                dnormfactor=1,
                taub=taub,
                dphase=damper_phase,
                M=M,
                Nb=Nb,
                gamma=ring.gamma,
                Q=tune,
                particle="electron",
            )
            beta_mbtrack2 = (
                ring.optics.local_beta[1]
                if plane == "vertical"
                else ring.optics.local_beta[0]
            )
            beta_delphi = ring.L / tune / 2 / np.pi
            coefZ *= beta_mbtrack2 / beta_delphi
            (
                eigvals,
                eigvecs,
                lmaxold,
                nmaxold,
                max_freq_old,
                matdamperold,
                matZold,
            ) = eigenmodesDELPHI_converged(
                nx=nx,
                M=M,
                omegaksi=omega_ksi,
                omega0=ring.omega0,
                tunefrac=(tune - np.floor(tune)),
                a=a,
                b=b,
                taub=taub,
                g=g0,
                Z=Z,
                kmax=kmax,
                freqZ=freqs,
                coefdamper=coefdamper,
                coefZ=coefZ,
                omegas=Q_s * ring.omega0,
                flageigenvect=False,
                optimized_convergence=True,
                lmaxold=-1,
                nmaxold=-1,
                crit=5e-2,
                abseps=1e-3,
            )
            for k in range(len(eigvals)):
                result = pd.DataFrame(
                    {
                        "bunch_current": Ib,
                        "eigvals_re": np.real(eigvals[k]) / ring.omega0,
                        "eigvals_im": np.imag(eigvals[k]) / ring.omega0,
                        "chromaticity": Qp,
                        "rms_bunch_length": taub/4,
                        "coupled_bunch_mode": nx,
                        "n_bunches": M,
                        "id_state": id_state,
                        "damper_gain": damper_gain,
                        "damper_phase": damper_phase,
                        "impedance_multiplier": impedance_multiplier
                    },
                    index=[k]
                )
                results = pd.concat([results, result], ignore_index=True)
            eigvecs_list.append(eigvecs)
    results.to_csv(
        path_or_buf=
        "/home/dockeruser/delphi/data/delphi(sigma_z={:.1e},id={:},plane={:},Ib={:.1e},M={:},Q_s={:.1e},n_max={:},damper_gain={:.1e},impedance_multiplier={:.1f}).csv"
        .format(taub/4, id_state, plane, Ib, M, Q_s, n_max, damper_gain, impedance_multiplier),
        sep="\t",
    )
    # np.save('delphi_eigvecs(sigma_z={:.1e},plane={:},Ib={:}).npy'.format(
    # sigma_z, plane, Ib), np.array(eigvecs_list, dtype=object))
    return results


if __name__ == "__main__":
    os.environ['KMP_DUPLICATE_LIB_OK']='True'
    parser = get_parser_for_delphi()
    args = parser.parse_args()
    if args.scan_type == "sb_chromaticity":
        run_chroma_scan(
            impedance_filename=args.filename,
            id_state=args.id_state,
            sigma_z=args.sigma_z,
            plane=args.plane,
            Ib=args.current,
            min_value=args.min_value,
            max_value=args.max_value,
            n_scan_points=args.n_scan_points,
            Q_s=args.Q_s,
            M=1,
            sigma_z_filename=args.sigmas_filename,
            damper_gain=args.damper_gain,
            impedance_multiplier=args.impedance_multiplier
        )
    elif args.scan_type == "sb_current":
        print("running a single bunch current scan")
        run_bunch_current_scan(
            impedance_filename=args.filename,
            id_state=args.id_state,
            sigma_z=args.sigma_z,
            plane=args.plane,
            Qp=args.chromaticity,
            min_value=args.min_value,
            max_value=args.max_value,
            n_scan_points=args.n_scan_points,
            Q_s=args.Q_s,
            sigma_z_filename=args.sigmas_filename,
            damper_gain=args.damper_gain,
            impedance_multiplier=args.impedance_multiplier
        )
    elif args.scan_type == "mb_chromaticity":
        run_chroma_scan(
            impedance_filename=args.filename,
            id_state=args.id_state,
            sigma_z=args.sigma_z,
            plane=args.plane,
            Ib=args.current,
            min_value=args.min_value,
            max_value=args.max_value,
            n_scan_points=args.n_scan_points,
            M=args.n_max+1,
            n_max=args.n_max,
            n_min=0,
            Q_s=args.Q_s,
            sigma_z_filename=args.sigmas_filename,
            damper_gain=args.damper_gain,
            impedance_multiplier=args.impedance_multiplier

        )
    elif args.scan_type == "mb_current":
        run_bunch_current_scan(
            impedance_filename=args.filename,
            id_state=args.id_state,
            sigma_z=args.sigma_z,
            plane=args.plane,
            Qp=args.chromaticity,
            min_value=args.min_value,
            max_value=args.max_value,
            n_scan_points=args.n_scan_points,
            M=args.n_max+1,
            n_max=args.n_max,
            n_min=0,
            Q_s=args.Q_s,
            sigma_z_filename=args.sigmas_filename,
            damper_gain=args.damper_gain,
            impedance_multiplier=args.impedance_multiplier  

        )
    else:
        print("Scan type is not supported.")
