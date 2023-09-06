from DELPHI import *
from SOLEILII_parameters.SOLEILII_TDR_parameters import *
from scipy.constants import pi
import pandas as pd
from machine_data.soleil import v2366
from tqdm import tqdm
from utils import get_parser_for_delphi


def read_impedance(filename):
    f = pd.read_csv('Zydip.dat', sep='\t', header=None)
    freqs = f[0][f[0] > 0]
    Zre = f[1][f[0] > 0]
    Zim = f[2][f[0] > 0]
    Z = np.vstack((Zre, Zim)).T
    return freqs, Z


def run_bunch_current_scan(impedance_filename, ID_state='open', sigma_z=9e-12, plane='vertical', Qp=0.0, min_value=0.1, max_value=5.0, n_scan_points=50):
    freqs, Z = read_impedance(impedance_filename)
    ring = v2366(IDs=ID_state, load_lattice=False)

    taub = 4*sigma_z
    g0, a, b = longdistribution_decomp(taub, typelong='Gaussian')
    Ib = 1e-3*np.linspace(min_value, max_value, n_scan_points)
    omega_ksi = Qp/ring.ac*ring.omega0
    kmax = 5
    results = pd.DataFrame(columns=['BunchCurrent',
                                    'eigvals_re',
                                    'eigvals_im',
                                    'Qp',
                                    'FinalBunchLength'])
    eigvecs_list = []
    for bunch_current in tqdm(Ib):
        Nb = bunch_current/e/ring.f0
        coefdamper, coefZ = computes_coef(f0=ring.f0,
                                          dmax=0,
                                          b=b,
                                          g0=g0,
                                          dnormfactor=np.infty,
                                          taub=taub,
                                          dphase=0,
                                          M=1,
                                          Nb=Nb,
                                          gamma=ring.gamma,
                                          Q=ring.tune[1],
                                          particle='electron')
        eigvals, eigvecs, lmaxold, nmaxold, max_freq_old, matdamperold, matZold = eigenmodesDELPHI_converged(
            nx=0,
            M=1,
            omegaksi=omega_ksi,
            omega0=ring.omega0,
            tunefrac=(ring.tune[1]-np.floor(ring.tune[1])),
            a=a,
            b=b,
            taub=taub,
            g=g0,
            Z=Z,
            kmax=kmax,
            freqZ=freqs,
            coefdamper=coefdamper,
            coefZ=coefZ,
            omegas=ring.synchrotron_tune(1.8e6)*ring.omega0,
            flageigenvect=True,
            optimized_convergence=True,
            lmaxold=-1,
            nmaxold=-1,
            crit=5e-2,
            abseps=1e-3
        )
        for k in range(len(eigvals)):
            result = pd.DataFrame({
                'BunchCurrent': bunch_current,
                'eigvals_re': np.real(eigvals[k])/ring.omega0,
                'eigvals_im': np.imag(eigvals[k])/ring.omega0,
                'Qp': Qp,
                'BunchLength': sigma_z
            }, index=[k])
            results = pd.concat([results, result], ignore_index=True)
    results.to_csv(path_or_buf='delphi(sigma_z={:.1e},plane={:},Qp={:}).csv'.format(
        sigma_z, plane, Qp), sep='\t')

    np.save('delphi_eigvecs(sigma_z={:.1e},plane={:},Ib={:}).npy'.format(
        sigma_z, plane, Ib), np.array(eigvecs_list))
    return results


def run_chroma_scan(impedance_filename, ID_state='open', sigma_z=9e-12, plane='vertical', Ib=1.2e-3, min_value=0.1, max_value=5.0, n_scan_points=50):
    freqs, Z = read_impedance(impedance_filename)
    ring = v2366(IDs=ID_state, load_lattice=False)
    taub = 4*sigma_z
    g0, a, b = longdistribution_decomp(taub, typelong='Gaussian')
    kmax = 5
    results = pd.DataFrame(columns=['BunchCurrent',
                                    'eigvals_re',
                                    'eigvals_im',
                                    'Qp',
                                    'FinalBunchLength'])
    chroma = np.linspace(min_value, max_value, n_scan_points)
    Nb = Ib/e/ring.f0
    eigvecs_list = []
    for Qp in tqdm(chroma):
        omega_ksi = Qp/ring.ac*ring.omega0
        coefdamper, coefZ = computes_coef(f0=ring.f0,
                                          dmax=0,
                                          b=b,
                                          g0=g0,
                                          dnormfactor=np.infty,
                                          taub=taub,
                                          dphase=0,
                                          M=1,
                                          Nb=Nb,
                                          gamma=ring.gamma,
                                          Q=ring.tune[1],
                                          particle='electron')
        eigvals, eigvecs, lmaxold, nmaxold, max_freq_old, matdamperold, matZold = eigenmodesDELPHI_converged(
            nx=0,
            M=1,
            omegaksi=omega_ksi,
            omega0=ring.omega0,
            tunefrac=(ring.tune[1]-np.floor(ring.tune[1])),
            a=a,
            b=b,
            taub=taub,
            g=g0,
            Z=Z,
            kmax=kmax,
            freqZ=freqs,
            coefdamper=coefdamper,
            coefZ=coefZ,
            omegas=ring.synchrotron_tune(1.8e6)*ring.omega0,
            flageigenvect=True,
            optimized_convergence=True,
            lmaxold=-1,
            nmaxold=-1,
            crit=5e-2,
            abseps=1e-3
        )
        for k in range(len(eigvals)):
            result = pd.DataFrame({
                'BunchCurrent': Ib,
                'eigvals_re': np.real(eigvals[k])/ring.omega0,
                'eigvals_im': np.imag(eigvals[k])/ring.omega0,
                'Qp': Qp,
                'BunchLength': sigma_z
            }, index=[k])
            results = pd.concat([results, result], ignore_index=True)
        eigvecs_list.append(eigvecs)
    results.to_csv(path_or_buf='delphi(sigma_z={:.1e},plane={:},Ib={:}).csv'.format(
        sigma_z, plane, Ib), sep='\t')
    np.save('delphi_eigvecs(sigma_z={:.1e},plane={:},Ib={:}).npy'.format(
        sigma_z, plane, Ib), np.array(eigvecs_list))
    return results


if __name__ == '__main__':
    parser = get_parser_for_delphi()
    args = parser.parse_args()
    if args.scan_type == 'sb_chromaticity':
        run_chroma_scan(impedance_filename=args.filename,
                        ID_state=args.ID_state,
                        sigma_z=args.sigma_z,
                        plane=args.plane,
                        Ib=args.current,
                        min_value=args.min_value,
                        max_value=args.max_value,
                        n_scan_points=args.n_scan_points)
    elif args.scan_type == 'sb_current':
        run_bunch_current_scan(impedance_filename=args.filename,
                               ID_state=args.ID_state,
                               sigma_z=args.sigma_z,
                               plane=args.plane,
                               Qp=args.chromaticity,
                               min_value=args.min_value,
                               max_value=args.max_value,
                               n_scan_points=args.n_scan_points)
