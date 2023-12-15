from DELPHI import *
from SOLEILII_parameters.SOLEILII_TDR_parameters import *
from scipy.constants import pi
import pandas as pd
import h5py as hp
from machine_data.soleil import v2366
from tqdm import tqdm
from utils import get_parser_for_delphi


def read_impedance(filename):
    f = pd.read_csv(filename, sep='\t', header=None)
    freqs = f[0][f[0] > 0]
    Zre = f[1][f[0] > 0]
    Zim = f[2][f[0] > 0]
    Z = np.vstack((Zre, Zim)).T
    return freqs, Z


def run_bunch_current_scan(impedance_filename='/home/dockeruser/delphi/Zydip.dat',
                           ID_state='open',
                           sigma_z=9e-12,
                           plane='vertical',
                           Qp=0.0,
                           min_value=0.1,
                           max_value=5.0,
                           n_scan_points=50,
                           M=1,
                           Q_s=0.0021):
    tune = ring.tune[0] if plane=='horizontal' else ring.tune[1]
    freqs, Z = read_impedance(impedance_filename)
    ring = v2366(IDs=ID_state, load_lattice=False)

    taub = 4*sigma_z
    g0, a, b = longdistribution_decomp(taub, typelong='Gaussian')
    Ib = 1e-3*np.linspace(min_value, max_value, n_scan_points)
    omega_ksi = Qp/ring.ac*ring.omega0
    kmax = 15
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
                                          M=M,
                                          Nb=Nb,
                                          gamma=ring.gamma,
                                          Q=ring.tune[1],
                                          particle='electron')
        eigvals, eigvecs, lmaxold, nmaxold, max_freq_old, matdamperold, matZold = eigenmodesDELPHI_converged(
            nx=M-1,
            M=M,
            omegaksi=omega_ksi,
            omega0=ring.omega0,
            tunefrac=(tune-np.floor(tune)),
            a=a,
            b=b,
            taub=taub,
            g=g0,
            Z=Z,
            kmax=kmax,
            freqZ=freqs,
            coefdamper=coefdamper,
            coefZ=coefZ,
            omegas=Q_s*ring.omega0,
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
        eigvecs_list.append(eigvecs)
    results.to_csv(path_or_buf='delphi(sigma_z={:.1e},plane={:},Qp={:},M={:},Q_s={:.1e}).csv'.format(
        sigma_z, plane, Qp, M, Q_s), sep='\t')
    # np.save('delphi_eigvecs(sigma_z={:.1e},plane={:},Qp={:}).npy'.format(
    # sigma_z, plane, Qp), np.array(eigvecs_list, dtype=object))
    return results


def run_chroma_scan(impedance_filename='/home/dockeruser/delphi/Zydip.dat',
                    ID_state='open',
                    sigma_z=9e-12,
                    plane='vertical',
                    Ib=1.2e-3,
                    min_value=0.1,
                    max_value=5.0,
                    n_scan_points=50,
                    M=1,
                    Q_s=0.0021):
    tune = ring.tune[0] if plane=='horizontal' else ring.tune[1]
    freqs, Z = read_impedance(impedance_filename)
    ring = v2366(IDs=ID_state, load_lattice=False)
    taub = 4*sigma_z
    g0, a, b = longdistribution_decomp(taub, typelong='Gaussian')
    kmax = 15
    results = pd.DataFrame(columns=['BunchCurrent',
                                    'eigvals_re',
                                    'eigvals_im',
                                    'Qp',
                                    'BunchLength'])
    results_h5py = hp.File('delphi_resultssigma_z={:.1e},plane={:},Ib={:}).h5'.format(
        sigma_z, plane, Ib), 'w')

    bunch_current_ds = results_h5py.create_dataset(
        'BunchCurrent', (n_scan_points,), dtype='f')
    eigvals_re_ds = results_h5py.create_dataset(
        'eigvals_re', (n_scan_points, kmax), dtype='f')
    eigvals_im_ds = results_h5py.create_dataset(
        'eigvals_im', (n_scan_points, kmax), dtype='f')
    Qp_ds = results_h5py.create_dataset('Qp', (n_scan_points,), dtype='f')
    bunch_length_ds = results_h5py.create_dataset(
        'BunchLength', (n_scan_points,), dtype='f')
    eigvecs_group = results_h5py.create_group('Eigenvectors')

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
                                          Q=tune,
                                          particle='electron')
        eigvals, eigvecs, lmaxold, nmaxold, max_freq_old, matdamperold, matZold = eigenmodesDELPHI_converged(
            nx=M-1,
            M=M,
            omegaksi=omega_ksi,
            omega0=ring.omega0,
            tunefrac=(tune-np.floor(tune)),
            a=a,
            b=b,
            taub=taub,
            g=g0,
            Z=Z,
            kmax=kmax,
            freqZ=freqs,
            coefdamper=coefdamper,
            coefZ=coefZ,
            omegas=Q_s*ring.omega0,
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
    results.to_csv(path_or_buf='delphi(sigma_z={:.1e},plane={:},Ib={:},M={:},Q_s={:.1e}).csv'.format(
        sigma_z, plane, Ib, M, Q_s), sep='\t')
    # np.save('delphi_eigvecs(sigma_z={:.1e},plane={:},Ib={:}).npy'.format(
    # sigma_z, plane, Ib), np.array(eigvecs_list, dtype=object))
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
                        n_scan_points=args.n_scan_points,
                        Q_s=args.Q_s)
    elif args.scan_type == 'sb_current':
        print("running a single bunch current scan")
        run_bunch_current_scan(impedance_filename=args.filename,
                               ID_state=args.ID_state,
                               sigma_z=args.sigma_z,
                               plane=args.plane,
                               Qp=args.chromaticity,
                               min_value=args.min_value,
                               max_value=args.max_value,
                               n_scan_points=args.n_scan_points,
                               Q_s=args.Q_s)
    elif args.scan_type == 'mb_chromaticity':
        run_chroma_scan(impedance_filename=args.filename,
                        ID_state=args.ID_state,
                        sigma_z=args.sigma_z,
                        plane=args.plane,
                        Ib=args.current,
                        min_value=args.min_value,
                        max_value=args.max_value,
                        n_scan_points=args.n_scan_points,
                        M=416,
                        Q_s=args.Q_s)
    elif args.scan_type == 'mb_current':
        run_bunch_current_scan(impedance_filename=args.filename,
                               ID_state=args.ID_state,
                               sigma_z=args.sigma_z,
                               plane=args.plane,
                               Qp=args.chromaticity,
                               min_value=args.min_value,
                               max_value=args.max_value,
                               n_scan_points=args.n_scan_points,
                               M=416,
                               Q_s=args.Q_s)
    else:
        print('Scan type is not supported.')
