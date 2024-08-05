import os
import numpy as np

def get_command_string(id_state,
                        plane, 
                        scan_type,
                        current, 
                        min_val,
                        max_val,
                        n_scan_points,
                        Q_s=2.0e-3,
                        damper_gain=0,
                        impedance_multiplier=1,
                        sub_mode='slurm'):
    job_name = f'id={id_state}_plane={plane}_current={current}_ximp={impedance_multiplier:.1f},damper_gain={damper_gain:.1e}'
    s = 'python submission.py '+f'--damper_gain {damper_gain:.1e} '+\
        f'--id_state {id_state:} ' + f'--filename Zydip_ID{id_state:}.dat '+\
            f'--plane {plane:} '+f'--scan_type {scan_type} '+\
                f'--current {current:} ' + f'--min_value {min_val:} '+\
                    f'--max_value {max_val:} '+f'--n_scan_points {n_scan_points:} '+\
                        f'--Q_s {Q_s} '+f'--job_name {job_name:} '+\
                            f'--n_max 415 '+\
                            f'--impedance_multiplier {impedance_multiplier:1f} '+\
                            f'--sub_mode {sub_mode:} '
    return s
if __name__ == "__main__":
    min_current = 0
    max_current = 1.2
    n_scan_points = 101
    sub_mode = "slurm"
    chromaticity = 1.6
    sigma_z1 = 8.52e-12
    n_max = 415
    plane = "vertical"
    # os.system(
    #     f"python submission.py --id_state {id_state:}"
    #     + f" --filename {filename:}"
    #     + f" --sigma_z {sigma_z1:}"
    #     f" --plane {plane:}"
    #     + f" --chromaticity {chromaticity:}"
    #     + f" --Q_s {Q_s:}"
    #     + f" --min_value {min_value:}"
    #     + f" --max_value {max_value:}"
    #     + f" --n_scan_points {n_scan_points:}"
    #     + f" --n_max {n_max:}"
    #     + f" --job_name {jobname:}"
    #     + f" --job_time 10000"
    #     + f" --sub_mode {sub_mode:}"
    # )
    for chromaticity in np.linspace(0, 5, 26):
        # os.system(
        #     f"python submission.py --id_state open --filename Zydip_IDopen.dat --sigma_z 9e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:} --max_value {max_current:} --n_scan_points {n_scan_points} --Q_s 2.117e-3 --job_name Over_tcbi_nominal_b --job_time 10000 --sub_mode slurm  --n_max 415"
        # )
        # os.system(
        #     f"python submission.py --id_state open --filename Zydip_IDopen.dat --sigmas_filename sigmas_zlong.txt --sigma_z 14e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:}.0 --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 2.117e-3 --job_name Over_tcbi_nominal_z --job_time 10000 --sub_mode slurm  --n_max 415"
        # )
        # # os.system(
        # #     f"python submission.py --id_state open --filename Zydip_IDopen.dat --sigma_z 40e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:}.0 --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 0.7e-3 --job_name Over_tcbi_nominal_h --job_time 10000 --sub_mode slurm  --n_max 415"
        # # )

        os.system(
            f"python submission.py --id_state close --filename Zydip_IDclose.dat --sigma_z 9e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:} --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 2.117e-3 --job_name Cver_tcbi_nominal_b --job_time 10000 --sub_mode slurm  --n_max 415"
        )
        os.system(
            f"python submission.py --id_state close --filename Zydip_IDclose.dat --sigmas_filename sigmas_zlong.txt --sigma_z 14e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:} --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 2.117e-3 --job_name Cver_tcbi_nominal_z --job_time 10000 --sub_mode slurm  --n_max 415"
        )
        os.system(
            f"python submission.py --id_state close --filename Zydip_IDclose.dat --sigma_z 40e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:} --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 0.7e-3 --job_name Cver_tcbi_nominal_h --job_time 10000 --sub_mode slurm  --n_max 415"
        )

        # os.system(
        #     f"python submission.py --id_state close --filename Zxdip_IDclose.dat --sigma_z 9e-12 --plane horizontal --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:}.0 --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 2.117e-3 --job_name Chor_tcbi_nominal_b --job_time 10000 --sub_mode slurm  --n_max 415"
        # )
        # os.system(
        #     f"python submission.py --id_state close --filename Zxdip_IDclose.dat --sigmas_filename sigmas_zlong.txt --sigma_z 14e-12 --plane horizontal --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:}.0 --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 2.117e-3 --job_name Chor_tcbi_nominal_z --job_time 10000 --sub_mode slurm  --n_max 415"
        # )
        # # os.system(
        # #     f"python submission.py --id_state close --filename Zxdip_IDclose.dat --sigma_z 40e-12 --plane horizontal --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:}.0 --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 0.7e-3 --job_name Chor_tcbi_nominal_h --job_time 10000 --sub_mode slurm  --n_max 415"
        # # )
