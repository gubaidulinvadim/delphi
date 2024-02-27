import os
import numpy as np

if __name__ == "__main__":
    min_current = 0
    max_current = 6.25
    n_scan_points = 101
    sub_mode = "slurm"
    chromaticity = 1.6
    sigma_z1 = 8.52e-12
    n_max = 31
    M = 32
    plane = "vertical"
    # os.system(
    #     f"python submission.py --ID_state {id_state:}"
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
        os.system(
            f"python submission.py --ID_state open --filename Zydip_IDopen.dat --sigma_z 8.52e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:} --max_value {max_current:} --n_scan_points {n_scan_points} --Q_s 2.117e-3 --job_name Over_tcbi_nominal_b --job_time 10000 --sub_mode slurm  --n_max {n_max:} --M {M:}"
        )
        os.system(
            f"python submission.py --ID_state open --filename Zydip_IDopen.dat --sigma_z 14e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:}.0 --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 2.117e-3 --job_name Over_tcbi_nominal_z --job_time 10000 --sub_mode slurm  --n_max {n_max:} --M {M:}"
        )
        os.system(
            f"python submission.py --ID_state open --filename Zydip_IDopen.dat --sigma_z 40e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:}.0 --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 0.7e-3 --job_name Over_tcbi_nominal_h --job_time 10000 --sub_mode slurm  --n_max {n_max:} --M {M:}"
        )

        os.system(
            f"python submission.py --ID_state close --filename Zydip_IDclose.dat --sigma_z 8.52e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:} --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 2.117e-3 --job_name Cver_tcbi_nominal_b --job_time 10000 --sub_mode slurm  --n_max {n_max:} --M {M:}"
        )
        os.system(
            f"python submission.py --ID_state close --filename Zydip_IDclose.dat --sigma_z 14e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:} --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 2.117e-3 --job_name Cver_tcbi_nominal_z --job_time 10000 --sub_mode slurm  --n_max {n_max:} --M {M:}"
        )
        os.system(
            f"python submission.py --ID_state close --filename Zydip_IDclose.dat --sigma_z 40e-12 --plane vertical --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:} --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 0.7e-3 --job_name Cver_tcbi_nominal_h --job_time 10000 --sub_mode slurm  --n_max {n_max:} --M {M:}"
        )

        os.system(
            f"python submission.py --ID_state close --filename Zxdip_IDclose.dat --sigma_z 8.52e-12 --plane horizontal --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:}.0 --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 2.117e-3 --job_name Chor_tcbi_nominal_b --job_time 10000 --sub_mode slurm  --n_max {n_max:} --M {M:}"
        )
        os.system(
            f"python submission.py --ID_state close --filename Zxdip_IDclose.dat --sigma_z 14e-12 --plane horizontal --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:}.0 --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 2.117e-3 --job_name Chor_tcbi_nominal_z --job_time 10000 --sub_mode slurm  --n_max {n_max} --M {M:}"
        )
        os.system(
            f"python submission.py --ID_state close --filename Zxdip_IDclose.dat --sigma_z 40e-12 --plane horizontal --scan_type mb_current --chromaticity {chromaticity:} --min_value {min_current:}.0 --max_value {max_current:} --n_scan_points {n_scan_points:} --Q_s 0.7e-3 --job_name Chor_tcbi_nominal_h --job_time 10000 --sub_mode slurm  --n_max {n_max} --M {M:}"
        )
