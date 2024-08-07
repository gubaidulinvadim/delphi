import os
if __name__ == '__main__':
    Q_s = 2e-3
    os.system(f'python submission.py --ID_state open --filename Zydip_IDopen.dat --sigma_z 8.52e-12 --plane vertical --scan_type sb_chromaticity --current 6.25 --min_value 0.1 --max_value 5 --n_scan_points 50 --Q_s 2.117e-3 --job_name Over_chroma32_b --job_time 10000 --sub_mode slurm')
    os.system(f'python submission.py --ID_state open --filename Zydip_IDopen.dat --sigmas_filename sigmas_Zlong.npy --plane vertical --scan_type sb_chromaticity --current 6.25 --min_value 0.1 --max_value 5 --n_scan_points 50 --Q_s 2.117e-3 --job_name Over_chroma32_z --job_time 10000 --sub_mode slurm')
    os.system(f'python submission.py --ID_state open --filename Zydip_IDopen.dat --sigmas_filename sigmas_HC.npy --plane vertical --scan_type sb_chromaticity --current 6.25 --min_value 0.1 --max_value 5 --n_scan_points 50 --Q_s 0.7e-3 --job_name Over_chroma32_h --job_time 10000 --sub_mode slurm')

    os.system(f'python submission.py --ID_state close --filename Zydip_IDclose.dat --sigma_z 8.52e-12 --plane vertical --scan_type sb_chromaticity --current 6.25 --min_value 0.1 --max_value 5 --n_scan_points 50 --Q_s 2.117e-3 --job_name Cver_chroma32_b --job_time 10000 --sub_mode slurm')
    os.system(f'python submission.py --ID_state close --filename Zydip_IDclose.dat --sigmas_filename sigmas_Zlong.npy --plane vertical --scan_type sb_chromaticity --current 6.25 --min_value 0.1 --max_value 5 --n_scan_points 50 --Q_s 2.117e-3 --job_name Cver_chroma32_z --job_time 10000 --sub_mode slurm')
    os.system(f'python submission.py --ID_state close --filename Zydip_IDclose.dat --sigmas_filename sigmas_HC.npy --plane vertical --scan_type sb_chromaticity --current 6.25 --min_value 0.1 --max_value 5 --n_scan_points 50 --Q_s 0.7e-3 --job_name Cver_chroma32_h --job_time 10000 --sub_mode slurm')

    os.system(f'python submission.py --ID_state close --filename Zxdip_IDclose.dat --sigma_z 8.52e-12 --plane horizontal --scan_type sb_chromaticity --current 6.25 --min_value 0.1 --max_value 5 --n_scan_points 50 --Q_s 2e-3 --job_name Chor_chroma32_b --job_time 10000 --sub_mode slurm')
    os.system(f'python submission.py --ID_state close --filename Zxdip_IDclose.dat --sigmas_filename sigmas_Zlong.npy --plane horizontal --scan_type sb_chromaticity --current 6.25 --min_value 0.1 --max_value 5 --n_scan_points 50 --Q_s 2.117e-3 --job_name Chor_chroma32_z --job_time 10000 --sub_mode slurm')
    os.system(f'python submission.py --ID_state close --filename Zxdip_IDclose.dat --sigmas_filename sigmas_HC.npy --plane horizontal --scan_type sb_chromaticity --current 6.25 --min_value 0.1 --max_value 5 --n_scan_points 50 --Q_s 0.7e-3 --job_name Chor_chroma32_h --job_time 10000 --sub_mode slurm')
