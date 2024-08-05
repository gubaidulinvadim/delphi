import os
def get_command_string(id_state,
                        plane, 
                        scan_type,
                        current, 
                        min_val,
                        max_val,
                        n_scan_points,
                        Q_s=2.0e-3,
                        damper_gain=0,
                        sub_mode='slurm'):
    job_name = f'id={id_state}_plane={plane}_current={current}_scan={scan_type}'
    s = 'python submission.py '+f'--damper_gain {damper_gain} '+\
        f'--ID_state {id_state:} ' + f'--filename Zydip_ID{id_state:}.dat '+\
            f'--plane {plane:} '+f'--scan_type {scan_type} '+\
                f'--current {current:} ' + f'--min_value {min_val:} '+\
                    f'--max_value {max_val:} '+f'--n_scan_points {n_scan_points:} '+\
                        f'--Q_s {Q_s} '+f'--job_name {job_name:} '+\
                            f'--sub_mode {sub_mode:} '
    return s
if __name__ == '__main__':
    # os.system('python submission.py --damper_gain 0 --ID_state close --filename Zxdip_IDclose.dat --sigmas_filename sigmas_HC.npy --plane horizontal --scan_type sb_chromaticity --current 1.2 --min_value 0.1 --max_value 5 --n_scan_points 50 --Q_s 0.7e-3 --job_name Chor_chroma416_h --job_time 10000 --sub_mode {mode}')
    mode = 'slurm'
    for damper_gain in [1./100, 1./250, 1./500, 1./1000, 1./2000]:
        for current in [6.25e-3]:
            for id_state in ['open', 'close']:
                s = get_command_string(id_state=id_state,
                                    plane='vertical',
                                    scan_type='sb_chromaticity',
                                    current=current,
                                    min_val=0.1,
                                    max_val=5.0,
                                    n_scan_points=25,
                                    damper_gain=damper_gain,
                                    sub_mode=mode)
                os.system(s + f'--sigma_z 14e-12')
                # os.system(s + f'--sigmas_filename sigmas_Zlong.npy')
    
