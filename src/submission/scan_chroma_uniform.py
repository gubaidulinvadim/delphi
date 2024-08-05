import os
def get_command_string(id_state,
                        plane, 
                        scan_type,
                        current, 
                        min_val,
                        max_val,
                        n_scan_points,
                        Q_s=2.117e-3,
                        damper_gain=0,
                        impedance_multiplier=1,
                        sub_mode='slurm'):
    job_name = f'id={id_state}_plane={plane}_current={current}_ximp={impedance_multiplier:.1f},damper_gain={damper_gain:.1e}_sb'
    s = 'python submission.py '+f'--damper_gain {damper_gain:.1e} '+\
        f'--id_state {id_state:} ' + f'--filename Zydip_ID{id_state:}.dat '+\
            f'--plane {plane:} '+f'--scan_type {scan_type} '+\
                f'--current {current:} ' + f'--min_value {min_val:} '+\
                    f'--max_value {max_val:} '+f'--n_scan_points {n_scan_points:} '+\
                        f'--Q_s {Q_s} '+f'--job_name {job_name:} '+\
                            f'--n_max 0 '+\
                            f'--impedance_multiplier {impedance_multiplier:1f} '+\
                            f'--sub_mode {sub_mode:} '
    return s
if __name__ == '__main__':
    mode = 'ccrt'
    for damper_gain in [0, 1/50, 1/100, 1/200, 1/400]:
        for current in np.linspace(0, 5, 26):
            for id_state in ['close']:
                for ximp in [1, 2]:
                    s = get_command_string(id_state=id_state,
                                        plane='vertical',
                                        scan_type='sb_chromaticity',
                                        chromaticity=chromaticity,
                                        min_val=0.1,
                                        max_val=1.2,
                                        n_scan_points=101,
                                        damper_gain=damper_gain,
                                        impedance_multiplier=ximp,
                                        sub_mode=mode)
                    os.system(s + f'--sigma_z 9e-12')
                    os.system(s + f'--sigmas_filename sigmas_zlong.txt')
    
