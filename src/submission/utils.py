# common.py
import os
import itertools

def get_command_string(config):
    job_name = f'id={config["id_state"]}_plane={config["plane"]}_current={config["current"]}_ximp={config["impedance_multiplier"]:.1f},damper_gain={config["damper_gain"]:.1e}'
    return (
        'python submission.py '
        f'--damper_gain {config["damper_gain"]:.1e} '
        f'--id_state {config["id_state"]} '
        f'--filename Zydip_ID{config["id_state"]}.dat '
        f'--plane {config["plane"]} '
        f'--scan_type {config["scan_type"]} '
        f'--current {config["current"]} '
        f'--min_value {config["min_val"]} '
        f'--max_value {config["max_val"]} '
        f'--n_scan_points {config["n_scan_points"]} '
        f'--Q_s {config["Q_s"]} '
        f'--job_name {job_name} '
        f'--n_max 415 '
        f'--impedance_multiplier {config["impedance_multiplier"]:.1f} '
        f'--sub_mode {config["sub_mode"]} '
    )

def execute_commands(config):
    for damper_gain, current, id_state, ximp, scan_type in itertools.product(
            config['damper_gains'], 
            config['currents'], 
            config['id_states'], 
            config['impedance_multipliers'], 
            config['scan_types']):
        
        config.update({
            'damper_gain': damper_gain,
            'current': current,
            'id_state': id_state,
            'impedance_multiplier': ximp,
            'scan_type': scan_type
        })
        
        command = get_command_string(config)
        os.system(command + ' --sigma_z 9e-12')
        os.system(command + ' --sigmas_filename sigmas_zlong.txt')
        os.system(command + ' --sigmas_filename sigmas_hc.txt')