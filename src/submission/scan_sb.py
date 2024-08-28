from utils import get_command_string, execute_commands

if __name__ == '__main__':
    config = {
        'plane': 'vertical',
        'min_val': 0.0,
        'max_val': 7.0,
        'n_scan_points': 71,
        'Q_s': 2.0e-3,
        'sub_mode': 'ccrt',
        'damper_gains': [0, 1/50, 1/100, 1/200],
        'currents': [0],  # Add appropriate values if needed
        'id_states': ['close'],
        'impedance_multipliers': [1, 2],
        'scan_types': ['sb_current']
    }
    execute_commands(config)