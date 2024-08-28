from utils import get_command_string, execute_commands

# scan_fb_chroma_uniform.py
if __name__ == '__main__':
    config = {
        'plane': 'vertical',
        'min_val': 0.0,
        'max_val': 3.0,
        'n_scan_points': 7,
        'Q_s': 2.0e-3,
        'sub_mode': 'ccrt',
        'damper_gains': [0],
        'currents': [0.3e-3, 0.6e-3, 1.2e-3],
        'id_states': ['close'],
        'impedance_multipliers': [1],
        'scan_types': ['mb_chromaticity']
    }
    execute_commands(config)