import argparse


def get_parser_for_delphi():
    parser = argparse.ArgumentParser(
        description="A script to run DELPHI, a linearized Vlasov solver. Please cite the following if you are using the code\n"
        + "Mounet, N. (2020). Landau damping in the transverse plane. In CERN Yellow Reports: Conference Proceedings (Issue 9).\n"
        + "https://doi.org/https://doi.org/10.23732/CYRCP-2020-009.45\n"
    )
    parser.add_argument(
        "--filename",
        action="store",
        metavar="IMPEDANCE_FILE",
        type=str,
        default="Zydip.dat",
        help='Impedance filename. Three columns (frequencies, real part and imag part) file with tabs as separators. Defaults to "Zydip.dat".',
    )
    parser.add_argument(
        "--id_state",
        action="store",
        metavar="ID_STATE",
        type=str,
        default="open",
        help="State of the IDs. Default to open.",
    )
    parser.add_argument(
        "--sigma_z",
        action="store",
        metavar="SIGMA_Z",
        type=float,
        default=9e-12,
        help="Bunch length (rms) in seconds. Defaults to 9 ps.",
    )
    parser.add_argument(
        "--plane",
        action="store",
        metavar="PLANE",
        type=str,
        default="vertical",
        help='Plane to calculate tune shifts. Defaults to "vertical".',
    )
    parser.add_argument(
        "--scan_type",
        action="store",
        metavar="SCAN_TYPE",
        type=str,
        default="sb_current",
        help='Parameter scan type. Available scans are ["sb_chromaticity", "sb_current", "mb_chromaticity", "mb_current"]. Defaults to "sb_current".',
    )
    parser.add_argument(
        "--current",
        action="store",
        type=float,
        metavar="I_b",
        default=1.2e-3,
        help="Bunch current for chromaticity scan. Defaults to 1.2 mA.",
    )
    parser.add_argument(
        "--chromaticity",
        action="store",
        type=float,
        metavar="Q_p",
        default=0.0,
        help="Chromaticity value for bunch current scan. Default to 0.0 for TMCI scan.",
    )
    parser.add_argument(
        "--max_value",
        action="store",
        metavar="MAX_VALUE",
        type=float,
        default=5.0,
        help="Maximal value for parameter scan. Defaults to 3.0.",
    )
    parser.add_argument(
        "--min_value",
        action="store",
        metavar="MIN_VALUE",
        type=float,
        default=0.1,
        help="Maximal value for parameter scan. Defaults to 0.01.",
    )
    parser.add_argument(
        "--n_scan_points",
        action="store",
        metavar="N_SCAN_POINTS",
        type=int,
        default=50,
        help="Number of values to be scanned. Defaults to 50.",
    )
    parser.add_argument(
        "--Q_s",
        action="store",
        metavar="Q_S",
        type=float,
        default=2.1e-3,
        help="Synchrotron tune. Defaults to 0.0021.",
    )
    parser.add_argument(
        "--sigmas_filename",
        action="store",
        metavar="SIGMAS_FILENAME",
        type=str,
        default=None,
        help='A numpy file with sigmas for different bunch current. The array size should match the scan size. If "None" then sigma_z is used for all currents.Default to None. ',
    )
    parser.add_argument(
        "--n_max",
        action="store",
        type=int,
        default=0,
        help="Coupled bunch mode number. Defaults to 0. A scan will be performed for all modes n < n_max+1.",
    )
    parser.add_argument(
        "--M",
        action="store",
        type=int,
        default=416,
        help="Number of bunches in the beam. Defaults to 416.",
    )

    parser.add_argument(
        "--damper_gain",
        action="store",
        type=float,
        default=0,
        help="Damping time in turns. Defaults to 0",
    )

    parser.add_argument(
        "--impedance_multiplier",
        action="store",
        type=float,
        default=1.0,
        help="Impedance scaling factor. Defaults to 1",
    )
    return parser
