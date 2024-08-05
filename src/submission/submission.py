import os, sys

sys.path.append('../')
from simulation.utils import get_parser_for_delphi


def write_tmp_submission_script_ccrt(
    job_name,
    job_time,
    filename,
    id_state,
    sigma_z,
    plane,
    scan_type,
    Ib,
    Qp,
    min_value,
    max_value,
    n_scan_points,
    Q_s,
    sigmas_filename,
    n_max,
    damper_gain,
    impedance_multiplier
):
    MOUNT_FOLDER = "/ccc/work/cont003/soleil/gubaiduv/delphi"
    dst_folder = "/home/dockeruser/delphi"
    IMAGE_NAME = "delphi"
    SCRIPT_NAME = "/home/dockeruser/delphi/src/simulation/run_delphi.py"
    with open(job_name, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("#MSUB -m work,scratch\n")
        f.write("#MSUB -q milan\n")
        f.write("#MSUB -Q long\n")
        f.write("#MSUB -n 1\n")
        f.write("#MSUB -c 4\n")
        f.write("#MSUB -T {:}\n".format(job_time))
        f.write("#MSUB -A soleil\n")
        f.write("#MSUB -@ gubaidulinvadim@gmail.com:begin,end,requeue\n")
        f.write("#MSUB -o /ccc/cont003/home/soleil/gubaiduv/{0:}.err\n".format(
            job_name))
        f.write("#MSUB -e /ccc/cont003/home/soleil/gubaiduv/{0:}.out\n".format(
            job_name))
        f.write("module purge\n")
        f.write(
            f"ccc_mprun -C {IMAGE_NAME} -E '--ctr-mount src={MOUNT_FOLDER:},dst={dst_folder:}' -- python {SCRIPT_NAME:} --filename {filename:} --id_state {id_state:} --sigma_z {sigma_z:} --plane {plane:} --scan_type {scan_type:} --current {Ib:} --chromaticity {Qp:} --min_value {min_value:} --max_value {max_value:} --n_scan_points {n_scan_points:} --Q_s {Q_s:} --sigmas_filename {sigmas_filename:} --n_max {n_max:} --damper_gain {damper_gain:} --impedance_multiplier {impedance_multiplier:}\n"
            )


def write_submission_script_slurm(
    job_name,
    job_time,
    filename,
    id_state,
    sigma_z,
    plane,
    scan_type,
    Ib,
    Qp,
    min_value,
    max_value,
    n_scan_points,
    Q_s,
    sigmas_filename,
    n_max,
    damper_gain,
    impedance_multiplier
):
    MOUNT_FOLDER = (
        "/lustre/scratch/sources/physmach/gubaidulin/delphi:/home/dockeruser/delphi"
    )
    IMAGE_NAME = "/lustre/scratch/sources/physmach/gubaidulin/delphi.sif"
    SCRIPT_NAME = "/home/dockeruser/delphi/src/simulation/run_delphi.py"
    with open(job_name, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --partition sumo\n")
        f.write("#SBATCH -n 4\n")
        f.write("#SBATCH --time={:}\n".format(job_time))
        f.write("#SBATCH --export=ALL\n")
        f.write("#SBATCH --mail-user='gubaidulinvadim@gmail.com'\n")
        f.write("#SBATCH --mail-type=begin,end,requeue\n")
        f.write(
            "#SBATCH --error=/home/sources/physmach/gubaidulin/err/{0:}.err\n".
            format(job_name))
        f.write(
            "#SBATCH --output=/home/sources/physmach/gubaidulin/out/{0:}.err\n"
            .format(job_name))
        f.write("module load tools/singularity/current\n")
        f.write(
            "singularity exec -e --no-home -B {0:} {1:} python {2:}  --filename {3:} --id_state {4:} --sigma_z {5:} --plane {6:} --scan_type {7:} --current {8:} --chromaticity {9:} --min_value {10:} --max_value {11:} --n_scan_points {12:} --Q_s {13:} --sigmas_filename {14:} --n_max {15:} --damper_gain {16:} --impedance_multiplier {17:}\n"
            .format(
                MOUNT_FOLDER,
                IMAGE_NAME,
                SCRIPT_NAME,
                filename,
                id_state,
                sigma_z,
                plane,
                scan_type,
                Ib,
                Qp,
                min_value,
                max_value,
                n_scan_points,
                Q_s,
                sigmas_filename,
                n_max,
                damper_gain,
                impedance_multiplier
            ))
    return job_name


if __name__ == "__main__":
    parser = get_parser_for_delphi()
    parser.add_argument(
        "--job_name",
        action="store",
        metavar="JOB_NAME",
        type=str,
        default="job",
        help=
        'Name of the job and associated .our and .err files. Defaults to "job".',
    )
    parser.add_argument(
        "--job_time",
        action="store",
        metavar="JOB_TIME",
        type=int,
        default=259_100,
        help="Time allocated to the job. Defaults to 10000.",
    )
    parser.add_argument(
        "--sub_mode",
        action="store",
        metavar="SUB_MODE",
        type=str,
        default="ccrt",
        help=
        'Submission mode. Accepted values are ["local", "ccrt", "slurm"], defaults to "ccrt"',
    )
    args = parser.parse_args()
    if args.sub_mode == "ccrt":
        write_tmp_submission_script_ccrt(
            args.job_name,
            args.job_time,
            args.filename,
            args.id_state,
            args.sigma_z,
            args.plane,
            args.scan_type,
            args.current,
            args.chromaticity,
            args.min_value,
            args.max_value,
            args.n_scan_points,
            args.Q_s,
            args.sigmas_filename,
            args.n_max,
            args.damper_gain,
            args.impedance_multiplier
        )
        os.system("ccc_msub {:}".format(args.job_name))
        os.system("rm -rf {:}".format(args.job_name))

    elif args.sub_mode == "slurm":
        write_submission_script_slurm(
            args.job_name,
            args.job_time,
            args.filename,
            args.id_state,
            args.sigma_z,
            args.plane,
            args.scan_type,
            args.current,
            args.chromaticity,
            args.min_value,
            args.max_value,
            args.n_scan_points,
            args.Q_s,
            args.sigmas_filename,
            args.n_max,
            args.damper_gain,
            args.impedance_multiplier
        )
        os.system("sbatch {:}".format(args.job_name))
        os.system("rm -rf {:}".format(args.job_name))
    elif args.sub_mode == "local":
        pass
