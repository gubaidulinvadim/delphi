from utils import get_parser_for_delphi
import os


def write_tmp_submission_script_ccrt(
    job_name,
    job_time,
    filename,
    ID_state,
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
):
    MOUNT_FOLDER = "/ccc/work/cont003/soleil/gubaiduv/delphi:/home/dockeruser/delphi"
    IMAGE_NAME = "delphi"
    SCRIPT_NAME = "delphi/run_delphi.py"
    with open(job_name, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("#MSUB -m work,scratch\n")
        f.write("#MSUB -q milan\n")
        f.write("#MSUB -Q normal\n")
        f.write("#MSUB -n 1\n")
        f.write("#MSUB -c 32\n")
        f.write("#MSUB -T {:}\n".format(job_time))
        f.write("#MSUB -A soleil\n")
        f.write("#MSUB -@ gubaidulinvadim@gmail.com:begin,end,requeue\n")
        f.write(
            "#MSUB -o /ccc/cont003/home/soleil/gubaiduv/{0:}.err\n".format(job_name)
        )
        f.write(
            "#MSUB -e /ccc/cont003/home/soleil/gubaiduv/{0:}.out\n".format(job_name)
        )
        f.write("module purge\n")
        f.write(
            "pcocc run --mount {0:} -I {1:} --entry-point -- python {2:} --filename {3:} --ID_state {4:} --sigma_z {5:} --plane {6:} --scan_type {7:} --current {8:} --chromaticity {9:} --min_value {10:} --max_value {11:} --n_scan_points {12:} --Q_s {13:} --sigmas_filename {14:} --n_max {15:}\n".format(
                MOUNT_FOLDER,
                IMAGE_NAME,
                SCRIPT_NAME,
                filename,
                ID_state,
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
            )
        )


def write_submission_script_slurm(
    job_name,
    job_time,
    filename,
    ID_state,
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
):
    MOUNT_FOLDER = (
        "/lustre/scratch/sources/physmach/gubaidulin/delphi:/home/dockeruser/delphi"
    )
    IMAGE_NAME = "/lustre/scratch/sources/physmach/gubaidulin/delphi.sif"
    SCRIPT_NAME = "/home/dockeruser/delphi/run_delphi.py"
    with open(job_name, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --partition sumo\n")
        f.write("#SBATCH -n 8\n")
        f.write("#SBATCH --time={:}\n".format(job_time))
        f.write("#SBATCH --export=ALL\n")
        f.write("#SBATCH --mail-user='gubaidulinvadim@gmail.com'\n")
        f.write("#SBATCH --mail-type=begin,end,requeue\n")
        f.write(
            "#SBATCH --error=/home/sources/physmach/gubaidulin/err/{0:}.err\n".format(
                job_name
            )
        )
        f.write(
            "#SBATCH --output=/home/sources/physmach/gubaidulin/out/{0:}.err\n".format(
                job_name
            )
        )
        f.write("module load tools/singularity/current\n")
        f.write(
            "singularity exec -e --no-home -B {0:} {1:} python {2:}  --filename {3:} --ID_state {4:} --sigma_z {5:} --plane {6:} --scan_type {7:} --current {8:} --chromaticity {9:} --min_value {10:} --max_value {11:} --n_scan_points {12:} --Q_s {13:} --sigmas_filename {14:} --n_max {15:} \n".format(
                MOUNT_FOLDER,
                IMAGE_NAME,
                SCRIPT_NAME,
                filename,
                ID_state,
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
            )
        )
    return job_name


if __name__ == "__main__":
    parser = get_parser_for_delphi()
    parser.add_argument(
        "--job_name",
        action="store",
        metavar="JOB_NAME",
        type=str,
        default="job",
        help='Name of the job and associated .our and .err files. Defaults to "job".',
    )
    parser.add_argument(
        "--job_time",
        action="store",
        metavar="JOB_TIME",
        type=int,
        default=10000,
        help="Time allocated to the job. Defaults to 10000.",
    )
    parser.add_argument(
        "--sub_mode",
        action="store",
        metavar="SUB_MODE",
        type=str,
        default="ccrt",
        help='Submission mode. Accepted values are ["local", "ccrt", "slurm"], defaults to "ccrt"',
    )
    args = parser.parse_args()
    if args.sub_mode == "ccrt":
        write_tmp_submission_script_ccrt(
            args.job_name,
            args.job_time,
            args.filename,
            args.ID_state,
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
        )
        os.system("ccc_msub {:}".format(args.job_name))
        os.system("rm -rf {:}".format(args.job_name))

    elif args.sub_mode == "slurm":
        write_submission_script_slurm(
            args.job_name,
            args.job_time,
            args.filename,
            args.ID_state,
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
        )
        os.system("sbatch {:}".format(args.job_name))
        os.system("rm -rf {:}".format(args.job_name))
    elif args.sub_mode == "local":
        pass
