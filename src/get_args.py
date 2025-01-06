import argparse
import configparser
import os
from astropy.time import Time


def get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ra",
        type=float,
        default=0.0,
        help="Telescope center RA [degrees] (default: 0.0)",
    )
    parser.add_argument(
        "--dec",
        type=float,
        default=0.0,
        help="Telescope center Dec [degrees] (default: 0.0)",
    )
    parser.add_argument(
        "--pa",
        type=float,
        default=0.0,
        help="Telescope position angle [degrees] (default: 0.0)",
    )
    parser.add_argument(
        "--observation_time",
        type=str,
        default="2021-11-20T15:00:00Z",
        help="planned time of observation in UTC (default: 2021-11-20 15:00:00)",
    )
    parser.add_argument(
        "--lim_target_mag",
        type=float,
        default="19.",
        help="magnitude of the faintest targets (obsolete) (default:19)",
    )

    parser.add_argument(
        "--design_dir",
        type=str,
        default=".",
        help="directory for storing PFS designs (default: .)",
    )

    parser.add_argument(
        "--guidestar_mag_max",
        type=float,
        default=19.0,
        help="maximum magnitude for guide star candidates (default: 19.)",
    )
    parser.add_argument(
        "--guidestar_neighbor_mag_min",
        type=float,
        default=21.0,
        help="minimum magnitude for objects in the vicinity of guide star candidates (default: 21.)",
    )
    parser.add_argument(
        "--guidestar_minsep_deg",
        type=float,
        default=1.0 / 3600,
        help="radius of guide star candidate vicinity (default: 1/3600)",
    )

    parser.add_argument(
        "--use_gurobi",
        action="store_true",
        help="use Gurobi",
    )
    parser.add_argument(
        "--cobra_coach_dir",
        type=str,
        default=".",
        help="path for temporary cobraCoach files (default: .)",
    )

    parser.add_argument(
        "--cobra_coach_module_version",
        type=str,
        default="final_20210920_mm",
        help="version of the bench decription file (default: final_20210920_mm)",
    )

    parser.add_argument(
        "--targetdb_conf",
        type=str,
        default="targetdb_config.ini",
        help="Config file for targetDB (default: targetdb_config.ini)",
    )
    parser.add_argument(
        "--gaiadb_conf",
        type=str,
        default="gaiadb_config_hilo.ini",
        help="Config file for Subaru's Gaia DB (default: gaiadb_config_hilo.ini",
    )
    parser.add_argument(
        "--target_mag_max",
        type=float,
        default=19.0,
        help="Maximum (faintest) magnitude for stars in fibers (default: 19.)",
    )
    parser.add_argument(
        "--target_mag_min",
        type=float,
        default=0.0,
        help="Minimum (brightest) magnitude for stars in fibers (default: 0)",
    )
    parser.add_argument(
        "--target_mag_filter",
        type=str,
        default="g",
        help="Photometric band (grizyj of PS1) to apply magnitude cuts (default: g)",
    )
    parser.add_argument(
        "--fluxstd_min_prob_f_star",
        type=float,
        default=0.0,
        help="Minimum acceptable prob_f_star (default: 0)",
    )
    parser.add_argument(
        "--telescope_elevation",
        type=float,
        default=60.0,
        help="Telescope elevation in degree (default: 60)",
    )
    parser.add_argument(
        "--n_fluxstd",
        type=int,
        default=50,
        help="Number of FLUXSTD stars to be allocated. (default: 50)",
    )
    parser.add_argument(
        "--pfs_instdata_dir",
        type=str,
        default="/Users/monodera/Dropbox/NAOJ/PFS/Subaru-PFS/pfs_instdata/",
        help="Location of pfs_instdata (default: /Users/monodera/Dropbox/NAOJ/PFS/Subaru-PFS/pfs_instdata/)",
    )

    args = parser.parse_args()

    if args.observation_time.lower() == "now":
        print("converting to the current time")
        args.observation_time = (
            Time.now().iso
        )  # astropy.time.Time.now() uses datetime.utcnow()

    return args