# Fiber Assignment for PFS SSP Cosmology Targets

This repository contains the code for performing fiber assignment for the Prime Focus Spectrograph (PFS) Subaru Strategic Program (SSP) cosmology targets.

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Introduction
The goal of this project is to efficiently assign fibers to cosmology targets for the PFS SSP. This involves optimizing the assignment process to maximize the number of targets observed while adhering to various constraints.

### Workflow
1. **Data Preparation**: Collect and preprocess the list of cosmology targets, calibration targets, and ancillary targets into .ecsv format file.
   + ID (unique)
   + R.A. [deg.]
   + Dec. [deg.]
   + Exposure Time [sec.]
   + Priority [1(highest) - 15(lowest)]
   + Magnitude [ABmag]. (optional)
   + Redshift (optional)
   + Object Type (optional)
   + stage (optional, default=0 for cosmology+calibration targets, =1 for ancillary targets)
     (above is the requirement by netflow, below is what we should add for SSP run)
   + AlreadyObserved (default=0, i.e. False; =1, True)
   + PartiallyObserved (???)
   + otime (otime value used when running netflow, overwritten by the last fiber design run before it's observed)
   + TileID (default=-1, set this when AlreadyObserved is set to True)
     
   + TileIDExtraExposure (TileID for the extra exposure in the 2nd pass)
   + otimeExtraExposure (otime used if extra exposure is assigned to the target in the 2nd pass)
   + BenchFile (bench file name used when running netflow, remove if not updated often)
   + BenchFileExtraExposure (bench file for the extra exposure in the 2nd pass, remove if not updated often)
2. **Predefine Pointing Centers (i.e. tiling strategy)**: Predefine the pointing centers for all visits
   + raTel [deg.]
   + decTel [deg.]
   + TileID [Unique value in the whole survey]
6. **Tile Selection**: Select next tile(s) 
5. **Verfication of inputs**: Check the updates of bench file (focal plane status), calibration targets, otime, TileID
7. **Fiber Assignment**: Use netflow to assign fibers to the targets in the slected tile(s), generate pfsDesign file
8. **Update FA Log**: update the fiber assignment note file
9. **Update Target List**: Update target list columns after observation
10. **Verify the Reproducibility**: verify the reproducibility the next day

***Track/record the status of the instrument through time (mainly bench files)***

## Installation
**NOTE** CollisionSimulator2 requires installation of opdb from https://github.com/Subaru-PFS/spt_operational_database
To install the necessary dependencies, run the following command:
(need to update the requirements.txt)
```bash
pip install -r requirements.txt
```

git clone https://github.com/jingssrs/pfs_co_fa

## Usage
To perform fiber assignment, use the following command (??? revise later):
```bash
```

## Contributing
We welcome contributions to improve the fiber assignment algorithm. Please fork the repository and submit a pull request with your changes.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
