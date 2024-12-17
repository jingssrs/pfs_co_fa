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
1. **Predefine Pointing Centers (i.e. tiling strategy)**: Predefine the pointing centers for each visit
2. **Data Preparation**: Collect and preprocess the list of cosmology targets.
3. **Optimization**: Use the fiber assignment algorithm to optimize the assignment of fibers to targets.
4. **Output Generation**: Generate the output file with the assigned fibers.
5. **Verfify the Reproducibility**: 

***Track current observational state of the targets and the current health of the instrument through time***


### Packages
- `numpy`: For numerical operations.
- `pandas`: For data manipulation and analysis.
- `scipy`: For optimization algorithms.
- `matplotlib`: For plotting and visualization.

## Installation
To install the necessary dependencies, run the following command:
```bash
pip install -r requirements.txt
```

git clone https://github.com/jingssrs/pfs_co_fa

## Usage
To perform fiber assignment, use the following command:
```bash
python fiber_assignment.py --input targets.csv --output assignment.csv
```
Replace `targets.csv` with your input file containing the list of targets and `assignment.csv` with the desired output file name.

## Contributing
We welcome contributions to improve the fiber assignment algorithm. Please fork the repository and submit a pull request with your changes.

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
