from utils.common_imports import np, time, plt, fits, Table, new_colors
import pandas as pd
from astropy.table import Column
from astropy.io import ascii


def plot_radec(dt, title, save_fig=True, output_dir='../output/figures/'):
    '''
    Plot the ra-dec distribution of the data.
    '''
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.set_xlabel('RA [deg]')
    ax.set_ylabel('Dec [deg]')
    ax.set_title(title)
    
    ax.plot(dt['R.A.'], dt['Dec.'], ",", alpha=1., ls='none')

    if save_fig:
        fig.savefig(output_dir + title + '_radec.png')
    plt.show()


def write_data_table(data_table, output_dir='../data_proc/', prefix=None, fmt='fits'):
    '''
    Write the calibration/cosmology data table into file that is ready for netflow.
    '''
    if prefix is None:
        raise ValueError('Give prefix (cos, sky, or star) of data file!')

    # Write the data table
    output_fn = output_dir + prefix + '_targets.' + fmt
    if fmt == 'fits':
        data_table.write(output_fn, format='fits', overwrite=True)
    elif fmt == 'ecsv':
        data_table.write(output_fn, format='ascii.ecsv', overwrite=True)
    else:
        raise ValueError('The format is not supported!')
    
    print('Data table is written to ' + output_fn)
    


    