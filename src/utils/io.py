from .common_imports import np, time, plt, fits, Table, new_colors
import pandas as pd
from astropy.table import Column
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from matplotlib.path import Path
import astropy.units as u

# TODO: same as in pointings.ipynb
def plot_filter_fov(raP, decP, raSci, decSci, PA=0, n_sci_fov_least=3000, filter_fov=True):
    '''plot one PFS FoV (hexagon) centered at the pointing center
    
    NOTE
    ==========
    flag_fov_reserved is obtained by using a threshold of targets in the FoV

    Parameters
    ==========
    raP, decP, PA : float
        ra, dec, PA of the pointing center

    raSci, decSci: numpy array, float
        ra, dec of the scientific targets
        only used to check the number of scientific targets in the FoV
    
    n_sci_fov_least: int
        the least number of scientific targets in the FoV

    filter_fov: Boolean
        if True, plot/select the FoV only when there are enough scientific targets in the FoV

        
    Returns
    =======
    plot a hexagon at the pointing center with diameter=1.38 deg
    
    flag_fov_reserved: Boolean, used to remove the pointing w/o enough scientific targets 
    '''
    
    center = SkyCoord(raP*u.deg, decP*u.deg)
    # PA=0 along y-axis, PA=90 along x-axis, PA=180 along -y-axis...
    hexagon = center.directional_offset_by([0+PA, 60+PA, 120+PA, 180+PA, 240+PA, 300+PA, 360+PA]*u.deg, 1.38/2.*u.deg)
    ra_h = hexagon.ra.deg
    dec_h = hexagon.dec.deg

    ra_h_in = np.where(np.fabs(ra_h-center.ra.deg)>180)
    if len(ra_h_in[0])>0:
        if ra_h[ra_h_in[0][0]]>180:ra_h[ra_h_in[0]]-=360
        elif ra_h[ra_h_in[0][0]]<180:ra_h[ra_h_in[0]]+=360
        #pdb.set_trace()

    # note ra range is [-180, 180] for HSC wide autumn field
    ra_h[ra_h>300] -= 360

    # scientific targets
    point = np.vstack((raSci, decSci)).T
    

    if filter_fov:
        polygon = Path([(ra_h[t],dec_h[t]) for t in range(len(ra_h))])
        index_ = np.where(polygon.contains_points(point)==True)[0]

        if(len(index_)<n_sci_fov_least):
            flag_fov_reserved = False
        else:
            flag_fov_reserved = True
            plt.plot(ra_h, dec_h, color='r', lw=0.5, ls='-', alpha=1., zorder=5)
    else:
        flag_fov_reserved = True
        plt.plot(ra_h, dec_h, color='r', lw=0.5, ls='-', alpha=1., zorder=5)
    
    return flag_fov_reserved

# TODO: modify the read_pointings function (return xxxx)
# read the pointing centers from the file
def read_pointings(file):
    """
    Read pre-defined pointings from a file
    """
    pointings = ascii.read(file)

    return pointings

def select_pointings(pointings, ra_range, dec_range):
    '''
    select pointings within the range
    '''
    ra_peaks = pointings['ppc_ra']
    dec_peaks = pointings['ppc_dec']
    # wrap the ra to [-180, 180] for the HSC wide autumn field
    ra_peaks[ra_peaks>300] -= 360

    mask_peaks = (ra_peaks>ra_range[0]) & (ra_peaks<ra_range[1]) & (dec_peaks>dec_range[0]) & (dec_peaks<dec_range[1])
    print("Select %d pointings." % np.sum(mask_peaks))
    
    return pointings[mask_peaks]


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


def plot_tgt_done(outfn_list, figname, figsize=(8, 8), plot_diffcolor=True):
    '''
    plot the targets that have been assigned to the fibers
    '''
    for i, fn_i in enumerate(outfn_list):
        tgt_id_done, tgt_ra_done, tgt_dec_done = np.loadtxt(fn_i, usecols=(0, 4, 5), unpack=True, dtype='str')
        tgt_ra_done, tgt_dec_done = tgt_ra_done.astype('float'), tgt_dec_done.astype('float')

        mask_cos = [tgt_id_done[i][1:4] == 'Cos' for i in range(len(tgt_id_done))]
        mask_star = [tgt_id_done[i][1:4] == 'Sta' for i in range(len(tgt_id_done))]
        mask_sky = [tgt_id_done[i][1:4] == 'Sky' for i in range(len(tgt_id_done))]
        mask_anc = [tgt_id_done[i][1:4] == 'Anc' for i in range(len(tgt_id_done))]

        if(plot_diffcolor):
            if(i==0): 
                plt.figure(figsize = figsize)
                plt.plot(tgt_ra_done[mask_cos], tgt_dec_done[mask_cos], 'k.', ms=0.5, alpha=0.5, label='cosmology')
                plt.plot(tgt_ra_done[mask_star], tgt_dec_done[mask_star], 'r*', ms=1.5, alpha=1., label='star')
                plt.plot(tgt_ra_done[mask_sky], tgt_dec_done[mask_sky], 'b^', ms=1.5, alpha=1., label='sky')
                plt.plot(tgt_ra_done[mask_anc], tgt_dec_done[mask_anc], 'gs', ms=1.5, alpha=1., label='ancillary')
            else:
                plt.plot(tgt_ra_done[mask_cos], tgt_dec_done[mask_cos], 'k.', ms=0.5, alpha=0.5)
                plt.plot(tgt_ra_done[mask_star], tgt_dec_done[mask_star], 'r*', ms=1.5, alpha=1.)
                plt.plot(tgt_ra_done[mask_sky], tgt_dec_done[mask_sky], 'b^', ms=1.5, alpha=1.)
                plt.plot(tgt_ra_done[mask_anc], tgt_dec_done[mask_anc], 'gs', ms=1.5, alpha=1.)
        else:
            if(i==0): 
                plt.figure(figsize = figsize)
                plt.plot(tgt_ra_done, tgt_dec_done, 'k.', ms=0.5, alpha=0.5, label='targets done')
            else:
                plt.plot(tgt_ra_done, tgt_dec_done, 'k.', ms=0.5, alpha=0.5)

    plt.legend(loc='upper right', fontsize=15, frameon=True)
    plt.xlabel('RA', fontsize=15)
    plt.ylabel('DEC', fontsize=15)
    plt.savefig(figname, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()