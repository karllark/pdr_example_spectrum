#!/usr/bin/env python

from __future__ import (absolute_import, print_function, division)

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from astropy.table import Table

def set_params(lw=1.5, universal_color='#262626', fontsize=16):
    '''Configure some matplotlib rcParams.

    Parameters
    ----------
    lw : scalar
        Linewidth of plot and axis lines. Default is 1.5.
    universal_color : str, matplotlib named color, rgb tuple
        Color of text and axis spines. Default is #262626, off-black
    fontsize : scalar
        Font size in points. Default is 12
    '''
    rc('font', size=fontsize)
    rc('lines', linewidth=lw, markeredgewidth=lw*0.5)
    rc('patch', linewidth=lw, edgecolor='#FAFAFA')
    rc('axes', linewidth=lw, edgecolor=universal_color,
       labelcolor=universal_color, 
       axisbelow=True)
    rc('image', origin='lower') # fits images
    rc('xtick.major', width=lw*0.75)
    rc('xtick.minor', width=lw*0.5)
    rc('xtick', color=universal_color)
    rc('ytick.major', width=lw*0.75)
    rc('ytick.minor', width=lw*0.5)
    rc('ytick', color=universal_color)
    rc('grid', linewidth=lw)
    rc('legend', loc='best', numpoints=1, scatterpoints=1, handlelength=1.5,
        fontsize=fontsize, columnspacing=1, handletextpad=0.75)

def initialize_parser():
    '''For running from command line, initialize argparse with common args
    '''
    ftypes = ['png', 'jpg', 'jpeg', 'pdf', 'ps', 'eps', 'rgba',
              'svg', 'tiff', 'tif', 'pgf', 'svgz', 'raw']
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--savefig', action='store', 
                        default=False, choices=ftypes,
                        help='Save figure to a file')
    return parser

if __name__ == '__main__':

    parser = initialize_parser()
    args = parser.parse_args()

    pdr_data_nir = Table.read('data/Table_PDR_gas_lines_R3000_NIR.dat',
                              format='ascii.commented_header',
                              header_start=-1)
    pdr_data_mir = Table.read('data/Table_PDR_gas_lines_R3000_MIR.dat',
                              format='ascii.commented_header',
                              header_start=-1)

    cloudy_data = Table.read('data/cloudy_continuumcut.dat',
                             format='ascii.commented_header',
                             header_start=-1)

    pah_data = Table.read('data/PAH_template_JWST.txt',
                          format='ascii.commented_header',
                          header_start=-1)

    dust_data = Table.read('data/DIRTY_bm1_neqgrain_tau_001.00_theta_090_global_lum.table.fits')
    print(dust_data.colnames)

    xsize = 11.0
    ysize = 11.0
    fig, ax = plt.subplots(nrows=5,figsize=(xsize,ysize), sharex=True)

    set_params()

    ptype = 'linear'

    ax[0].plot(cloudy_data['wave'], cloudy_data['diffuse_spec'], 'k-')
    ax[0].set_xscale('log')
    ax[0].set_xlim(1.0,29.0)
    ax[0].set_yscale(ptype)
    ax[0].set_ylim(1e-3,50)

    ax[1].plot(pdr_data_nir['wave'], pdr_data_nir['atomic_pdr'], 'k-')
    ax[1].plot(pdr_data_mir['wave'], pdr_data_mir['atomic_pdr'], 'k-')
    ax[1].set_yscale(ptype)
    ax[1].set_ylim(1e-4,0.04)

    ax[2].plot(pdr_data_nir['wave'], pdr_data_nir['mol_pdr'], 'k-')
    ax[2].plot(pdr_data_mir['wave'], pdr_data_mir['mol_pdr'], 'k-')
    ax[2].set_yscale(ptype)
    ax[2].set_ylim(1e-4,0.08)

    ax[3].plot(pah_data['wave'], pah_data['pah_spec'], 'k-')
    ax[3].set_yscale(ptype)
    ax[3].set_ylim(1e-2,1e3)

    ax[4].plot(dust_data['wavelength'],
               dust_data['wavelength']*(dust_data['Flux_rt_s']
                                        + dust_data['Flux_de_d_1']), 'k-')
    ax[4].plot(dust_data['wavelength'],
               dust_data['wavelength']*dust_data['Flux_rt_s'], 'k--')
    ax[4].plot(dust_data['wavelength'],
               dust_data['wavelength']*dust_data['Flux_de_d_1'], 'k:')
    ax[4].plot(dust_data['wavelength'],
               dust_data['wavelength']*dust_data['Flux_de_d_3'], 'k-.')
    ax[4].plot(dust_data['wavelength'],
               dust_data['wavelength']*dust_data['Flux_de_d_5'], 'k:')
    ax[4].set_yscale('log')
    ax[4].set_ylim(1e35,2e37)
    ax[4].set_xlabel('wavelength [$\mu$m]')

    fig.tight_layout(h_pad=0.15)

 
    # save the plot
    basename = 'pdr_example'
    if args.savefig:
        fig.savefig('{}.{}'.format(basename, args.savefig))
    else:
        plt.show()
