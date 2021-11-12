#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 21:59:04 2021

@author: rodrigo
"""

import matplotlib as mpl
import matplotlib.colors as co
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

mpl.style.use('default')


def read_JAEA(isotope, yield_type, ZNA, file=None):
    """
    Reads data that is obtained from JAEA at:
    https://wwwndc.jaea.go.jp/cgi-bin/FPYfig
    """

    if file is None:
        file = f'{isotope}_{yield_type}_{ZNA}'

    if ZNA == 'ZN':
        noa = 'N'
    elif ZNA == 'ZA':
        noa = 'A'

    df = pd.read_table(file,
                       sep='\s+',
                       skiprows=5,
                       header=None,
                       names=['Z', noa, 'Fission Yield'])

    Z = np.arange(df['Z'].min(), df['Z'].max()+1)
    NA = np.arange(df[noa].min(), df[noa].max()+1)

    fission_yields = df['Fission Yield'].to_numpy()
    fission_yields = np.reshape(fission_yields, (len(Z), len(NA)))

    organized_df = pd.DataFrame(fission_yields, index=Z, columns=NA)

    return organized_df


def log_tick_formatter(val, pos=None):
    """
    This function helps to avoid the bug in setting the 3D plot to logarithmic
    scale, which causes the labels to change but not the ticks and vales,
    so the result is bad.
    """
    return "{:.0e}".format(10**val)


def plot3D(isotope, yield_type, ZNA, config,
           df=None):

    if df is None:
        df = data[isotope][yield_type][ZNA]

    zmin, zmax = config['zmin'], config['zmax']
    color = config['color']

    fig = plt.figure(figsize=config['figsize'])
    ax = Axes3D(fig)

    ax.set_facecolor('w')

    Z = df.index
    NA = df.columns
    NA_mesh, Z_mesh = np.meshgrid(NA, Z)

    filtered_df = df.applymap(lambda value: np.nan if value < zmin else value)

    # BUG WARNING! https://github.com/matplotlib/matplotlib/issues/209
    # We want to plot with Z axis in logarithmic scale, but the usual
    # functionality of setting the axis to log scale is not working.
    # Our option is to apply the log10 to the Z values, renormalize the colors
    # in log scale but later use the custom tick formater function to show
    # labels as actual values (between 1E-10 and 0.1) instead of
    # log10 of values (-10 and -1).
    # Also the filtered_df has to be converted to numpy 2D array, or else the
    # index of the dataframe becomes weird when we apply np.log10 and
    # the plot functions behave incorrectly.
    log_fission_yields = np.log10(filtered_df.to_numpy())
    log_vmin, log_vmax = np.log10(zmin), np.log10(zmax)
    norm = co.Normalize(vmin=log_vmin, vmax=log_vmax)

    zticks = np.logspace(log_vmin, log_vmax, num=10)
    ax.zaxis.set_major_formatter(mt.FuncFormatter(log_tick_formatter))

    # Start plotting
    surf = ax.plot_surface(NA_mesh, Z_mesh, log_fission_yields,
                           cmap=color,
                           norm=norm,
                           linewidth=1,
                           antialiased=False)

    if ZNA == 'ZN':
        ax.set_xlabel('Neutron Number')
    elif ZNA == 'ZA':
        ax.set_xlabel('Nucleons')
    ax.set_ylabel('Atomic Number')
    ax.set_zlabel('Fission Yield', labelpad=10)

    ax.contourf(NA_mesh, Z_mesh, log_fission_yields,
                zdir='z',
                offset=log_vmin,
                cmap=color,
                norm=norm)

    cbar = plt.colorbar(surf, shrink=0.6, aspect=10)
    cbar.ax.set_yticklabels(["{:.0e}".format(val) for val in zticks])

    E = isotope.rstrip('0123456789')
    A = isotope[len(E):]
    AL = '{' + A + '}'

    plt.title(f'$^{AL}${E} {yield_type} Fission Yield at 0.0253 eV')

    plt.savefig(f'{isotope}_{yield_type}_{ZNA}_3D.pdf')


def heatmap(isotope, yield_type, ZNA, config,
            df=None,
            annotated=False,
            mark_DNP=False):

    if df is None:
        df = data[isotope][yield_type][ZNA]

    vmin, vmax = plot_config['zmin'], plot_config['zmax']
    color = plot_config['color']

    log_norm = co.LogNorm(vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(figsize=plot_config['figsize'])

    def get_edges(intervals):
        edges = intervals.to_numpy() - 0.5
        edges = np.append(edges, edges[-1] + 1)
        return edges

    NA_edges, Z_edges = get_edges(df.columns), get_edges(df.index)

    heatmap = ax.pcolor(NA_edges, Z_edges, df[df > vmin],
                        cmap=color,
                        norm=log_norm,
                        edgecolors='black')

    ax.xaxis.set_minor_locator(mt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mt.MultipleLocator(1))
    ax.set_axisbelow(True)
    ax.grid(which='both')

    if ZNA == 'ZN':
        ax.set_xlabel('Neutron Number')
    elif ZNA == 'ZA':
        ax.set_xlabel('Nucleons')
    ax.set_ylabel('Atomic Number')

    plt.colorbar(heatmap,
                 ticks=np.logspace(np.log10(vmin), np.log10(vmax), num=10))

    E = isotope.rstrip('0123456789')
    A = isotope[len(E):]
    AL = '{' + A + '}'

    plt.title(f'$^{AL}${E} {yield_type} Fission Yield at 0.0253 eV')

    save_name = f'{isotope}_{yield_type}_{ZNA}_HM'

    if annotated:
        save_name += '_annotated'
        value_to_color = mpl.cm.ScalarMappable(norm=log_norm, cmap='Greys')

        for na in df.columns:
            for z in df.index:
                val = df.at[z, na]
                if val > vmin:
                    plt.text(na, z, '%.2e' % val,
                             color=value_to_color.to_rgba(val),
                             fontsize=1.2,
                             horizontalalignment='center',
                             verticalalignment='center',
                             )

    if mark_DNP:
        save_name += '_DNP'
        for key, values in DNP[ZNA].items():
            diff_values = np.diff(values)
            split_index = np.where(diff_values > 1)[0]
            split_values = np.split(values, split_index+1)

            for interval in split_values:
                na_n = np.array(interval) - 0.5
                na_p = np.array(interval) + 0.5
                na_edges = np.unique(np.concatenate([na_n, na_p]))
                ax.fill_between(na_edges, key - 0.5, key + 0.5,
                                fc='none', ec='k', hatch='////')

    plt.savefig(save_name + '.pdf',
                bbox_inches='tight',
                transparent=True,)


# HNRC Vol2 Pg259
# DNP responsible for 95%+ of DN
DNP = {'ZA': {31: [79, 80, 81, 82],
              33: [84, 85, 86, 87],
              34: [87, 88, 89, 91],
              35: [87, 88, 89, 90, 91, 92],
              36: [92, 93, 94],
              37: [92, 93, 94, 95, 96, 97, 98],
              38: [97, 98],
              39: [98, 99],
              49: [127, 128, 129, 130, 131, 132],
              50: [134],
              51: [134, 135, 136],
              52: [136, 137, 138],
              53: [137, 138, 139, 140, 141],
              54: [141, 142],
              55: [141, 142, 143, 144, 145, 146],
              56: [147, 148],
              57: [147]},
       'ZN': {}
       }

for key, values in DNP['ZA'].items():
    N = np.array(values) - key
    DNP['ZN'][key] = N

# Standard data automatically available
isotopes = ['U233', 'U235', 'Pu239']
yield_types = ['IND', 'CUM']
ZNA = ['ZN', 'ZA']

# Nested dictionary that is accessed as
# data[isotopes][yield_types][ZNA]
# each element is a dataframe with Z index, and N or A columns
data = {I: {yt: {zna: read_JAEA(I, yt, zna)
                 for zna in ZNA} for yt in yield_types} for I in isotopes}

def plot_everything():
    for isotope in isotopes:
        for yield_type in yield_types:
            for zna in ZNA:
                plot3D(isotope, yield_type, zna, plot_config)
                heatmap(isotope, yield_type, zna, plot_config)


# Default plot configurations
plot_config = {'figsize': (10, 6),
               'zmin': 1E-10,
               'zmax': 1E-1,
               'color': mpl.cm.inferno}





isotope = 'U235'
yield_type = 'IND'
zna = 'ZN'

plot3D(isotope, yield_type, zna, plot_config)
heatmap(isotope, yield_type, zna, plot_config)

heatmap(isotope, yield_type, zna, plot_config, mark_DNP=True)
heatmap(isotope, yield_type, zna, plot_config, annotated=True)

isotope = 'Cf242'
yield_type = 'IND'
zna = 'ZN'

Cf242 = read_JAEA('Cf242', 'IND', 'ZN', file='Cf242')
plot3D(isotope, yield_type, ZNA, plot_config, df=Cf242)

plot_everything()
