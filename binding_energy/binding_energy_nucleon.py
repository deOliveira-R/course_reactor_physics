# -*- coding: utf-8 -*-
"""
Calculate the Binding Energy per Nucleon (BEN)
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

mpl.style.use('seaborn')


data = pd.read_csv("binding_energy.csv")

ax = data.plot(x='Nucleons',
               y='BEN',
               legend=None,
               marker='d',
               markersize=3,
               markeredgecolor='k',
               xticks=[1, 25, 50, 75, 100, 125, 150, 175, 200, 225],
               yticks=[0, 1, 2, 3, 4, 5, 6, 7, 8])

plt.vlines(117, 0, 8.519, ls='--', colors='green')

plt.xlim(xmin=0, xmax=250)
plt.ylim(ymin=0)
plt.ylabel('Binding Energy per Nucleon (MeV)')

data = data.set_index('Isotope')

print(data.loc['He-4', 'Nucleons'])

ax.annotate('$^2$H',
            (data.loc['H-2', 'Nucleons'], data.loc['H-2', 'BEN']),
            xycoords="data",
            xytext=(4, 1.1),
            textcoords='data',
#            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
            size=8,
            ha='left',
            va='center')

ax.annotate('$^3$H',
            (data.loc['H-3', 'Nucleons'], data.loc['H-3', 'BEN']),
            xycoords="data",
            xytext=(5, 2.8),
            textcoords='data',
#            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
            size=8,
            ha='left',
            va='center')

ax.annotate('$^3$He',
            (data.loc['He-3', 'Nucleons'], data.loc['He-3', 'BEN']),
            xycoords="data",
            xytext=(5, 2.5),
            textcoords='data',
#            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
            size=8,
            ha='left',
            va='center')

ax.annotate('$^4$He',
            (data.loc['He-4', 'Nucleons'], data.loc['He-4', 'BEN']),
            xycoords="data",
            xytext=(4, 7.3),
            textcoords='data',
#            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
            size=8,
            ha='center',
            va='center')

ax.annotate('$^6$Li',
            (data.loc['Li-6', 'Nucleons'], data.loc['Li-6', 'BEN']),
            xycoords="data",
            xytext=(8, 5.2),
            textcoords='data',
#            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
            size=8,
            ha='left',
            va='center')

ax.annotate('$^7$Li',
            (data.loc['Li-7', 'Nucleons'], data.loc['Li-7', 'BEN']),
            xycoords="data",
            xytext=(9, 5.6),
            textcoords='data',
#            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
            size=8,
            ha='left',
            va='center')

ax.annotate('$^{12}$C',
            (data.loc['C-12', 'Nucleons'], data.loc['C-12', 'BEN']),
            xycoords="data",
            xytext=(6, 8.3),
            textcoords='data',
            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
            bbox=dict(pad=-2, facecolor="none", edgecolor="none"),
            size=8,
            ha='center',
            va='center')

ax.annotate('$^{16}$O',
            (data.loc['O-16', 'Nucleons'], data.loc['O-16', 'BEN']),
            xycoords="data",
            xytext=(14, 8.7),
            textcoords='data',
            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
            bbox=dict(pad=-2, facecolor="none", edgecolor="none"),
            size=8,
            ha='center',
            va='center')

ax.annotate('$^{56}$Fe',
            (data.loc['Fe-56', 'Nucleons'], data.loc['Fe-56', 'BEN']),
            xytext=(56, 9.4),
            textcoords='data',
            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
            bbox=dict(pad=-2, facecolor="none", edgecolor="none"),
            size=8,
            ha='center',
            va='center')


ax.annotate('$^{235}$U\n' + f'({data.loc["U-235","BEN"]} MeV)',
            (data.loc['U-235', 'Nucleons'], data.loc['U-235', 'BEN']),
            xytext=(230, 6.5),
            textcoords='data',
            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
#            bbox=dict(pad=-1, facecolor="none", edgecolor="none"),
            size=8,
            ha='center',
            va='center')

ax.annotate('Average\nFission Product\n(8.519 MeV)',
            (117, 8.519),
            xytext=(140, 7),
            textcoords='data',
            arrowprops=dict(arrowstyle='-|>', edgecolor='k', facecolor='black'),
#            bbox=dict(pad=-1, facecolor="none", edgecolor="none"),
            size=8,
            ha='center',
            va='center')

plt.savefig('binding_energy.pdf', bbox='tight')
