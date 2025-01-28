'''
Common imports for the project
'''
# basic 
import numpy as np
import time
import os

# io 
from astropy.io import fits
from astropy.table import Table, column
from astropy.io import fits, ascii

# plot 
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc

# Set default font properties
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica'], 'size': 20})
rc('text', **{'usetex':'False'})
rc('axes', **{'grid':'True','labelsize':24,'linewidth':1.5,})
rc('xtick',**{'direction':'in','labelsize':20,'minor.visible':True,'major.size': 10,'major.width': 1.5, \
              'minor.size': 3.72,'major.width': 0.8})
rc('ytick',**{'direction':'in','labelsize':20,'minor.visible':True,'major.size': 10,'major.width': 1.5, \
              'minor.size': 3.72,'major.width': 0.8})
rc('legend',**{'frameon':False,'numpoints':1,'scatterpoints':1,'fontsize':18})
rc('figure',**{'autolayout':True})
rc('grid',**{'alpha':0.6,'linewidth':0.8,'linestyle':'--'})

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
new_colors = [mpl.colors.to_rgba(color_i) for color_i in new_colors]