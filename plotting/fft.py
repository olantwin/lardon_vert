import config as cf
import data_containers as dc

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import collections  as mc
import itertools as itr
import math
import colorcet as cc


adc_min = -10
adc_max = 10
cmap_fft = cc.cm.fire_r

def plot_fft_ps(ps, zmax=1., option=None, to_be_shown=False):

    fig = plt.figure(figsize=(9,4))
    gs = gridspec.GridSpec(nrows=2, 
                           ncols=2, 
                           height_ratios=[1, 10])

    ax_col = fig.add_subplot(gs[0,:])
    ax_v0  = fig.add_subplot(gs[1, 0])
    ax_v1  = fig.add_subplot(gs[1, 1], sharey=ax_v0)


    ax_v0.imshow(ps[0, :, :].transpose(), 
                 origin = 'lower', 
                 aspect = 'auto', 
                 cmap   = cmap_fft,
                 vmin   = 0., 
                 extent = [0, cf.n_ChanPerView, 0., zmax])

    ax_v0.set_title('View 0 ' + '['+cf.view_type[0]+']')
    ax_v0.set_ylabel('Signal Frequencies [MHz]')



    ax_v1.imshow(ps[1, :, :].transpose(), 
                 origin = 'lower', 
                 aspect = 'auto', 
                 cmap   = cmap_fft,
                 vmin   = 0., 
                 extent = [0, cf.n_ChanPerView, 0., zmax])

    ax_v1.set_title('View 1 ' + '['+cf.view_type[1]+']')
    ax_v1.set_ylabel('Signal Frequencies [MHz]')
    ax_v1.yaxis.tick_right()
    ax_v1.yaxis.set_label_position("right")


    ax_v0.set_xlabel('View Channel')
    ax_v1.set_xlabel('View Channel')



    ax_v0.yaxis.set_major_locator(plt.MaxNLocator(4))

    ax_col.set_title('Intensity')
        
                         
    cb = fig.colorbar(ax_v0.images[-1], cax=ax_col, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
        
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.02)


    if(option):
        option = "_"+option
    else:
        option = ""

    run_day = str(dc.evt_list[-1].date)
    evt_nb = str(dc.evt_list[-1].evt_nb)

    plt.savefig('ED/fft'+option+'_'+run_day+'_evt_'+evt_nb+'.png')

    if(to_be_shown):
        plt.show()
    plt.close()
