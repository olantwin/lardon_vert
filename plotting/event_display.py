import config as cf
import data_containers as dc

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import collections  as mc
import itertools as itr
import math
import colorcet as cc

import matplotlib.patches as patches


def_adc_min = -150 #-10
def_adc_max = 150 #10
cmap_ed = cc.cm.diverging_tritanopic_cwr_75_98_c20



def draw(view, ax=None, t_min = -1, t_max = -1, ch_min = -1, ch_max = -1, adc_min=def_adc_min, adc_max=def_adc_max):

    ax = plt.gca() if ax is None else ax
    
    ax.imshow(dc.data[view, :, :].transpose(), 
              origin = 'lower', 
              aspect = 'auto', 
              interpolation='none',
              cmap   = cmap_ed,
              vmin   = adc_min, 
              vmax   = adc_max)


    if(ch_min > -1):
        ax.set_xbound(lower=ch_min)
    if(ch_max > -1):
        ax.set_xbound(upper=ch_max)

    if(t_min > -1):
        ax.set_ybound(lower=t_min)
    if(t_max > -1):
        ax.set_ybound(upper=t_max)
    
    return ax


def plot_ed_zoom(view, t_min = -1, t_max = -1, ch_min = -1, ch_max = -1, option=None, to_be_shown=False):

    fig = plt.figure( figsize=(6,6))
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[1,30])
    
    ax_data = fig.add_subplot(gs[1, 0])
    ax_data = draw(view, ax=ax_data, t_min = t_min, t_max = t_max, ch_min = ch_min, ch_max = ch_max)
    ax_data.set_xlabel('View Channel')
    ax_data.set_ylabel('Time')
    ax_data.yaxis.set_major_locator(plt.MaxNLocator(4))
    ax_data.set_title('View '+str(view) + '['+cf.view_type[view]+']')

    ax_col = fig.add_subplot(gs[0, 0])
    ax_col.set_title('Collected Charge [ADC]')
        
                         
    cb = fig.colorbar(ax_data.images[-1], cax=ax_col, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')

    plt.tight_layout()

    if(option):
        option = "_"+option
    else:
        option = ""

    run_day = str(dc.evt_list[-1].date)
    evt_nb = str(dc.evt_list[-1].evt_nb)



    plt.savefig('ED/ed_zoom_v'+str(view)+option+'_'+run_day+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



def plot_ed(nameout="", option=None, to_be_shown=False, adc_min=def_adc_min, adc_max=def_adc_max):


    fig = plt.figure(figsize=(9,4))
    gs = gridspec.GridSpec(nrows=2, 
                           ncols=2, 
                           height_ratios=[1, 10])
    ax_col = fig.add_subplot(gs[0,:])
    ax_v0  = fig.add_subplot(gs[1, 0])
    ax_v1  = fig.add_subplot(gs[1, 1], sharey=ax_v0)
    

    """ draw data """

    """ View 0 """
    ax_v0 = draw(0, ax=ax_v0, adc_min=adc_min, adc_max=adc_max)
    ylow, yup = ax_v0.get_ylim()
    ax_v0.plot([32, 32], [ylow, yup], ls='--', c='k')
    ax_v0.set_title('View 0 ' + '['+cf.view_type[0]+']')
    ax_v0.set_ylabel('Time')

    """ View 1 """
    ax_v1 = draw(1, ax=ax_v1, adc_min=adc_min, adc_max=adc_max)
    ax_v1.set_title('View 1 ' + '['+cf.view_type[1]+']')
    ax_v1.set_ylabel('Time')
    ax_v1.yaxis.tick_right()
    ax_v1.yaxis.set_label_position("right")

        
    ax_v0.set_xlabel('View Channel')
    ax_v1.set_xlabel('View Channel')



    ax_v0.yaxis.set_major_locator(plt.MaxNLocator(4))

    ax_col.set_title('Collected Charge [ADC]')
        
                         
    cb = fig.colorbar(ax_v0.images[-1], cax=ax_col, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
        
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.02)

    """
    plt.subplots_adjust(wspace=0.02, 
                        hspace=0.25, 
                        top=0.92, 
                        bottom=0.08, 
                        left=0.1, 
                        right=0.9)
    """

    if(option):
        option = "_"+option
    else:
        option = ""

    run_day = str(dc.evt_list[-1].date)
    evt_nb = str(dc.evt_list[-1].evt_nb)

    plt.savefig('ED/ed_'+nameout+option+'_'+run_day+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()


def plot_ed_data(option=None, to_be_shown=False, adc_min=def_adc_min, adc_max=def_adc_max):
    plot_ed(nameout="data", option=option, to_be_shown=to_be_shown, adc_min=adc_min, adc_max=adc_max)




def plot_ed_data_and_hits(option=None, to_be_shown=False, adc_min=def_adc_min, adc_max=def_adc_max):


    fig = plt.figure(figsize=(9,4))
    gs = gridspec.GridSpec(nrows=2, 
                           ncols=2, 
                           height_ratios=[1, 10])
    ax_col = fig.add_subplot(gs[0,:])
    ax_v0  = fig.add_subplot(gs[1, 0])
    ax_v1  = fig.add_subplot(gs[1, 1], sharey=ax_v0)
    

    r0, r1 = [], []
    for h in dc.hits_list:
        if(h.view==0):
            r0.append(patches.Rectangle((h.channel-0.5,h.start),1,h.stop-h.start,linewidth=1,edgecolor='k',facecolor='none'))
        else:
            r1.append(patches.Rectangle((h.channel-0.5,h.start),1,h.stop-h.start,linewidth=1,edgecolor='k',facecolor='none'))



    """ draw data """

    """ View 0 """
    ax_v0 = draw(0, ax=ax_v0, adc_min=adc_min, adc_max=adc_max)
    ylow, yup = ax_v0.get_ylim()
    ax_v0.plot([32, 32], [ylow, yup], ls='--', c='k')
    ax_v0.set_title('View 0 ' + '['+cf.view_type[0]+']')
    ax_v0.set_ylabel('Time')
    for r in r0:
        ax_v0.add_patch(r)

    """ View 1 """
    ax_v1 = draw(1, ax=ax_v1, adc_min=adc_min, adc_max=adc_max)
    ax_v1.set_title('View 1 ' + '['+cf.view_type[1]+']')
    ax_v1.set_ylabel('Time')
    ax_v1.yaxis.tick_right()
    ax_v1.yaxis.set_label_position("right")
    for r in r1:
        ax_v1.add_patch(r)

        
    ax_v0.set_xlabel('View Channel')
    ax_v1.set_xlabel('View Channel')



    ax_v0.yaxis.set_major_locator(plt.MaxNLocator(4))

    ax_col.set_title('Collected Charge [ADC]')
        
                         
    cb = fig.colorbar(ax_v0.images[-1], cax=ax_col, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
        
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.02)

    """
    plt.subplots_adjust(wspace=0.02, 
                        hspace=0.25, 
                        top=0.92, 
                        bottom=0.08, 
                        left=0.1, 
                        right=0.9)
    """

    if(option):
        option = "_"+option
    else:
        option = ""

    run_day = str(dc.evt_list[-1].date)
    evt_nb = str(dc.evt_list[-1].evt_nb)

    plt.savefig('ED/ed_data_and_hits'+option+'_'+run_day+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()

