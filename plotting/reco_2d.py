import config as cf
import data_containers as dc

from .select_hits import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import collections  as mc
from matplotlib.legend_handler import HandlerTuple
import itertools as itr
import math
import colorcet as cc


adc_min = -5
adc_max = 20
cmap_ed = cc.cm.kbc_r
marker_size = 5
"""define default color cycle from colorcet"""
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=cc.glasbey_warm)

color_noise     = "#c7c7c7"
color_clustered = "#ffb77d"
color_matched1  = "#28568f"
color_matched2  = "#abdf7f"
color_track2d   = "#de425b"
color_track3d   = "#00a9b2"

def draw_hits(pos, time, z=[], ax=None, **kwargs):

    ax = plt.gca() if ax is None else ax

    if(len(pos) == 0):
        return ax
    
    if(len(z) > 0):
        ax.scatter(pos, time, c=z, **kwargs)
    else:
        ax.scatter(pos, time, **kwargs)
    return ax

def draw_all_hits(ax_v0, ax_v1, sel='True', adc=False, charge=False, **kwargs):

    axs = [ax_v0, ax_v1]

    for iview in range(2):
        z = []
        if(adc==True):
            z = get_hits_adc(iview, sel)
        elif(charge==True):
            z = get_hits_charge(iview, sel)

        axs[iview] = draw_hits(pos=get_hits_pos(iview, sel), 
                               time=get_hits_z(iview, sel), 
                               z=z,
                               ax=axs[iview], **kwargs)
    return axs[0], axs[1]



def draw_tracks(pos, time, ax=None, legend="", **kwargs):

    ax = plt.gca() if ax is None else ax
    
    if(len(pos) == 0):
        return ax

    if(len(legend)>0):
        ax.plot(pos[0], time[0], label=legend, **kwargs)
        
    for tx,tz in zip(pos, time):
        ax.plot(tx,tz, **kwargs)
    
    return ax


def draw_all_tracks(ax_v0, ax_v1, sel='True', legend="", **kwargs):

    axs = [ax_v0, ax_v1]

    for iview in range(2):
        axs[iview] = draw_tracks(pos=get_2dtracks_pos(iview,sel), 
                                 time=get_2dtracks_z(iview,sel), 
                                 ax=axs[iview], 
                                 legend=legend,
                                 **kwargs)
    return axs[0], axs[1]



def template_data_view():

    fig = plt.figure(figsize=(9,4))
    gs = gridspec.GridSpec(nrows=2, ncols=2, 
                           height_ratios=[1,10])

    
    ax_col = fig.add_subplot(gs[0,:])
    ax_v0 = fig.add_subplot(gs[1, 0])
    ax_v1 = fig.add_subplot(gs[1, 1], sharey=ax_v0)

    ax_v0.set_title('View 0 ['+cf.view_type[0]+']')
    ax_v0.set_ylabel('Z [cm]')
    ax_v0.set_xlabel('X [cm]')
    ax_v0.set_xlim([0, cf.len_det_x])
    ax_v0.set_ylim([0., cf.Anode_Z])

    ax_v1.set_title('View 1 ['+cf.view_type[1]+']')
    ax_v1.set_ylabel('Z [cm]')
    ax_v1.yaxis.tick_right()
    ax_v1.yaxis.set_label_position("right")
    ax_v1.set_xlabel('Y [cm]')
    ax_v1.set_xlim([0, cf.len_det_y])
    ax_v1.set_ylim([0., cf.Anode_Z])

    plt.tight_layout()
    #plt.subplots_adjust(wspace=0.02)
    
    return fig, ax_col, ax_v0, ax_v1



def plot_2dview_hits(option=None, to_be_shown=False):
    if(dc.evt_list[-1].nHits[0] == 0):
        return

    fig, ax_col, ax_v0, ax_v1 = template_data_view()

    max_adc=50
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, adc=True, cmap=cmap_ed, s=marker_size, vmin=0, vmax=max_adc)
    
    """ color bar """
    ax_col.set_title('Hit Max ADC')

    cb = fig.colorbar(ax_v0.collections[0], cax=ax_col, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')


    if(option):
        option = "_"+option
    else:
        option = ""


    run_day = str(dc.evt_list[-1].date)
    evt_nb = str(dc.evt_list[-1].evt_nb)

    plt.tight_layout()

    plt.savefig('ED/hit_view'+option+'_'+run_day+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



def plot_2dview_clusters(option=None, to_be_shown=False):
    
    fig, ax_col, ax_v0, ax_v1 = template_data_view()

    """ clustered hits """

    for icl in range(dc.evt_list[-1].nClusters[0]):
        sel = 'x.view == 0 and x.cluster=='+str(icl)
        ax_v0 = draw_hits(pos=get_hits_pos(0, sel), time=get_hits_z(0, sel), ax=ax_v0, s=marker_size, marker='o')

    for icl in range(dc.evt_list[-1].nClusters[1]):
        sel = 'x.view == 1 and x.cluster=='+str(icl)        
        ax_v1 = draw_hits(pos=get_hits_pos(1, sel), time=get_hits_z(1, sel), ax=ax_v1, s=marker_size, marker='o')



    """ unclustered hits """
    sel = 'x.cluster==-1'
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_noise, s=marker_size, marker='o')

    ax_col.axis('off')

    if(option):
        option = "_"+option
    else:
        option = ""


    run_day = str(dc.evt_list[-1].date)
    evt_nb = str(dc.evt_list[-1].evt_nb)


    plt.savefig('ED/cluster_2dview'+option+'_run_'+run_day+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()






def plot_2dview_2dtracks(option=None, to_be_shown=False):
    fig, ax_leg, ax_v0, ax_v1 = template_data_view()
    
    """ all hits """
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, c="#e6e6e6", s=marker_size, marker='o', label='Noise Hits')

    
    
    """ 2D tracks """
    ax_v0, ax_v1 = draw_all_tracks(ax_v0, ax_v1, linewidth=1, legend='2D Track')


    """ legend """
    ax_leg.axis('off')

    
    if(option):
        option = "_"+option
    else:
        option = ""


    run_day = str(dc.evt_list[-1].date)
    evt_nb = str(dc.evt_list[-1].evt_nb)

    plt.tight_layout()

    plt.savefig('ED/alltrack2D'+option+'_'+run_day+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



def plot_2dview_hits_2dtracks(option=None, to_be_shown=False):
    fig, ax_leg, ax_v0, ax_v1 = template_data_view()
    
    """ unclustered hits """
    sel = 'x.cluster == -1'    
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_noise, s=marker_size, marker='o', label='Noise Hits')

    """ clustered hits """
    sel = 'x.cluster > -1'    
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_clustered, s=marker_size, marker='o', label='Hits Clustered')

    """ delta_ray hits attached to track """
    sel = 'x.matched <0 and x.matched > -9999'
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_matched2, s=marker_size, marker='o', label='Delta Rays')

    """ hits attached to track """
    sel = 'x.matched >= 0'
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_matched1, s=marker_size, marker='o', label='Hits Attached to Track')

    
    """ 2D tracks """
    ax_v0, ax_v1 = draw_all_tracks(ax_v0, ax_v1, legend='2D Track', c=color_track2d, linewidth=1)


    """ legend """
    ax_leg.axis('off')
    
    """ re-arrange the legend (line last), and merge blue and green entries """
    """ might not work anymore if the plotting order is changed """
    h, l = ax_v0.get_legend_handles_labels()

    if(False): #len(h)==5):
        leg = ax_leg.legend([h[1], h[2], (h[3], h[4]), h[0]], [l[1], l[2], 'Hits Attached to Track (1,2)', l[0]], loc='center', ncol=4, markerscale=4, handler_map={tuple: HandlerTuple(ndivide=None)})
    else:
        """otherwise this works well """
        leg = ax_leg.legend(*ax_v0.get_legend_handles_labels(),loc='center', ncol=5, markerscale=4, markerfirst=True)
    
    """ make line in the legend bigger """
    for line in leg.get_lines():
        line.set_linewidth(3)

    plt.tight_layout()
    
    if(option):
        option = "_"+option
    else:
        option = ""


    run_day = str(dc.evt_list[-1].date)
    evt_nb = str(dc.evt_list[-1].evt_nb)

    plt.savefig('ED/track2D'+option+'_'+run_day+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



def plot_2dview_hits_and_3dtracks(option=None, to_be_shown=False):

    fig, ax_leg, ax_v0, ax_v1 = template_data_view()
    
    """ unclustered hits """
    sel = 'x.cluster == -1'    
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_noise, s=marker_size, marker='o', label='Noise')

    """ clustered hits """
    sel = 'x.cluster > -1'    
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_clustered, s=marker_size, marker='o', label='Clustered')

    """ delta rays attached to track """
    sel = 'x.matched <0 and x.matched > -9999'
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_matched2, s=marker_size, marker='o', label='Delta Rays')

    """ hits attached to track """
    sel = 'x.matched >= 0'
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_matched1, s=marker_size, marker='o', label='Attached to Track')

    
    """ 2D tracks """
    ax_v0, ax_v1 = draw_all_tracks(ax_v0, ax_v1, legend='2D Track', c=color_track2d, linewidth=1)


    """ 3D tracks """
    sel = 't.matched >= 0'
    ax_v0, ax_v1 = draw_all_tracks(ax_v0, ax_v1, sel, c=color_track3d, linewidth=2, legend='3D Track')

    
    """ legend """
    ax_leg.axis('off')

    """ re-arrange the legend (lines last), and merge blue and green entries """
    h, l = ax_v0.get_legend_handles_labels()

    if(False): #len(h)==6):
        leg = ax_leg.legend([h[2], h[3], (h[4], h[5]), h[0], h[1]], [l[2], l[3], 'Hits Attached to Track (1,2)', l[0], l[1]], loc='center', ncol=5, markerscale=4, handler_map={tuple: HandlerTuple(ndivide=None)})
    else:
        leg = ax_leg.legend(*ax_v0.get_legend_handles_labels(),loc='center', ncol=6, markerscale=4, markerfirst=True)

    """ make lines in the legend bigger """
    for line in leg.get_lines():
        line.set_linewidth(3)


    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig('ED/track3D_proj'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



