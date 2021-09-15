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

def draw_all_hits(*axes, sel='True', adc=False, charge=False, **kwargs):

    axes = list(axes)
    for iview in range(cf.physical_views):
        z = []
        if(adc==True):
            z = get_hits_adc(iview, sel)
        elif(charge==True):
            z = get_hits_charge(iview, sel)

        axes[iview] = draw_hits(pos=get_hits_pos(iview, sel),
                               time=get_hits_z(iview, sel), 
                               z=z,
                               ax=axes[iview], **kwargs)
    return tuple(axes)


def draw_tracks(pos, time, ax=None, legend="", **kwargs):

    ax = plt.gca() if ax is None else ax
    
    if(len(pos) == 0):
        return ax

    if(len(legend)>0):
        ax.plot(pos[0], time[0], label=legend, **kwargs)
        
    for tx,tz in zip(pos, time):
        ax.plot(tx,tz, **kwargs)
    
    return ax


def draw_all_tracks(*axes, sel='True', legend="", **kwargs):

    axes = list(axes)
    for iview in range(cf.physical_views):
        axes[iview] = draw_tracks(pos=get_2dtracks_pos(iview,sel),
                                 time=get_2dtracks_z(iview,sel), 
                                 ax=axes[iview],
                                 legend=legend,
                                 **kwargs)
    return tuple(axes)



def template_data_view():
    # TODO duplicated from event display?

    fig = plt.figure(figsize=(3*cf.physical_views, 4))
    gs = gridspec.GridSpec(nrows=2, ncols=cf.physical_views,
                           height_ratios=[1, 10])

    
    ax_col = fig.add_subplot(gs[0,:])
    axes = list(range(cf.physical_views))
    axes[0]  = fig.add_subplot(gs[1, 0])
    for view in range(1, cf.physical_views):
        axes[view]  = fig.add_subplot(gs[1, view], sharey=axes[0])
        plt.setp(axes[view].get_yticklabels(), visible=False)

    for view in range(cf.physical_views):
        axes[view].set_title(f'View {view} [{cf.view_type[0]}]')
        axes[view].set_ylabel('Z [cm]')
        axes[view].set_xlabel(cf.title_of_view[view])
        axes[view].set_xlim([0, cf.len_of_view[view]])
        axes[view].set_ylim([0., cf.Anode_Z])

    axes[cf.physical_views-1].yaxis.tick_right()
    axes[cf.physical_views-1].yaxis.set_label_position("right")

    plt.tight_layout()
    #plt.subplots_adjust(wspace=0.02)
    
    return (fig, ax_col, *axes)



def plot_2dview_hits(option=None, to_be_shown=False):
    if(dc.evt_list[-1].nHits[0] == 0):
        return

    fig, ax_col, axes = template_data_view()

    max_adc=50
    axes = draw_all_hits(*axes, adc=True, cmap=cmap_ed, s=marker_size, vmin=0, vmax=max_adc)
    
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
    
    fig, ax_col, axes = template_data_view()

    """ clustered hits """

    for icl in range(dc.evt_list[-1].nClusters[0]):
        sel = 'x.view == 0 and x.cluster=='+str(icl)
        axes[0] = draw_hits(pos=get_hits_pos(0, sel), time=get_hits_z(0, sel), ax=axes[0], s=marker_size, marker='o')

    for icl in range(dc.evt_list[-1].nClusters[1]):
        sel = 'x.view == 1 and x.cluster=='+str(icl)        
        axes[1] = draw_hits(pos=get_hits_pos(1, sel), time=get_hits_z(1, sel), ax=axes[1], s=marker_size, marker='o')



    """ unclustered hits """
    sel = 'x.cluster==-1'
    axes = draw_all_hits(*axes, sel=sel, c=color_noise, s=marker_size, marker='o')

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
    fig, ax_leg, *axes = template_data_view()
    
    """ all hits """
    axes = draw_all_hits(axes, c="#e6e6e6", s=marker_size, marker='o', label='Noise Hits')

    
    
    """ 2D tracks """
    axes = draw_all_tracks(ax_v0, ax_v1, linewidth=1, legend='2D Track')


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
    fig, ax_leg, *axes = template_data_view()
    
    """ unclustered hits """
    sel = 'x.cluster == -1'    
    axes = draw_all_hits(*axes, sel=sel, c=color_noise, s=marker_size, marker='o', label='Noise Hits')

    """ clustered hits """
    sel = 'x.cluster > -1'    
    axes = draw_all_hits(*axes, sel=sel, c=color_clustered, s=marker_size, marker='o', label='Hits Clustered')

    """ delta_ray hits attached to track """
    sel = 'x.matched <0 and x.matched > -9999'
    axes = draw_all_hits(*axes, sel=sel, c=color_matched2, s=marker_size, marker='o', label='Delta Rays')

    """ hits attached to track """
    sel = 'x.matched >= 0'
    axes = draw_all_hits(*axes, sel=sel, c=color_matched1, s=marker_size, marker='o', label='Hits Attached to Track')

    
    """ 2D tracks """
    axes = draw_all_tracks(*axes, legend='2D Track', c=color_track2d, linewidth=1)


    """ legend """
    ax_leg.axis('off')
    
    """ re-arrange the legend (line last), and merge blue and green entries """
    """ might not work anymore if the plotting order is changed """
    h, l = axes[0].get_legend_handles_labels()

    if(False): #len(h)==5):
        leg = ax_leg.legend([h[1], h[2], (h[3], h[4]), h[0]], [l[1], l[2], 'Hits Attached to Track (1,2)', l[0]], loc='center', ncol=4, markerscale=4, handler_map={tuple: HandlerTuple(ndivide=None)})
    else:
        """otherwise this works well """
        leg = ax_leg.legend(*axes[0].get_legend_handles_labels(),loc='center', ncol=5, markerscale=4, markerfirst=True)
    
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

    fig, ax_leg, *axes = template_data_view()
    
    """ unclustered hits """
    sel = 'x.cluster == -1'    
    axes = draw_all_hits(*axes, sel=sel, c=color_noise, s=marker_size, marker='o', label='Noise')

    """ clustered hits """
    sel = 'x.cluster > -1'    
    axes = draw_all_hits(*axes, sel=sel, c=color_clustered, s=marker_size, marker='o', label='Clustered')

    """ delta rays attached to track """
    sel = 'x.matched <0 and x.matched > -9999'
    axes = draw_all_hits(*axes, sel=sel, c=color_matched2, s=marker_size, marker='o', label='Delta Rays')

    """ hits attached to track """
    sel = 'x.matched >= 0'
    axes = draw_all_hits(*axes, sel=sel, c=color_matched1, s=marker_size, marker='o', label='Attached to Track')

    
    """ 2D tracks """
    axes = draw_all_tracks(*axes, legend='2D Track', c=color_track2d, linewidth=1)


    """ 3D tracks """
    sel = 't.matched >= 0'
    axes = draw_all_tracks(*axes, sel=sel, c=color_track3d, linewidth=2, legend='3D Track')

    
    """ legend """
    ax_leg.axis('off')

    """ re-arrange the legend (lines last), and merge blue and green entries """
    h, l = axes[0].get_legend_handles_labels()

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



