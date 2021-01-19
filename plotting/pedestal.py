import config as cf
import data_containers as dc

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import collections  as mc
from matplotlib.legend_handler import HandlerTuple
import itertools as itr
import math
import colorcet as cc

col_v0      = '#FBA120' #c'
col_v1      = '#435497' #'orange'

ped_min = -3
ped_max = 3
rms_min = 0.
rms_max = 3.5

def plot_event_ped_and_rms(option=None, to_be_shown=False):

    fig = plt.figure(figsize=(9,4))
    gs = gridspec.GridSpec(nrows=2, 
                           ncols=2)


    ax_ped_v0 = fig.add_subplot(gs[0,0])
    ax_ped_v1 = fig.add_subplot(gs[0,1], sharey=ax_ped_v0)
    ax_rms_v0  = fig.add_subplot(gs[1, 0], sharex=ax_ped_v0)
    ax_rms_v1  = fig.add_subplot(gs[1, 1], sharey=ax_rms_v0, sharex=ax_ped_v1)
    

    """ draw data """
    ch = np.arange(64)

    """ Mean Pedestal """
    ax_ped_v0.scatter(ch, dc.ped_mean[0,:], s=3, marker='o', color=col_v0)
    ylow, yup = ax_ped_v0.get_ylim()
    ax_ped_v0.plot([32, 32], [ylow, yup], ls='--', c='k')
    ax_ped_v0.set_title('View 0 ' + '['+cf.view_type[0]+']')
    ax_ped_v0.set_ylabel('Mean [ADC]')
    ax_ped_v0.set_ybound(lower=ped_min, upper=ped_max)

    ax_ped_v1.scatter(ch, dc.ped_mean[1,:], s=3, marker='o', color=col_v1)
    ax_ped_v1.set_title('View 1 ' + '['+cf.view_type[1]+']')
    ax_ped_v1.set_ylabel('Mean [ADC]')
    ax_ped_v1.yaxis.tick_right()
    ax_ped_v1.yaxis.set_label_position("right")
    ax_ped_v1.set_ybound(lower=ped_min, upper=ped_max)

    ax_ped_v0.tick_params(labelbottom=False)
    ax_ped_v1.tick_params(labelbottom=False)

    """ RMS Pedestal """
    ax_rms_v0.scatter(ch, dc.ped_rms[0,:], s=3, marker='o', color=col_v0)
    #ax_rms_v0.set_title('View 0 ' + '['+cf.view_type[0]+']')
    ylow, yup = ax_rms_v0.get_ylim()
    ax_rms_v0.plot([32, 32], [ylow, yup], ls='--', c='k')

    ax_rms_v0.set_ylabel('RMS [ADC]')
    ax_rms_v0.set_ybound(lower=rms_min, upper=rms_max)

    ax_rms_v1.scatter(ch, dc.ped_rms[1,:], s=3, marker='o', color=col_v1)
    ax_rms_v1.set_ylabel('RMS [ADC]')
    ax_rms_v1.yaxis.tick_right()
    ax_rms_v1.yaxis.set_label_position("right")
    ax_rms_v1.set_ybound(lower=rms_min, upper=rms_max)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.02, hspace=0.06)


    if(option):
        option = "_"+option
    else:
        option = ""

    run_day = str(dc.evt_list[-1].date)
    evt_nb = str(dc.evt_list[-1].evt_nb)

    plt.savefig('ED/ped_'+option+'_'+run_day+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()


