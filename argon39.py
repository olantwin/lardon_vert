import config as cf
import data_containers as dc
import lar_param as lar

import numpy as np


def recompute_charge(i, j, m):

    nchan = 5
    ntick = 10

    h0, h1 = dc.hits_list[i], dc.hits_list[j]

    if(h0.view == 1):
        h0, h1 = h1, h0

    """view 0"""
    ch_v0_min = h0.channel - nchan
    ch_v0_max = h0.channel + nchan + 1
    if(ch_v0_min < 0): ch_v0_min = 0
    if(ch_v0_max >= cf.n_ChanPerView): ch_v0_max = cf.n_ChanPerView
    tb_v0_min = h0.start - ntick
    tb_v0_max = h0.stop + ntick + 1
    if(tb_v0_min < 0): tb_v0_min = 0
    if(tb_v0_max >= cf.n_Sample): tb_v0_max = cf.n_Sample

    sum_charge = 0.
    for c in range(ch_v0_min, ch_v0_max):
        for t in range(tb_v0_min, tb_v0_max):
            sum_charge += dc.data[0,c,t]
            m[0,c,t] = True

    #print("view 0  ch : ", ch_v0_min, " ", ch_v0_max, " t: ", tb_v0_min, tb_v0_max)
    #print("hit ", h0.channel, " ", h0.max_t)
    h0.set_charge_singlehit(sum_charge)


    """view 1"""
    ch_v1_min = h1.channel - nchan
    ch_v1_max = h1.channel + nchan + 1
    if(ch_v1_min < 0): ch_v1_min = 0
    if(ch_v1_max >= cf.n_ChanPerView): ch_v1_max = cf.n_ChanPerView
    tb_v1_min = h1.start - ntick
    tb_v1_max = h1.stop + ntick + 1
    if(tb_v1_min < 0): tb_v1_min = 0
    if(tb_v1_max >= cf.n_Sample): tb_v1_max = cf.n_Sample


    sum_charge = 0.
    for c in range(ch_v1_min, ch_v1_max):
        for t in range(tb_v1_min, tb_v1_max):
            sum_charge += np.fabs(dc.data[1,c,t])
            m[1,c,t] = True

    h1.set_charge_singlehit(sum_charge)
    #print("view 1  ch : ", ch_v1_min, " ", ch_v1_max, " t: ", tb_v1_min, tb_v1_max)
    #print("hit ", h1.channel, " ", h1.min_t)
    


def argon_39_finder(time_tol):
    pair = []

    nhits = np.sum(dc.evt_list[-1].nHits)
    matched = np.zeros(nhits, dtype=bool)
    jbest = -1
    dtbest = 9999
    qbest  = -1

    int_mask = np.zeros((cf.n_View, cf.n_ChanPerView, cf.n_Sample), dtype=bool)

    for i in range(nhits):

        if(np.sum(matched) == nhits): break

        if(matched[i] == True): continue

        h0 = dc.hits_list[i]


        jbest = -1

        for j in range(nhits):
            if(i==j): continue
            if(matched[j] == True): continue

            h1 = dc.hits_list[j]
            """ single hits must be in different views """
            if(h0.view == h1.view): continue

            if(h0.view == 1):
                h0, h1 = h1, h0


            if(int_mask[h0.view, h0.channel, h0.max_t]==True): continue
            if(int_mask[h1.view, h1.channel, h1.max_t]==True): continue

            dt = np.fabs(h0.max_t - h1.min_t)
            if(dt < time_tol):
                if(jbest < 0): 
                    jbest = j
                    dtbest = dt
                    qbest = h1.min_adc
                elif(dt < dtbest and h1.min_adc < qbest):
                    jbest = j
                    dtbest = dt
                    qbest = h1.min_adc
        if(jbest >= 0):
            pair.append((i,jbest))
            matched[i] = True
            matched[jbest] = True
            recompute_charge(i,jbest, int_mask)
            #print(" TEST ", np.sum(int_mask))
            #print("1st hit v", dc.hits_list[i].view, " ", dc.hits_list[i].charge_int, " ", dc.hits_list[i].charge_sh)
            #print("2nd hit v", dc.hits_list[jbest].view, " ", dc.hits_list[jbest].charge_int, " ", dc.hits_list[jbest].charge_sh)
    return pair
