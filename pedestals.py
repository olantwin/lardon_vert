import config as cf
import data_containers as dc

import numba as nb
import numpy as np
import bottleneck as bn



def ini_pedestal():
    ped_mean = bn.move_mean(dc.data, window=50, axis=-1)
    ped_std  = bn.move_std(dc.data, window=50, axis=-1)
    ped_std  = np.nan_to_num(ped_std, nan=9999.)

    idx = np.argmin(ped_std, axis=-1)

    dc.ped_mean = np.take_along_axis(ped_mean, np.expand_dims(idx, axis=-1), axis=-1)[:,:,0]
    dc.ped_rms  = np.take_along_axis(ped_std, np.expand_dims(idx, axis=-1), axis=-1)[:,:,0]




@nb.jit(nopython = True)
def compute_pedestal_RMS_nb(data, mask):
    """ do not remove the @jit above, the code would then take ~40 s to run """
    shape = data.shape

    """ to avoid cases where the rms goes to 0"""
    min_val = 1e-5

    res   = np.zeros(shape[:-1])
    for idx,v in np.ndenumerate(res):
        ch = data[idx]
        ma  = mask[idx]
        """ use the assumed mean method """
        K = ch[0]
        n, Ex, Ex2, = 0., 0., 0.
        for x,v in zip(ch,ma):
            if(v == True):
                n += 1
                Ex += x-K
                Ex2 += (x-K)*(x-K)

        """cut at min. 10 pts to compute a proper RMS"""
        if( n < 10 ):
            res[idx] = -1.
        else:
            val = np.sqrt((Ex2 - (Ex*Ex)/n)/(n-1))
            res[idx] = min_val if val < min_val else val
    return res



def compute_pedestal_RMS():
    """ the numpy way is slower and cannot handle well dead channels """
    dc.ped_rms = compute_pedestal_RMS_nb(dc.data, dc.mask)


def compute_pedestal_mean():
    """ the numpy way may be faster but do not handle dead channels """

    with np.errstate(divide='ignore', invalid='ignore'):
        dc.ped_mean = np.einsum('ijk,ijk->ij', dc.data, dc.mask)/dc.mask.sum(axis=2)
        """require at least 3 points to take into account the mean"""
        dc.ped_mean[dc.mask.sum(axis=2) < 3] = -999.



def store_raw_ped_rms():
    """ store the raw pedestal and rms """
    for i in range(cf.n_ChanTot):
        view, ch = dc.map_ped[i].get_ana_chan()        
        dc.map_ped[i].set_raw_pedestal(dc.ped_mean[view,ch], dc.ped_rms[view,ch])
        dc.ped_mean[view, ch] = 0.

def subtract_ped_mean():
    for i in range(cf.n_ChanTot):
        view, ch = dc.map_ped[i].get_ana_chan()        
        if(dc.ped_mean[view,ch] != -999.):
            dc.data[view, ch] -= dc.ped_mean[view,ch]

def compute_ped():
    compute_pedestal_mean()
    compute_pedestal_RMS()

def compute_and_subtract_ped():
    compute_ped()
    subtract_ped_mean()

def store_final_ped_rms():
    for i in range(cf.n_ChanTot):
        view, ch = dc.map_ped[i].get_ana_chan()
        dc.map_ped[i].set_evt_pedestal(dc.ped_mean[view,ch], dc.ped_rms[view,ch])

