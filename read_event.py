import config as cf
import data_containers as dc
from datetime import datetime

import numpy as np
import numba as nb


def read_event(data, iev):
    idx = iev * cf.evt_byte_size
    data.seek(idx,0)

    for i in range(cf.n_femb):
        head = data.read(2*2) #2
        dt = np.dtype(np.uint32)
        dt = dt.newbyteorder('>')
        seq = np.frombuffer(head, dtype=dt)

        data.read(7*2) #8
        
        if(i == 0):
            t = data.read(3824*2)
        else:
            t += data.read(3824*2)
        for k in range(3):
            data.read(8*2)
            t += data.read(3825*2)
        data.read(8*2)
        t += data.read(851*2)



    shape_and_store( read_evt_uint12( t ))
    return 0
    

def get_file_infos(data_name, bit_size):

    timestamp = data_name[len(cf.file_prename):-4]

    ts = int(timestamp[:-2])
    tms = int(timestamp[-2:])*10

    date = datetime.fromtimestamp(ts).strftime("%d_%m_%Y")
    n_evt = int(bit_size/cf.evt_byte_size)
    return date, ts, tms, n_evt


def shape_and_store(data_raw):
    """
    Difference between arrays of times and arrays of channels
    """
    data_raw = np.reshape(data_raw, (-1, 32)).T # split into 32 channel chunks
    data_raw = np.reshape(data_raw, (32, 8, 646)) # Reshape into
    data_raw = data_raw.swapaxes(0,1)
    data_raw = np.reshape(data_raw, (256, 646))

    """reshape the array and subtract reference pedestal"""    
    for idq in range(cf.n_ChanTot):

        view, vch = dc.map_ped[idq].get_ana_chan()
        dc.data[view,vch] = data_raw[idq]

    
""" CANNOT COMPILE DUE TO ENDIAN READING PART """
@nb.jit
def read_evt_uint12_nb(data):
    """ data is taken as 12-bit and written as such, need to do bit-wise operation to retrieve values """
    """ data order is 32 adjacent channels for time t, then 32 chan for t+1 etc """
    """ given the 12-bit thing, the 32 values are encoded in 24 uint16 bytes """
    """ but each group of these 24 values have an extra byte at front """
    
    dt = np.dtype(np.uint16)
    dt = dt.newbyteorder('>')    

    tt = np.frombuffer(data, dtype=dt)
    """in big endian """
    #tt = tt.newbyteorder('>')
    """ remove the extra byte """
    tt = tt.reshape(-1,25)[:,1:].flatten()

    assert np.mod(tt.shape[0],3)==0

    out=np.empty(tt.shape[0]//3*4,dtype=np.uint16)

    for i in nb.prange(tt.shape[0]//3):
        fst_uint8=np.uint16(data[i*3])
        mid_uint8=np.uint16(data[i*3+1])
        lst_uint8=np.uint16(data[i*3+2])

        out[i*4] = ((fst_uint8 & 0x0fff) << 0)
        out[i*4+1] = ((mid_uint8 & 0x00ff) << 4) + ((fst_uint8 & 0xf000) >> 12)
        out[i*4+2] = ((lst_uint8 & 0x000f) << 8) + ((mid_uint8 & 0xff00) >> 8)
        out[i*4+3] = ((lst_uint8 & 0xfff0) >> 4)

    return out

        
def read_evt_uint12(data):
    #adapted from 
    #https://stackoverflow.com/questions/44735756/python-reading-12-bit-binary-files
    dt = np.dtype(np.uint16)
    dt = dt.newbyteorder('>')    

    tt = np.frombuffer(data, dtype=dt)    
    tt = tt.reshape(-1,25)[:,1:].flatten()

    fst, mid,  lst = np.reshape(tt, (tt.shape[0] // 3, 3)).astype(np.uint16).T

    fst_uint12 = ((fst & 0x0fff) << 0)
    snd_uint12 = ((mid & 0x00ff) << 4) + ((fst & 0xf000) >> 12) 
    thd_uint12 = ((lst & 0x000f) << 8) + ((mid&0xff00) >> 8)
    fth_uint12 = ((lst & 0xfff0) >> 4)

    return np.reshape(np.concatenate((fst_uint12[:, None], snd_uint12[:, None], thd_uint12[:,None], fth_uint12[:,None]), axis=1), 4 * fst_uint12.shape[0])


