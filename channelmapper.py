import config as cf
import data_containers as dc

n_ChPerConnector = 8

serhans_map = [128,127,126,125,124,123,122,121,120,119,118,117,116,115,114,113,
          112,111,110,109,108,107,106,105,104,103,102,101,100, 99, 98, 97,
          48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33,
          64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49,
          209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,
          193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,
          129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,
          145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,
          12, 14, 16,
          18, 20, 22, 24, 26, 28, 30, 32,  1,  3, 78, 80, 81, 83, 85, 87,
          89, 91, 93, 95, 65, 67, 69, 71,234,236,238,240,242,244,246,248,
          250,252,254,256,225,227,174,
          176,177,179,181,183,185,187,189,191,161,163,165,
          173,172,171,170,169,168,167,166,164,162,
          2,  4,  5,  6,  7,  8,  9,175,178,180,182,184,186,188,190,192,
          10, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 66, 68, 70, 72,
          73, 74, 75, 76, 77, 79, 82, 84, 86, 88, 90, 92, 94, 96,226,228,
          229,230,231,232,233,235,237,239,241,243,245,247,249,251,253,255]
serhans_as_map = {chan : serhans_map[chan] for chan in range(256)}
inverse_mapper = sorted(serhans_as_map, key=lambda item: serhans_as_map[item])


def check():
    if(len(dc.map_ped) > 0):
        del dc.map_ped[:]

def ChannelMapper():    
    check()
    for idaq in range(cf.n_ChanTot):
        view, vchan = DAQToAna(idaq)
        dc.map_ped.append(dc.pdmap(view,vchan))


def DAQToAna(daqch):
    conn = int(daqch/n_ChPerConnector)
    ch_femb = 7 - (daqch-conn*n_ChPerConnector) + conn*n_ChPerConnector

    ch_global = inverse_mapper[ch_femb] - 1

    if(128 > ch_global >= 64):
        # Induction 1
        chan = ch_global - 64
        view = 1
    elif(182 > ch_global >=128):
        # Induction 2
        chan = ch_global - 128
        view = 2
    elif(ch_global >= 182):
        # Shield
        chan = ch_global - 182
        view = 3
        chan = 0
    else:
        # Readout
        chan = ch_global
        view = 0

    return view, chan
