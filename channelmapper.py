import config as cf
import data_containers as dc

n_ChPerConnector = 8

mapper = [65, 66, 67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,
          85,86,87,88,89,90,91,92,93,94,95,96,113,114,115,116,117,118,119,
          120,121,122,123,124,125,126,127,128,97,98,99,100,101,102,103,
          104,105,106,107,108,109,110,111,112,16,15,14,13,12,11,10,9,8,
          7,6,5,4,3,2,1,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,
          64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49,48,47,46,45,44,
          43,42,41,40,39,38,37,36,35,34,33]


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

    ch_global = mapper[ch_femb] - 1

    if(ch_global >= 64):
        chan = ch_global - 64
        view = 1
    else:
        chan = ch_global
        view = 0

    return view, chan
