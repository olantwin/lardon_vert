import config as cf
from tables import *
import numpy as np
import data_containers as dc

class Infos(IsDescription):
    date         = StringCol(10)
    run          = UInt8Col()
    #time_s       = UInt32Col()
    #time_ms      = UInt8Col()
    nEventTot     = UInt32Col()
    nEventMu      = UInt32Col()
    process_date = UInt32Col()

class Event(IsDescription):
    evt_nb_global = UInt32Col()
    evt_nb_local  = UInt32Col()
    time_s        = UInt32Col()
    time_ms       = UInt32Col()
    nHits         = UInt32Col(shape=(cf.n_View))
    nClusters     = UInt32Col(shape=(cf.n_View))
    nTracks2D     = UInt32Col(shape=(cf.n_View))
    nTracks3D     = UInt32Col()

class Pedestal(IsDescription):
    view      = UInt8Col()
    channel   = UInt16Col()
    raw_mean  = Float16Col()
    raw_rms   = Float16Col()
    filt_mean = Float16Col()
    filt_rms  = Float16Col()


class AllPedestal(IsDescription):
    mean_v0   = Float16Col(shape=(cf.n_ChanPerView))
    rms_v0    = Float16Col(shape=(cf.n_ChanPerView))
    mean_v1   = Float16Col(shape=(cf.n_ChanPerView))
    rms_v1    = Float16Col(shape=(cf.n_ChanPerView))

class SingleHits(IsDescription):
    channel_v0 = UInt16Col()
    channel_v1 = UInt16Col()

    tdc_max_v0 = UInt16Col()
    tdc_max_v1 = UInt16Col()
    tdc_min_v0 = UInt16Col()
    tdc_min_v1 = UInt16Col()

    x          = Float16Col()
    y          = Float16Col()
    z          = Float16Col()

    adc_max_v0 = Float16Col()
    adc_max_v1 = Float16Col()
    adc_min_v0 = Float16Col()
    adc_min_v1 = Float16Col()
    
    charge_int_v0 = Float16Col()
    charge_int_v1 = Float16Col()
    charge_max_v0 = Float16Col()
    charge_max_v1 = Float16Col()
    charge_min_v0 = Float16Col()
    charge_min_v1 = Float16Col()
    charge_pv_v0  = Float16Col()
    charge_pv_v1  = Float16Col()

class Hits(IsDescription):
    view    = UInt8Col()
    channel = UInt16Col()
    tdc_max = UInt16Col()
    tdc_min = UInt16Col()
    z       = Float16Col()
    x       = Float16Col()
    dt      = Float16Col()

    adc_max  = Float16Col()
    adc_min  = Float16Col()

    charge_int  = Float16Col()
    charge_max = Float16Col()
    charge_min = Float16Col()
    charge_pv  = Float16Col()
    cluster = Int16Col()

    ped_bef  = Float16Col()
    ped_aft  = Float16Col()

class Tracks2D(IsDescription):
    view    = UInt8Col()
    pos_ini = Float16Col()
    pos_end = Float16Col()
    z_ini   = Float16Col()
    z_end   = Float16Col()
    nHits   = UInt16Col()
    chi2    = Float16Col()

    slope_ini = Float16Col()
    slope_end = Float16Col()

    len_straight = Float16Col()
    len_path     = Float16Col()

    total_charge_int = Float16Col()
    total_charge_max = Float16Col()
    total_charge_min = Float16Col()
    total_charge_pv  = Float16Col()

class Tracks3D(IsDescription):
    x_ini   = Float16Col()
    y_ini   = Float16Col()
    z_ini   = Float16Col()
    x_end   = Float16Col()
    y_end   = Float16Col()
    z_end   = Float16Col()
    chi2    = Float16Col()

    theta_ini = Float16Col()
    theta_end = Float16Col()
    phi_ini   = Float16Col()
    phi_end   = Float16Col()

    nHits        = UInt16Col(shape=(cf.n_View))
    len_straight = Float16Col(shape=(cf.n_View))
    len_path     = Float16Col(shape=(cf.n_View))

    total_charge_int   = Float16Col(shape=(cf.n_View))
    total_charge_max   = Float16Col(shape=(cf.n_View))
    total_charge_min   = Float16Col(shape=(cf.n_View))
    total_charge_pv   = Float16Col(shape=(cf.n_View))

    z0_corr = Float16Col()
    t0_corr = Float16Col()
    
def new_event(h5file, event_nb):
    return h5file.create_group("/", 'event_'+str(event_nb), 'Event '+str(event_nb))    

def store_infos(h5file, date, run, nevtTot, nevtMu, time):
    table = h5file.create_table("/", 'infos', Infos, 'Infos')
    inf = table.row

    inf['date']    = date
    inf['run']     = run
    #inf['time_s']  = t_s
    #inf['time_ms'] = t_ms
    inf['nEventTot']  = nevtTot
    inf['nEventMu']   = nevtMu
    inf['process_date'] = time

    inf.append()
    table.flush()

def store_event(h5file, group):
    table = h5file.create_table(group, 'event', Event, "Event")

    evt = table.row

    evt['evt_nb_global'] = dc.evt_list[-1].evt_nb_glob
    evt['evt_nb_local']  = dc.evt_list[-1].evt_nb
    evt['time_s']        = dc.evt_list[-1].time_s
    evt['time_ms']       = dc.evt_list[-1].time_ms
    evt['nHits']         = dc.evt_list[-1].nHits
    evt['nClusters']     = dc.evt_list[-1].nClusters
    evt['nTracks2D']     = dc.evt_list[-1].nTracks2D
    evt['nTracks3D']     = dc.evt_list[-1].nTracks3D
    evt.append()
    table.flush()

def create_all_pedestals(h5file):
    table = h5file.create_table("/", 'AllPedestal', AllPedestal, 'all pedestals')

def store_all_pedestals(h5file):
    table = h5file.root.AllPedestal
    ped = table.row

    pedv0 = np.zeros(cf.n_ChanPerView)
    pedv1 = np.zeros(cf.n_ChanPerView)
    rmsv0 = np.zeros(cf.n_ChanPerView)
    rmsv1 = np.zeros(cf.n_ChanPerView)
    for x in dc.map_ped:
        v, ch = x.get_ana_chan()
        if(v==0):
            pedv0[ch] = x.raw_ped
            rmsv0[ch] = x.raw_rms
        else:
            pedv1[ch] = x.raw_ped
            rmsv1[ch] = x.raw_rms


    ped['mean_v0'] = pedv0
    ped['rms_v0']  = rmsv0
    ped['mean_v1'] = pedv1
    ped['rms_v1']  = rmsv1
    ped.append()
    table.flush()
    


def create_single_hits(h5file):
    table = h5file.create_table("/", 'SingleHits', SingleHits, 'single hits')

def store_single_hits(h5file, i, j):
    table = h5file.root.SingleHits

    hit = table.row
    hv0 = dc.hits_list[i]
    hv1 = dc.hits_list[j]

        
    hit['channel_v0'] = hv0.channel
    hit['channel_v1'] = hv1.channel
    hit['tdc_max_v0'] = hv0.max_t
    hit['tdc_max_v1'] = hv1.max_t
    hit['tdc_min_v0'] = hv0.min_t
    hit['tdc_min_v1'] = hv1.min_t

    hit['x']       = hv0.X
    hit['y']       = hv1.X
    hit['z']       = 0.5*(hv0.Z + hv1.Z)

    hit['adc_max_v0'] = hv0.max_adc
    hit['adc_min_v0'] = hv0.min_adc
    hit['adc_max_v1'] = hv1.max_adc
    hit['adc_min_v1'] = hv1.min_adc
        
    hit['charge_int_v0']  = hv0.charge_int
    hit['charge_max_v0']  = hv0.charge_max
    hit['charge_min_v0']  = hv0.charge_min
    hit['charge_pv_v0']   = hv0.charge_pv

    hit['charge_int_v1']  = hv1.charge_int
    hit['charge_max_v1']  = hv1.charge_max
    hit['charge_min_v1']  = hv1.charge_min
    hit['charge_pv_v1']   = hv1.charge_pv


    hit.append()
    table.flush()
    


def store_pedestal(h5file, group):
    table = h5file.create_table(group, 'pedestals', Pedestal, 'Pedestals')

    ped = table.row

    for x in dc.map_ped:
        ped['view']      = x.view
        ped['channel']   = x.vchan
        ped['raw_mean']  = x.raw_ped
        ped['raw_rms']   = x.raw_rms
        ped['filt_mean'] = x.evt_ped
        ped['filt_rms']  = x.evt_rms
        ped.append()
    table.flush()


def store_hits(h5file, group):
    table = h5file.create_table(group, 'hits', Hits, 'Hits')

    hit = table.row
    for h in dc.hits_list:
        hit['view']    = h.view
        hit['channel'] = h.channel
        hit['tdc_max'] = h.max_t
        hit['tdc_min'] = h.min_t
        hit['z']       = h.Z
        hit['x']       = h.X
        hit['dt']      = (h.stop - h.start)*cf.n_Sampling

        hit['adc_max'] = h.max_adc
        hit['adc_min'] = h.min_adc

        hit['charge_int']  = h.charge_int
        hit['charge_max']  = h.charge_max
        hit['charge_min']  = h.charge_min
        hit['charge_pv']   = h.charge_pv

        hit['cluster'] = h.cluster

        hit['ped_bef'] = h.ped_bef
        hit['ped_aft'] = h.ped_aft
        
        hit.append()
    table.flush()


def store_tracks2D(h5file, group):    
    table = h5file.create_table(group, 'tracks2D', Tracks2D, "Tracks 2D")       

    t2d_hits_v0 = h5file.create_group(group, 'tracks2D_v0', 'Tracks 2D View0')
    t2d_hits_v1 = h5file.create_group(group, 'tracks2D_v1', 'Tracks 2D View1')

    t2d = table.row
    i_v0 = 0
    i_v1 = 0

    for t in dc.tracks2D_list:
        t2d['view']      = t.view
        t2d['pos_ini']   = t.path[0][0]
        t2d['z_ini']     = t.path[0][1]
        t2d['pos_end']   = t.path[-1][0]
        t2d['z_end']     = t.path[-1][1]

        t2d['nHits']     = t.nHits
        t2d['chi2']      = t.chi2

        t2d['slope_ini'] = t.ini_slope
        t2d['slope_end'] = t.end_slope

        t2d['len_straight'] = t.len_straight
        t2d['len_path']     = t.len_path

        t2d['total_charge_int'] = t.tot_charge_int
        t2d['total_charge_max'] = t.tot_charge_max
        t2d['total_charge_min'] = t.tot_charge_min
        t2d['total_charge_pv'] = t.tot_charge_pv

        pts = [[p[0], p[1], q[0], q[1], q[2], q[3]] for p,q in zip(t.path,t.dQ)]

        if(t.view==0):
            h5file.create_array(t2d_hits_v0, 'track_%i'%(i_v0), np.asarray(pts), 'track hits')
            i_v0 += 1
        else:
            h5file.create_array(t2d_hits_v1, 'track_%i'%(i_v1), np.asarray(pts), 'track hits')
            i_v1 += 1

        t2d.append()
    table.flush()



def store_tracks3D(h5file, group):
    table = h5file.create_table(group, 'tracks3D', Tracks3D, "Tracks 3D")       

    t3d_hits_v0 = h5file.create_group(group, 'tracks3D_v0', 'Tracks 3D View0')
    t3d_hits_v1 = h5file.create_group(group, 'tracks3D_v1', 'Tracks 3D View1')

    t3d = table.row
    i = 0

    for t in dc.tracks3D_list:
        t3d['x_ini']     = t.ini_x
        t3d['y_ini']     = t.ini_y
        t3d['z_ini']     = t.ini_z
        t3d['x_end']     = t.end_x
        t3d['y_end']     = t.end_y
        t3d['z_end']     = t.end_z

        t3d['chi2']      = t.chi2
        t3d['nHits']     = [t.nHits_v0, t.nHits_v1]
        t3d['len_straight'] = [t.len_straight_v0, t.len_straight_v1]
        t3d['len_path']     = [t.len_path_v0, t.len_path_v1]

        t3d['total_charge_int'] = [t.tot_charge_int_v0, t.tot_charge_int_v1]
        t3d['total_charge_max'] = [t.tot_charge_max_v0, t.tot_charge_max_v1]
        t3d['total_charge_min'] = [t.tot_charge_min_v0, t.tot_charge_min_v1]
        t3d['total_charge_pv']  = [t.tot_charge_pv_v0,  t.tot_charge_pv_v1]


        t3d['theta_ini'] = t.ini_theta
        t3d['theta_end'] = t.end_theta
        t3d['phi_ini']   = t.ini_phi
        t3d['phi_end']   = t.end_phi

        t3d['z0_corr']    = t.z0_corr
        t3d['t0_corr']    = t.t0_corr

        pts_v0 = [[p[0], p[1], p[2], q, r, s, t, u] for p,q,r,s,t,u in zip(t.path_v0,t.dQds_int_v0, t.dQds_max_v0, t.dQds_min_v0, t.dQds_pv_v0, t.ds_v0)]
        pts_v1 = [[p[0], p[1], p[2], q, r, s, t, u] for p,q,r,s,t,u in zip(t.path_v1, t.dQds_int_v1, t.dQds_max_v1, t.dQds_min_v1, t.dQds_pv_v1, t.ds_v1)]
        

        h5file.create_array(t3d_hits_v0, 'track_%i'%(i), np.asarray(pts_v0), 'track hits')
        h5file.create_array(t3d_hits_v1, 'track_%i'%(i), np.asarray(pts_v1), 'track hits')
        i += 1

        t3d.append()
    table.flush()
