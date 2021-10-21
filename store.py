import config as cf
from tables import *
import numpy as np
import data_containers as dc

class Infos(IsDescription):
    date         = StringCol(10)
    run          = UInt8Col()
    commit        = UInt8Col()
    nviews        = UInt8Col()
    #time_s       = UInt32Col()
    #time_ms      = UInt8Col()
    nEventTot     = UInt32Col()
    nEventMu      = UInt32Col()
    nEventSH      = UInt32Col()
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
    nSingleHits   = UInt32Col()


class Pedestal(IsDescription):
    ini_mean_v0   = Float16Col(shape=(cf.n_ChanPerView))
    ini_rms_v0    = Float16Col(shape=(cf.n_ChanPerView))
    ini_mean_v1   = Float16Col(shape=(cf.n_ChanPerView))
    ini_rms_v1    = Float16Col(shape=(cf.n_ChanPerView))
    ini_mean_v2   = Float16Col(shape=(cf.n_ChanPerView))
    ini_rms_v2    = Float16Col(shape=(cf.n_ChanPerView))

    fin_mean_v0   = Float16Col(shape=(cf.n_ChanPerView))
    fin_rms_v0    = Float16Col(shape=(cf.n_ChanPerView))
    fin_mean_v1   = Float16Col(shape=(cf.n_ChanPerView))
    fin_rms_v1    = Float16Col(shape=(cf.n_ChanPerView))
    fin_mean_v2   = Float16Col(shape=(cf.n_ChanPerView))
    fin_rms_v2    = Float16Col(shape=(cf.n_ChanPerView))


class SingleHits(IsDescription):
    channel_v0 = UInt16Col()
    channel_v1 = UInt16Col()
    channel_v2 = UInt16Col()

    tdc_max_v0 = UInt16Col()
    tdc_max_v1 = UInt16Col()
    tdc_max_v2 = UInt16Col()
    tdc_min_v0 = UInt16Col()
    tdc_min_v1 = UInt16Col()
    tdc_min_v2 = UInt16Col()

    x          = Float16Col()
    y          = Float16Col()
    z          = Float16Col()

    adc_max_v0 = Float16Col()
    adc_max_v1 = Float16Col()
    adc_max_v2 = Float16Col()
    adc_min_v0 = Float16Col()
    adc_min_v1 = Float16Col()
    adc_min_v2 = Float16Col()
    
    charge_int_v0 = Float16Col()
    charge_int_v1 = Float16Col()
    charge_int_v2 = Float16Col()
    charge_max_v0 = Float16Col()
    charge_max_v1 = Float16Col()
    charge_max_v2 = Float16Col()
    charge_min_v0 = Float16Col()
    charge_min_v1 = Float16Col()
    charge_min_v2 = Float16Col()
    charge_pv_v0  = Float16Col()
    charge_pv_v1  = Float16Col()
    charge_pv_v2  = Float16Col()

    charge_all_v0 = Float16Col()
    charge_all_v1 = Float16Col()
    charge_all_v2 = Float16Col()


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

class FFT(IsDescription):
    PS_v0 = Float16Col(shape=(64, 324))
    PS_v1 = Float16Col(shape=(64, 324))
    PS_v2 = Float16Col(shape=(64, 324))


def create_output(h5file):
    h5file.create_table("/", 'infos', Infos, "General Reconstruction Informations")
    h5file.create_table("/", 'events', Event, "Event Reconstruction Informations")
    h5file.create_table("/", 'pedestals_fake', Pedestal, "Pedestals Empty Events")
    h5file.create_table("/", 'pedestals_muon', Pedestal, "Pedestals Data Events")
    h5file.create_table("/", 'hits', Hits, "Pedestals")

    h5file.create_table("/", 'tracks2D', Tracks2D, "Tracks 2D")
    h5file.create_vlarray("/", 'trk2D_v0', Float32Atom(shape=(6)), "2D Track V0 (x, z, qint, qmax, qmin, qpv)")  
    h5file.create_vlarray("/", 'trk2D_v1', Float32Atom(shape=(6)), "2D Track V1 (x, z, qint, qmax, qmin, qpv)")
    h5file.create_vlarray("/", 'trk2D_v2', Float32Atom(shape=(6)), "2D Track V2 (x, z, qint, qmax, qmin, qpv)")

    h5file.create_table("/", 'tracks3D', Tracks3D, "Tracks 3D")
    h5file.create_vlarray("/", 'trk3D_v0', Float32Atom(shape=(8)), "3D Track V0 (x, y, z, qint, qmax, qmin, qpv, ds)")
    h5file.create_vlarray("/", 'trk3D_v1', Float32Atom(shape=(8)), "3D Track V1 (x, y, z, qint, qmax, qmin, qpv, ds)")
    h5file.create_vlarray("/", 'trk3D_v2', Float32Atom(shape=(8)), "3D Track V2 (x, y, z, qint, qmax, qmin, qpv, ds)")

    h5file.create_table("/", 'singleHits', SingleHits, "Single 3D Hits")
    h5file.create_table("/", 'fft_fake', FFT, "FFT Empty Events")
    h5file.create_table("/", 'fft_muon', FFT, "FFT Data Events")



def store_infos(h5file, date, run, nevtTot, nevtMu, nevtSH, time):
    inf = h5file.root.infos.row
    inf['date']    = date
    inf['run']     = run
    #inf['time_s']  = t_s
    #inf['time_ms'] = t_ms
    inf['nEventTot']  = nevtTot
    inf['nEventMu']   = nevtMu
    inf['nEventSH']   = nevtSH
    inf['process_date'] = time
    
    inf.append()

    

def store_event(h5file):
    evt = h5file.root.events.row
    evt['evt_nb_global'] = dc.evt_list[-1].evt_nb_glob
    evt['evt_nb_local']  = dc.evt_list[-1].evt_nb
    evt['time_s']        = dc.evt_list[-1].time_s
    evt['time_ms']       = dc.evt_list[-1].time_ms
    evt['nHits']         = dc.evt_list[-1].nHits
    evt['nClusters']     = dc.evt_list[-1].nClusters
    evt['nTracks2D']     = dc.evt_list[-1].nTracks2D
    evt['nTracks3D']     = dc.evt_list[-1].nTracks3D
    evt['nSingleHits']   = dc.evt_list[-1].nSingleHits
    evt.append()

def store_pedestals_data(h5file):
    ped = h5file.root.pedestals_muon.row
    
    ini_mean = np.zeros((cf.n_View, cf.n_ChanPerView))
    ini_rms  = np.zeros((cf.n_View, cf.n_ChanPerView))
    fin_mean = np.zeros((cf.n_View, cf.n_ChanPerView))
    fin_rms  = np.zeros((cf.n_View, cf.n_ChanPerView))


    for x in dc.map_ped:
        ini_mean[x.view, x.vchan] = x.raw_ped
        ini_rms[x.view, x.vchan]  = x.raw_rms
        fin_mean[x.view, x.vchan] = x.evt_ped
        fin_mean[x.view, x.vchan] = x.evt_ped


    for view in range(cf.physical_views):
        ped[f'ini_mean_v{view}'] = ini_mean[view]
        ped[f'ini_rms_v{view}']  = ini_rms[view]
        ped[f'fin_mean_v{view}'] = fin_mean[view]
        ped[f'fin_rms_v{view}']  = fin_rms[view]

    ped.append()

def store_pedestals_fake(h5file, ini_mean, ini_rms, fin_mean, fin_rms):
    ped = h5file.root.pedestals_fake.row
    

    for view in range(cf.physical_views):
        ped[f'ini_mean_v{view}'] = ini_mean[view]
        ped[f'ini_rms_v{view}']  = ini_rms[view]
        ped[f'fin_mean_v{view}'] = fin_mean[view]
        ped[f'fin_rms_v{view}']  = fin_rms[view]

    ped.append()

def store_hits(h5file):
    hit = h5file.root.hits.row
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



def store_tracks2D(h5file):
    t2d = h5file.root.tracks2D.row
    views = h5file.root.trk2D_v0, h5file.root.trk2D_v1, h5file.root.trk2D_v2

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

        t2d.append()

        views[t.view].append(pts)

def store_tracks3D(h5file):
    t3d = h5file.root.tracks3D.row
    v0  = h5file.root.trk3D_v0
    v1  = h5file.root.trk3D_v1


    for t in dc.tracks3D_list:
        t3d['x_ini']     = t.ini_x
        t3d['y_ini']     = t.ini_y
        t3d['z_ini']     = t.ini_z
        t3d['x_end']     = t.end_x
        t3d['y_end']     = t.end_y
        t3d['z_end']     = t.end_z

        t3d['chi2']      = t.chi2
        # TODO add sensible values for third view
        t3d['nHits']     = [t.nHits_v0, t.nHits_v1, 0, 0]
        t3d['len_straight'] = [t.len_straight_v0, t.len_straight_v1, 0, 0]
        t3d['len_path']     = [t.len_path_v0, t.len_path_v1, 0, 0]

        t3d['total_charge_int'] = [t.tot_charge_int_v0, t.tot_charge_int_v1, 0, 0]
        t3d['total_charge_max'] = [t.tot_charge_max_v0, t.tot_charge_max_v1, 0, 0]
        t3d['total_charge_min'] = [t.tot_charge_min_v0, t.tot_charge_min_v1, 0, 0]
        t3d['total_charge_pv']  = [t.tot_charge_pv_v0,  t.tot_charge_pv_v1, 0, 0]


        t3d['theta_ini'] = t.ini_theta
        t3d['theta_end'] = t.end_theta
        t3d['phi_ini']   = t.ini_phi
        t3d['phi_end']   = t.end_phi

        t3d['z0_corr']    = t.z0_corr
        t3d['t0_corr']    = t.t0_corr

        pts_v0 = [[p[0], p[1], p[2], q, r, s, t, u] for p,q,r,s,t,u in zip(t.path_v0,t.dQds_int_v0, t.dQds_max_v0, t.dQds_min_v0, t.dQds_pv_v0, t.ds_v0)]
        pts_v1 = [[p[0], p[1], p[2], q, r, s, t, u] for p,q,r,s,t,u in zip(t.path_v1, t.dQds_int_v1, t.dQds_max_v1, t.dQds_min_v1, t.dQds_pv_v1, t.ds_v1)]
        
        t3d.append()
        v0.append(pts_v0)
        v1.append(pts_v1)

def store_single_hits(h5file, i, j):
    hit = h5file.root.singleHits.row


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

    hit['charge_all_v0']  = hv0.charge_sh
    hit['charge_all_v1']  = hv1.charge_sh

    hit.append()


def store_fft_data(h5file, ps):
    fft = h5file.root.fft_muon.row
    fft['PS_v0'] = ps[0]
    fft['PS_v1'] = ps[1]
    fft.append()

def store_fft_fake(h5file, ps):
    fft = h5file.root.fft_fake.row
    fft['PS_v0'] = ps[0]
    fft['PS_v1'] = ps[1]
    fft.append()

