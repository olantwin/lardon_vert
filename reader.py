import sys
import os
import glob

import numpy as np
import time 
import tables as tables

import config as cf
import data_containers as dc

import pedestals as ped
import channelmapper as cmap
import read_event as read
#import plot_event as plot_ev
import noise_filter as noise
import hitfinder as hf
import clustering as clus
import track_2d as trk2d
import track_3d as trk3d
import store as store
#import read_mc as mc
import argon39 as ar

import plotting as plot


def need_help():
    print("Usage: python reader.py ")
    print(" FOR ONLY ONE BINARY FILE: ")
    print(" -file <file name> ")
    print(" TO RECONSTRUCT A WHOLE DAY/RUN: ")
    print(" -day <MM_DD_YYYY> -run <number of runXXtri folder>")
    print(" -n   <number of event to process in a file>  [default (or -1) is all]")
    print(" -ntot   <Total number of event to process>  [default (or -1) is all]")
    print(" -out <output name optn>")
    print(" -v :: activate verbose option ")
    print(" -h print this message")
    
    sys.exit()
    

if len(sys.argv) == 1:
    need_help()
else:
    for index, arg in enumerate(sys.argv):
        if arg in ['-h'] :
            need_help()
            

""" Reconstruction parameters """
""" not usefull in 50 L data """
lowpasscut       = 0.5 #0.06 #0.1 #MHz    
freqlines        = []#0.0234, 0.0625, 0.0700] #in MHz


signal_thresh    = [4., 3., 3.5]
signal_thresh_2  = [2.5, 1.5, 2.]

""" Different values for the ADC ROI search 
[collection, Induction positive, Induction Negative] """
adc_thresh       = [10., 10., -20.]
coherent_groups  = [64]

outname_option = ""
nevent_per_file = -1 
ntotevent = -1
nevent_mu = 0
file_in = ""
day = ""
run = -1

verbose = False

for index, arg in enumerate(sys.argv):
    if arg in ['-file'] and len(sys.argv) > index + 1:
        file_in = sys.argv[index + 1]
    elif arg in ['-day'] and len(sys.argv) > index + 1:
        day = sys.argv[index + 1]
    elif arg in ['-run'] and len(sys.argv) > index + 1:
        run = int(sys.argv[index + 1])
    elif arg in ['-n'] and len(sys.argv) > index + 1:
        nevent_per_file = int(sys.argv[index + 1])
    elif arg in ['-ntot'] and len(sys.argv) > index + 1:
        ntotevent = int(sys.argv[index + 1])
    elif arg in ['-out'] and len(sys.argv) > index + 1:
        outname_option = sys.argv[index + 1]
    elif arg in ['-v']:
        verbose = True
    



if(len(file_in) == 0 and len(day)==0 and run < 0):
    need_help()
data_list = []

if(len(file_in) > 0 and len(day)==0 and run < 0):
    name_in = cf.data_path + "/" + file_in
    if(os.path.exists(name_in) is False):
        print(" ERROR ! there is no ", name_in, " ! ")
        sys.exit()
    outname = "reco_tracks"
    data_list.append(name_in)

elif(len(file_in) == 0 and len(day) > 0 and run > 0):
    name_in = cf.data_path + "/Rawdata_" + day + "/run%02dtri"%run
    if(os.path.exists(name_in) is False):
        print(" ERROR ! there is no ", name_in, " ! ")
        sys.exit()
    outname = day + "_run" + str(run) + "_reco_tracks"
    print("--> ", outname)
    data_list.extend(glob.glob(name_in+"/*.bin"))
else:
    need_help()



if(outname_option):
    outname = outname + "_" + outname_option
        
    
name_out = cf.store_path + "/" +  outname + ".h5"
print(name_out)

output = tables.open_file(name_out, mode="w", title="Reconstruction Output")
store.create_output(output)

#name_out_sh = cf.store_path + "/" +  outname + "_singlehits.h5"
#output_sh = tables.open_file(name_out_sh, mode="w", title="Reconstruction Output")
#store.create_single_hits(output_sh)

#name_out_ped = cf.store_path + "/" +  outname + "_pedestals.h5"
#output_ped = tables.open_file(name_out_ped, mode="w", title="Pedestals")
#store.create_all_pedestals(output_ped)

""" Build DAQ Channel <-> Analysis Channel correspondance """
cmap.ChannelMapper()

""" Set the eventual broken channel to remove from the reconstruction """
for ibrok in cf.daq_broken_channels:
    view, vch = cmap.DAQToAna(ibrok)
    dc.alive_chan[view, vch, : ] = False


for view, vch in cf.crp_broken_channels:    
    dc.alive_chan[view, vch, : ] = False

""" set the plotting style """
plot.set_style()

tstart = time.time()
ievent_tot = 0


nfake = 0
ps_fake = np.zeros((2,64,324))
ps_muon = np.zeros((2,64,324))

n_events_single_hits = 0

print(" =^=^=^=^=^=^=^=^=^=^=^=^ ")
print(" RUNNING OVER ", len(data_list), " FILES ")
print(" =^=^=^=^=^=^=^=^=^=^=^=^ ")



for fbin in data_list:
    data = open(fbin, "rb")

    """ Data binary files have no Run Header """
    fsize = os.path.getsize(fbin)
    date, ts, tms, nb_evt = read.get_file_infos(fbin[fbin.rfind("/")+1:], fsize)

    if(verbose):
        print(date, " ", ts, " ", tms, " ", nb_evt)

    nevent = nevent_per_file
    if(nevent > nb_evt or nevent < 0):
        nevent = nb_evt

    if(ntotevent > 0 and ievent_tot >= ntotevent):
        #print("bye bye")
        break

    print(" File ", fbin[fbin.rfind("/")+1:])

    if(verbose):
        print(" --->> Will process ", nevent, " events [ out of ", nb_evt, "] ")



    for ievent in range(nevent):
        ievent_tot += 1
        if(ntotevent > 0 and ievent_tot >= ntotevent):
            break

        if(verbose):
            print("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-")
            print(" READING EVENT ", ievent, " / ", ievent_tot)
            print("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-")

        dc.reset_event()

        if( read.read_event( data, ievent) < 0):
            print(" there is a pbm with the data, event ", ievent)
            continue

        ev = dc.event(date, ievent, ievent_tot, ts, tms)
        dc.evt_list.append(ev)

        if(verbose):
            dc.evt_list[-1].dump()

        ini_ped_mean = np.mean(dc.data, axis=2)
        ini_ped_rms = np.std(dc.data, axis=2)


        ped.ini_pedestal()
        ped.subtract_ped_mean()
        ped.store_raw_ped_rms()

    
    
        #plot.plot_wvf_multi_current([(0,8),(0,22),(0,45),(1,8),(1,22),(1,45)], option="raw", to_be_shown=True)
    
        #plot.plot_ed_data(option="rawpres", to_be_shown=False)
        #plot.plot_ed_data(option="rawzoompres", to_be_shown=False, adc_min=-10, adc_max=10)
        #plot.plot_event_ped_and_rms(option='initial', to_be_shown=True)
            
        """ FFT seems completely useless for the 50L data """
        ps = noise.FFTLowPass(lowpasscut, freqlines)            
        #plot.plot_fft_ps(ps, to_be_shown=False)



        """ 1st ROI attempt based on ADC cut + broken channels """
        """ useless in current 50L ana chain, as an estimate of ped RMS is already done """
        #noise.define_ROI_ADC(adc_thresh)

                        
        """ Make ROI based on ped rms """
        noise.define_ROI(signal_thresh, 2)    
        ped.subtract_ped_mean()


        #plot.plot_ed_data(option="ped_roi", to_be_shown=False)

        """Apply coherent filter(s) """
        noise.coherent_filter(coherent_groups)


        """ Update ROI regions """
        """ this part is not needed for 50L data """
        """
        noise.define_ROI(signal_thresh, 2)
        noise.median_filter(400)            
        """


        """ Update ROI regions """
        noise.define_ROI(signal_thresh_2, 2)

        ped.compute_and_subtract_ped()

        """ store final estimate of pedestal mean and rms """
        ped.store_final_ped_rms()
        #plot.plot_ed_data(option="filtpres", to_be_shown=False, adc_min=-10, adc_max=10)

        """ parameters : pad left (n ticks) pad right (n ticks), min dt, thr1 (collection), thr2(overlap collection), thr3 (induction) """
        hf.hit_finder(5, 10, 6, 3., 6., 2.)
        
        if(verbose):
            print(" found ", dc.evt_list[-1].nHits[0], " and ", dc.evt_list[-1].nHits[1], " hits")


        """ nothing found, store empty event case """
        if(np.sum(dc.evt_list[-1].nHits) == 0):

            filt_ped_mean = np.mean(dc.data, axis=2)
            filt_ped_rms = np.std(dc.data, axis=2)
            ps_fake += ps
            nfake += 1 

            store.store_pedestals_fake(output, ini_ped_mean, ini_ped_rms, filt_ped_mean, filt_ped_rms)
            store.store_event(output)
            continue


            

        """
        if(dc.evt_list[-1].nHits[0] > 0 and dc.evt_list[-1].nHits[0] < 4 and dc.evt_list[-1].nHits[1] > 0 and dc.evt_list[-1].nHits[1] < 4):
            for h in dc.hits_list:
                h.dump()
            plot.plot_ed_data_and_hits(option="ar39test", to_be_shown=True, adc_min=-10, adc_max=10)
        """

        #plot.plot_ed_data_and_hits(option="pres", to_be_shown=False, adc_min=-10, adc_max=10)
        #plot.plot_2dview_hits(option="test", to_be_shown=False)
    
    
        """ 1st search for most of the tracks"""
        """parameters : eps (cm), min pts, y axis squeeze"""
        #clus.dbscan(4, 3, 1.)

        """ parameters : search radius (cm), min nb of hits in cluster """
        #clus.mst(5, 3)
        
        [x.set_cluster(0) for x in dc.hits_list]
        dc.evt_list[-1].nClusters[0] = 1
        dc.evt_list[-1].nClusters[1] = 1
        
        """
        plot.plot_2dview_clusters(option="nofake", to_be_shown=False)
        """
        
        """parameters : min nb hits, rcut, chi2cut, y error, slope error, pbeta"""
        trk2d.find_tracks(5, 6., 8., 0.5, 1., 3., 1)

                

        """ parameters are : min distance in between 2 tracks end points, slope error tolerance, extrapolated distance tolerance + filter input parameters"""
        trk2d.stitch_tracks(10., 10., 4., 0.5, 1., 3., 1)


        """ no 2D tracks found, search for single hits then """
        """ stop the reco here anyway """
        if(np.sum(dc.evt_list[-1].nTracks2D) == 0 and dc.evt_list[-1].nHits[0] < 6 and dc.evt_list[-1].nHits[1] < 6):

            matched_hits = ar.argon_39_finder(3)
            
            if(len(matched_hits)>0):
                n_events_single_hits += 1
                for i, j in matched_hits:
                    store.store_single_hits(output, i, j)
            dc.evt_list[-1].nSingleHits = len(matched_hits)
            store.store_event(output)
            continue



        #plot.plot_2dview_hits_2dtracks(option="test", to_be_shown=False)
    

        """ parameters are : z start/end agreement cut (cm), v0/v1 charge balance, distance to detector boundaries for time correction (cm) """
        trk3d.find_tracks(3., 0.25, 2.)
        

        #plot.plot_2dview_hits_and_3dtracks(option="filter", to_be_shown=True)
        if(verbose):
            print(" Found ", dc.evt_list[-1].nTracks3D, " 3D Tracks ")


        store.store_event(output)
        if(dc.evt_list[-1].nTracks3D > 0):
            #plot.plot_ed_data_and_hits(option="3D", to_be_shown=True, adc_min=-10, adc_max=10)
            #plot.plot_3d("test", False)
            if(verbose):
                [x.dump() for x in dc.tracks3D_list]
                dc.evt_list[-1].dump_reco()

                        
            store.store_pedestals_data(output)
            store.store_hits(output)
            store.store_tracks2D(output)
            store.store_tracks3D(output)

            ps_muon += ps
            nevent_mu += 1

        """
        tstore = time.time()              
        print("time to store %.3f"%(time.time()-tstore))
        """

    data.close()
store.store_infos(output, day, run, ievent_tot, nevent_mu, n_events_single_hits, time.time())

ps_muon /= nevent_mu
store.store_fft_data(output, ps_muon)
ps_fake /= nfake
store.store_fft_fake(output, ps_fake)

output.close()
#output_sh.close()
#output_ped.close()
tottime = time.time() - tstart
print(" Has reconstructed ", nevent_mu, " events with at-least one 3D track out of ", ievent_tot, " events")
print(" TOTAL RUNNING TIME %.3f s == %.3f s/evt"% (tottime, tottime/ievent_tot))

