data_path = "/eos/user/l/lzambell/analysis/lardon_vert/50L_spring20/data"
store_path = "/eos/user/l/lzambell/analysis/lardon_vert/50L_spring20/reco"

file_prename = "WIB00step18_FEMB_B8_"

n_View = 4
n_Sample = 646
n_ChanPerView = 64
n_ChanTot = n_View * n_ChanPerView
n_Sampling = 0.5 #in mu-seconds
ChanPitch = 0.5 #cm
LAr_Temperature = 87.5 #(?)
E_drift = 0.5 #kV/cm
Anode_Z = 52. #cm
len_det_x = 32. #cm
len_det_y = 32. #cm

view_type = {0:'Collection',1:'Induction', 2: 'Induction', 3: 'Shield'}


evt_byte_size = 129528 # TODO update
n_femb = 8 #Card number?

ADCperfC = 33.7078 #from paper preprint (5.4e-3 ADCtick/e)
AreaCorr = 1.269*4.
#ADCperfC = 48.016 #from paper preprint (5.4e-3 ADCtick/e)

""" BROKEN CHANNELS TO BE REMOVED FROM THE ANALYSIS"""

""" provided in DAQ Channel """
daq_broken_channels = [] #3508, 3507, 3505, 3504]

crp_broken_channels = []


""" or provided in (view, channel) tuple """
