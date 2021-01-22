import config as cf
import data_containers as dc
import pedestals as ped

import numpy as np
import numexpr as ne 

import bottleneck as bn

from sklearn import linear_model



def define_ROI_ADC(thresh):
    thresh_coll = thresh[0]
    dc.mask[0,:,:] = ne.evaluate( "where((data > thresh_coll) | ~alive_chan, 0, 1)", global_dict={'data':dc.data[0,:,:], 'alive_chan':dc.alive_chan[0,:,:]}).astype(bool)
    thresh_ind_p = thresh[1]
    dc.mask[1,:,:] = ne.evaluate( "where((data > thresh_ind_p) | ~alive_chan, 0, 1)", global_dict={'data':dc.data[1,:,:], 'alive_chan':dc.alive_chan[1,:,:]}).astype(bool)
    thresh_ind_m = thresh[2]
    dc.mask[1,:,:] = ne.evaluate( "where((data < thresh_ind_m) | ~alive_chan |~mask, 0, 1)", global_dict={'data':dc.data[1,:,:], 'alive_chan':dc.alive_chan[1,:,:], 'mask':dc.mask[1,:,:]}).astype(bool)


def define_ROI(thresh, iteration):
    #ne.set_num_threads(4) #does not speed up things
    thresh_coll = thresh[0]
    thresh_ind_p = thresh[1]
    thresh_ind_m = thresh[2]

    """ Update the mask based on pedestal RMS """    
    for it in range(iteration):

        dc.ped_rms = dc.ped_rms[:,:,None]
        dc.ped_mean = dc.ped_mean[:,:,None]


        dc.mask[0,:,:] = ne.evaluate( "where((data > mean + thresh_coll*rms) | (~alive_chan), 0, 1)", global_dict={'data':dc.data[0,:,:], 'alive_chan':dc.alive_chan[0,:,:], 'rms':dc.ped_rms[0,:,:], 'mean':dc.ped_mean[0,:,:]}).astype(bool)

        dc.mask[1,:,:] = ne.evaluate( "where((data > mean + thresh_ind_p*rms) | (~alive_chan), 0, 1)", global_dict={'data':dc.data[1,:,:], 'alive_chan':dc.alive_chan[1,:,:], 'rms':dc.ped_rms[1,:,:], 'mean':dc.ped_mean[1,:,:]}).astype(bool)

        dc.mask[1,:,:] = ne.evaluate( "where((data < mean - thresh_ind_m*rms) | (~alive_chan) | ~mask, 0, 1)", global_dict={'data':dc.data[1,:,:], 'alive_chan':dc.alive_chan[1,:,:], 'mask':dc.mask[1,:,:], 'rms':dc.ped_rms[1,:,:], 'mean':dc.ped_mean[1,:,:]}).astype(bool)
        
        dc.ped_rms = np.squeeze(dc.ped_rms, axis=2)
        dc.ped_mean = np.squeeze(dc.ped_mean, axis=2)

        ped.compute_pedestal_RMS()
        ped.compute_pedestal_mean()

def coherent_filter(groupings):
    """
    1. Computes the mean along group of channels for non ROI points
    2. Subtract mean to all points
    """

    for group in groupings:
        if( (cf.n_ChanPerView % group) > 0):
            print(" Coherent Noise Filter in groups of ", group, " is not a possible ! ")
            return

        nslices = int(cf.n_ChanPerView / group)
        
        dc.data = np.reshape(dc.data, (cf.n_View, nslices, group, cf.n_Sample))
        dc.mask = np.reshape(dc.mask, (cf.n_View, nslices, group, cf.n_Sample))


    
        """sum data if mask is true"""
        with np.errstate(divide='ignore', invalid='ignore'):
            """sum the data along the N channels (subscript k) if mask is true,
            divide by nb of trues"""
            mean = np.einsum('ijkl,ijkl->ijl', dc.data, dc.mask)/dc.mask.sum(axis=2)

            """require at least 3 points to take into account the mean"""
            mean[dc.mask.sum(axis=2) < 3] = 0.
        

        """Apply the correction to all data points"""
        dc.data -= mean[:,:,None,:]
        
        """ restore original data shape """
        dc.data = np.reshape(dc.data, (cf.n_View, cf.n_ChanPerView, cf.n_Sample))
        dc.mask = np.reshape(dc.mask, (cf.n_View, cf.n_ChanPerView, cf.n_Sample))



def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def FFTLowPass(lowpass, freqlines) :

    """it's 324 (n/2+1) points from 0Hz to 1./(2*sampling) = 1.25MHz (nyquist freq)"""
    
    n    = int(cf.n_Sample/2) + 1
    rate = 1./cf.n_Sampling #in MHz
    freq = np.linspace(0, rate/2., n)

    """define gaussian low pass filter"""
    gauss_cut = np.where(freq < lowpass, 1., gaussian(freq, lowpass, 0.02))

    """go to frequency domain"""
    fdata = np.fft.rfft(dc.data)


    regr = linear_model.LinearRegression()
    regi = linear_model.LinearRegression()

    """remove specific frequencies"""

    # smooth the frequencies removed from prev and next freq. value (linear fit)
    # still introduce artefacts - not recommended to use

    for f in freqlines:
        fbin = int(f * cf.n_Sample * cf.n_Sampling)
        for iview in range(cf.n_View):
            for ichan in range(cf.n_ChanPerView):
                ptsx = []
                ptsyr = []
                ptsyi = []

                for i in reversed(range(4)):
                    ptsx.append(fbin - 5 - i)
                    ptsyr.append(fdata[iview, ichan, fbin - 5 - i].real)
                    ptsyi.append(fdata[iview, ichan, fbin - 5 - i].imag)
                for i in range(4):
                    ptsx.append(fbin + 5 + i)
                    ptsyr.append(fdata[iview, ichan, fbin + 5 + i].real)
                    ptsyi.append(fdata[iview, ichan, fbin + 5 + i].imag)
                ptsx = np.asarray(ptsx).reshape(-1,1)
                ptsyr = np.asarray(ptsyr)
                ptsyi = np.asarray(ptsyi)

                regr.fit(ptsx, ptsyr)
                regi.fit(ptsx, ptsyi)


                xtorm = np.asarray(range(fbin-4,fbin+5)).reshape(-1,1)
                yvalr = regr.predict(xtorm)
                yvali = regi.predict(xtorm)

                for ib in range(9):
                    fdata[iview,ichan,fbin-4+ib] = complex(yvalr[ib], yvali[ib]) 
          
        #gauss_cut[max(fbin-2,0):min(fbin+3,n)] = 0.2
        #gauss_cut[max(fbin-1,0):min(fbin+2,n)] = 0.1


    """get power spectrum (before cut)"""
    #ps = 10.*np.log10(np.abs(fdata)+1e-1) 

    
    """Apply filter"""
    fdata *= gauss_cut[None, None, :]

    """go back to time"""
    dc.data = np.fft.irfft(fdata)


    """get power spectrum after cut"""
    ps = 10.*np.log10(np.abs(fdata)+1e-1) 
    return ps

    
def centered_median_filter(array, size):
    """ pads the array such that the output is the centered sliding median"""
    rsize = size - size // 2 - 1
    array = np.pad(array, pad_width=((0, 0), (0, 0) , (0, rsize)),
                   mode='constant', constant_values=np.nan)
    return bn.move_median(array, size, min_count=1, axis=-1)[:, :, rsize:]


def median_filter(window):
    """ apply median filter on data to remove microphonic noise """

    """ mask the data with nan where potential signal is (ROI)"""
    med = centered_median_filter(np.where(dc.mask[:,:,:]==True, dc.data[:,:,:], np.nan), window)

    """ in median computation, if everything is masked (all nan) nan is returnedso changed these cases to 0 """
    med = np.nan_to_num(med, nan=0.)

    """ apply correction to data points """
    dc.data[:,:,:] -= med
