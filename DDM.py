'''This module contains a collection of function to perform differential dynamic microscopy analysis of videos'''

from ipywidgets import interactive
import matplotlib.pyplot as plt
from scipy import fftpack
import pandas as pd
import numpy as np
import ipywidgets

def browse_images_FFT(video,interval=1,muperpix=1):
    '''
    Expect pims video object
    
    construct widget to explore DDM signal
    '''
    frames=len(video)
    
    try:
        video=video[:,:,1]
    except:
        pass
        
        
    def view_image(framenum,delta):
        plt.figure(figsize=(12,8))
        
        maxpix=min(video.frame_shape[:2])

        
        im1=video[framenum][:maxpix,:maxpix].astype(np.float)
        im2=video[framenum+delta][:maxpix,:maxpix].astype(np.float)
        imdiff=im2-im1

        F1=fftpack.fft2(imdiff)
        F2 = fftpack.fftshift( F1 )
        psd2D = np.abs( F2 )**2
        psd1D = azimuthalAverage(psd2D)
        
        
        plt.subplot(2,3,1)
        plt.title('(a) image 1 $I(r,t)$, time: ' + str(framenum*interval) + ' s')
        plt.imshow( im1)
        plt.set_cmap('gray')
        plt.axis('off')
        
        plt.subplot(2,3,2)
        plt.title('(b) image 2 $I(r,t+ \Delta t)$, time: ' + str((framenum+delta)*interval) + ' s')
        plt.imshow( im2)
        plt.set_cmap('gray')
        plt.axis('off')
        
        plt.subplot(2,3,3)
        plt.title('(c) image2-image1, delta t:'+ str((delta)*interval) +
              ' s\n $I( q ,t+ \Delta t ) - I(q,t)$')
        plt.imshow( imdiff)
        plt.set_cmap('gray')
        plt.axis('off')
        
        plt.subplot(2,3,4)
        plt.title('(d)FFT - $g(q,t)$ ')
        plt.imshow( np.log10( psd2D ))
        plt.subplot(2,3,5)
        plt.plot( psd1D )
        plt.yscale('log')
    
        plt.xlabel('Spatial Frequency')
        plt.ylabel('Power Spectrum')
        plt.subplot(2,3,6)
        plt.plot( psd1D )
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Spatial Frequency')
        plt.ylabel('Power Spectrum')
        plt.show()
        
    w=interactive(view_image, framenum=(0,frames),delta=(1,frames-100))
    return w

def azimuthalAverage(image, center=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr

    return radial_prof

def _calculate_iq_delta(video,delta,naverage=None):
    '''
    Expect pims video object and delta frame interval to consider
    
    return mean azimutal average over naverage diff frame with dinstance delta
    '''
    
    frames=len(video)
    progress_dt.value=str(delta)
    min(video.frame_shape[:2])    
    iq=[]
    
    try:
        video=video[:,:,1]
    except:
        pass
    
    if naverage is None:
        naverage=frames-delta
    
    framelist=range(0,frames-delta,int((frames-delta)/naverage))
    
    for frame in framelist:
        
        maxpix=min(video.frame_shape[:2])

        
        im1=video[frame][:maxpix,:maxpix].astype(np.float)
        im2=video[frame+delta][:maxpix,:maxpix].astype(np.float)
        imdiff=im2-im1

        F1=fftpack.fft2(imdiff)
        F2 = fftpack.fftshift( F1 )
        psd2D = np.abs( F2 )**2
        psd1D = azimuthalAverage(psd2D)
        
        iq.append(psd1D)
        
        progress_frame_pair.value=str(frame)
        
    return np.mean(iq,0)

def _calculate_iq_delta_multi(video,deltalist,naverage=None):
    '''
    Expect pims video object and delta frame interval list to consider
    
    return mean azimutal average over naverage diff frame with dinstance delta for each delta frame in the provided list
    '''

    frames=len(video)
    iq_matrix=[]
        
    for delta in deltalist:
        iq_matrix.append(_calculate_iq_delta(video,delta,naverage))
    
    return np.array(iq_matrix)

def calculate_DDM(video,naverage=None,numdt=None, interval=1,muperpix=1):
    '''
    Expect pims video object
    naverage : number of frame pair to average over for each delta frame
    numdt : number of delta frame value to consider (will produce a log space list)
    interval : time between frames in sec
    muperpix: physical lengh of each pixel in micron
    '''
    global progress_frame_pair
    
    progress_frame_pair=ipywidgets.Text(description='current differential frame realization', value='')
    
    global progress_dt
    
    progress_dt=ipywidgets.Text(description='current dt', value='0')
    
    display(progress_frame_pair)
    display(progress_dt)
    
    
    if numdt is None:
        numdt=20
    
    if naverage is None:
        naverage=10
    
    frames=len(video)
    
    frame_len=min(video.frame_shape[:2])

    
    
    deltalist=np.unique(np.logspace(0, np.log10(frames-int(frames/2)), num=numdt, endpoint=True, base=10.0, dtype=int))
    
    deltalistsec=np.transpose(np.array(deltalist)*interval)

    iq=_calculate_iq_delta_multi(video,deltalist,naverage)
    
    xax=np.array(range(iq.shape[1]))*2*3.14/(muperpix*frame_len)
    
    result=pd.DataFrame(iq,index=deltalistsec, columns=xax)
    
    video.result_FFT=result    
    return video

from lmfit import  Model
from IPython.display import clear_output
import lmfit

def DDM_func(x, off, amp, tau):
    '''
    single exponential decay function expected for bronian motion of non interacting particles
    '''
    return (off + amp * (1 - np.exp(-x/tau)))

DDMmod = Model(DDM_func)


def calculate_viscosity(video,radius=1,muperpix=1,qmin=None,qmax=None):
    '''
    expects pims video object
    radius: of the diffusing objects in um
    muperpix : calibration
    qmin: index of min q to consider in the fit
    qmax: index of max q to consider in the fit
    assumes room temperature
    '''
    if qmin is None:
        qmin=1
        
        
    frames=len(video)
    frame_len=min(video.frame_shape[:2])
    
    result=video.result_FFT
    xax=result.columns
    result=result[np.sort(list(xax))]
    deltalistsec=np.array(result.index)
    iq2=result.values
    xax=result.columns
    
    if qmax is None:
        qmax=len(xax)
    
    tau=[]
    xax=[]
    
    for i in range(qmin,qmax):
        
        try:
            result = DDMmod.fit(iq2[1:,i], x=deltalistsec[1:], amp=np.max(iq2[1:,i])-np.min(iq2[1:,i]),
                                off=min(iq2[1:,i]),
                                tau=np.mean(deltalistsec[1:]))
        
            tau.append(result.params['tau'].value)
            xax.append(i*2*3.14/(muperpix*frame_len))
        except:
            pass
            
    PL_model=lmfit.models.PowerLawModel()
    pars = PL_model.guess(np.array(tau), x=xax)
    out  = PL_model.fit(np.array(tau), params=pars, x=xax)
    amplitude=out.params['amplitude'].value
    
    plt.plot(xax,np.transpose(np.array(tau)),'o')
    plt.plot(xax,out.best_fit)
    #xscale('log')
    #yscale('log')
    plt.xlabel('q [um-1]')
    plt.ylabel('tau [s]')
    
    Dm=1/amplitude
    viscosity=1.3806488E-23 * (273+25) / (6*3.14*Dm*radius*1e-18)
    print(viscosity)
    
    return viscosity

def calculate_radius(video,viscosity=0.001,muperpix=1,qmin=None,qmax=None):
    '''
    expects pims video object
    radius: of the diffusing objects in um
    muperpix: calibration
    qmin: index of min q to consider in the fit
    qmax: index of max q to consider in the fit
    assumes room temperature
    '''
    if qmin is None:
        qmin=1
        
        
    frames=len(video)
    frame_len=min(video.frame_shape[:2])
    
    result=video.result_FFT
    xax=result.columns
    result=result[np.sort(list(xax))]
    deltalistsec=np.array(result.index)
    iq2=result.values
    xax=result.columns
    
    if qmax is None:
        qmax=len(xax)
    
    tau=[]
    xax=[]
    
    for i in range(qmin,qmax):
        
        try:
            result = DDMmod.fit(iq2[1:,i], x=deltalistsec[1:], amp=np.max(iq2[1:,i])-np.min(iq2[1:,i]),
                                off=min(iq2[1:,i]),
                                tau=np.mean(deltalistsec[1:]))
        
            tau.append(result.params['tau'].value)
            xax.append(i*2*3.14/(muperpix*frame_len))
        except:
            pass
            
    PL_model=lmfit.models.PowerLawModel()
    pars = PL_model.guess(np.array(tau), x=xax)
    out  = PL_model.fit(np.array(tau), params=pars, x=xax)
    amplitude=out.params['amplitude'].value
    
    plt.plot(xax,np.transpose(np.array(tau)),'o')
    plt.plot(xax,out.best_fit)
    #xscale('log')
    #yscale('log')
    plt.xlabel('q [um-1]')
    plt.ylabel('tau [s]')
    
    Dm=1/amplitude
    radius=1.3806488E-23 * (273+25) / (6*3.14*Dm*viscosity*1e-18)
    print(radius)
    
    return radius


def explore_iq_dt(video, interval=1):
    '''
    Expect pims video object
    interval: time between frames in sec
    construct widget to explore DDM results
    '''
    frames=len(video)
    frame_len=min(video.frame_shape[:2])
    
    result=video.result_FFT
    xax=result.columns
    result=result[np.sort(list(xax))]
    deltalistsec=result.index
    iq2=result.values
    xax=result.columns
    
    def view_plot(selectq,deltat):
        plt.figure(figsize=(15,3))
        
        gamma=[]
        aq=[]
        bq=[]
        
        plt.subplot(1,3,1)
        plt.plot(xax,iq2[deltat,:],'o')
        plt.yscale('log')
        plt.xscale('log')
        plt.ylabel('I(q) [a.u.]')
        plt.xlabel('q [um-1]')
        plt.axvline(xax[selectq], color='k', linestyle='--')
        plt.title(r'I(q) fixed $\Delta t$')
        
        
        plt.subplot(1,3,2)
        plt.plot(deltalistsec,iq2[:,selectq],'o')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('g(q,dt) [a.u]')
        plt.xlabel('Delta t [s]')


                 
        result = DDMmod.fit(iq2[1:,selectq], x=deltalistsec[1:], amp=max(iq2[1:,selectq])-min(iq2[1:,selectq]),
                            off=min(iq2[1:,selectq]),
                            tau=np.mean(deltalistsec[1:]));



        #clear_output()
        tau=result.params['tau'].value
        off=result.params['off'].value
        amp=result.params['amp'].value

        plt.plot(deltalistsec[1:], result.best_fit, 'r-')
        plt.axvline(deltat*interval, color='k', linestyle='--')

        plt.subplot(1,3,3)
        plt.plot(deltalistsec,-((iq2[:,selectq]-(off+amp))/amp),'o')
        plt.plot(deltalistsec,np.exp(-deltalistsec/tau))
        plt.axvline(deltat*interval, color='k', linestyle='--')
        plt.xscale('log')
        plt.ylabel('f(q,dt) [a.u.]')
        plt.xlabel('Delta t')
      
        
        
    w=interactive(view_plot, selectq=(1,len(iq2[1,:]-1)),deltat=(1,len(deltalistsec)-1))
    return w