'''
    File name: detectRadioburst.py
    Author: Peijin Zhang 张沛锦
    Date : 2022-4-30
    
    # some useful functions to detect solar radio bursts
    # 自动识别射电暴程序

    Notes:

    1. input dspec should be in the shape of (time, frequency)
'''

import sys
sys.path.insert(1, '../') # make sure to use the code in this repo

import matplotlib.dates as mdates
import datetime
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
import scipy
from scipy import interpolate,optimize

from scipy.optimize import curve_fit

import matplotlib as mpl
# try to use the precise epoch
mpl.rcParams['date.epoch']='1970-01-01T00:00:00'
try:
    mdates.set_epoch('1970-01-01T00:00:00')
except:
    pass
from . import radioTools as rt
from skimage.transform import probabilistic_hough_line

# frequency [MHz]


def read_fits_lofar(fname):
    # for lofar fits file
    hdu = fits.open(fname)
    dyspec = np.array(hdu[0].data)
    f_fits = np.array(hdu[1].data['FREQ'][0])
    t_fits = np.array(hdu[1].data['TIME'][0])
    return (dyspec,t_fits,f_fits,hdu)


def idx_val_pos(f_fits, target):
    return np.abs(np.array(f_fits)-target).argmin()

def cut_low(dyspec,f_fits,f_low_cut_val=21):
    idx_cut = idx_val_pos(f_fits, f_low_cut_val)
    dyspec=dyspec[:,idx_cut:]
    f_fits = f_fits[idx_cut:]
    return (dyspec,f_fits)
    
    
def preproc(dyspec,gauss_sigma=1.5, background_normalize=True):
    # const background removal and gaussian smooth
    if background_normalize:
        data_fits_new_tmp = (dyspec / np.nanmean(
                np.sort(dyspec, 0)[
                int(dyspec.shape[0] * 0.1):int(dyspec.shape[0] * 0.3), :], 0))-1
    else:
        data_fits_new_tmp = dyspec
    data_fits_new  = scipy.ndimage.gaussian_filter(data_fits_new_tmp,
        gauss_sigma, order=0, output=None,  cval=0.0, truncate=5.0,mode='nearest')
    return (data_fits_new_tmp,data_fits_new)
    
def binarization(data_fits_new,N_order=6,peak_r=0.99):
    # with high order local-max method 
    bmap = np.ones_like(data_fits_new)
    N_pad = N_order
    local_max_arr = np.pad(data_fits_new,((N_pad,N_pad),(0,0)))
    for idx in range(N_pad-1):
        bmap=bmap* ((peak_r*local_max_arr[N_pad+idx+1:-N_pad+idx+1,:]<
                            local_max_arr[N_pad+idx:-N_pad+idx,:]) & 
                    (peak_r*local_max_arr[N_pad-idx-1:-N_pad-idx-1,:]<
                            local_max_arr[N_pad-idx:-N_pad-idx,:]) )
    return bmap

def hough_detect(bmap,dyspec,threshold=50,line_gap=10,line_length=25,
            theta=np.linspace(np.pi/2-np.pi/8,np.pi/2-1/180*np.pi,300)):
    lines = probabilistic_hough_line(bmap, threshold=threshold,line_gap=line_gap,line_length=line_length,
                                 theta=theta)
    return lines

    
norm = np.linalg.norm
def point_to_line_distance(p1,p2,p3):
    d = np.abs(norm(np.cross(p2-p1, p1-p3)))/norm(p2-p1)
    return d

def point_to_point_distance(p1,p2):
    return norm(p1-p2)

def line_grouping(lines,min_dist=3, threshmode=2): # pix
    """
    Parameters
    ----------
    lines : list
        list of lines detected by hough transform.
    min_dist : float, optional  # pix
        minimum distance between two lines to be grouped together.
    threshmode : int, optional
        DESCRIPTION. 1: only distance of point to line
                        2: distance of point to line and angle between lines
    """
    # group the detected lines into group in regard of events
    lines = sorted(lines, key=lambda i: (i[0][1]+i[1][1])/2)
    line_sets = [[lines[0]]]
    for idx,line in enumerate(lines[0:-1]):
        (A,B),(C,D) = np.array([lines[idx], lines[idx+1] ])

        # use the longer line as the reference
        l_AB = norm(A-B)
        l_CD = norm(C-D)
        if l_AB<l_CD:
            A,B,C,D = C,D,A,B
        
        # point to line distance
        Line_dist_thresh = np.min([point_to_line_distance(A,B,C),point_to_line_distance(A,B,D)])< min_dist
        Line_dist_thresh2 = np.min([point_to_line_distance(A,B,C),point_to_line_distance(A,B,D)])< min_dist*1.5

        # two lines segment not too far away
        Point_dist_thresh = np.min([point_to_point_distance(A,C),point_to_point_distance(A,D),
                                    point_to_point_distance(B,C),point_to_point_distance(B,D)])< np.max(
                                        [point_to_point_distance(A,B),point_to_point_distance(C,D)])
       
        if threshmode==1:
            final_thresh = Line_dist_thresh
        elif threshmode==2:
            final_thresh = Line_dist_thresh & Point_dist_thresh
        elif threshmode==3:
            if B[0]>C[0] and B[1]>C[1]:
                final_thresh = Line_dist_thresh2 & Point_dist_thresh

        if final_thresh:
            # the line join
            line_sets[len(line_sets)-1].append(lines[idx+1])
        else:
            # new set
            line_sets.append([lines[idx+1]])
    
    return line_sets



def get_info_from_linegroup(line_sets,t_fits,f_fits):

    # mapping from t and f to index of x and y
    t_idx_arr = np.arange(0, t_fits.shape[0])
    f_idx_arr = np.arange(0, f_fits.shape[0])
    t_interf = interpolate.interp1d(t_idx_arr, t_fits)
    f_interf = interpolate.interp1d(f_idx_arr, f_fits)
    
    v_beam = []
    f_range_burst = []
    t_range_burst = []
    model_curve_set = []
    t_set_arr_set = []
    f_set_arr_set = []
    t_set_arr = []
    f_set_arr = []
    t_model_arr = []
    f_model_arr = []

    for lines in line_sets:
        if len(lines)==1:
            continue
        try:
            x_set=[]
            y_set=[]
            for line in lines:
                x_set.append(line[0][1])
                x_set.append(line[1][1])
                y_set.append(line[0][0])
                y_set.append(line[1][0])

            t_set_arr = (t_interf(x_set) - np.min(t_fits))*24*3600
            f_set_arr = f_interf(y_set)

            popt, pcov = optimize.curve_fit(rt.freq_drift_f_t,
                t_set_arr, f_set_arr, p0=(0.1,np.min(t_set_arr)-3./3600/24), method="lm")

            t_model_arr  = np.linspace(rt.freq_drift_t_f(np.min(f_set_arr),*popt) ,
                                       rt.freq_drift_t_f(np.max(f_set_arr),*popt),50) 
            f_model_arr = rt.freq_drift_f_t(t_model_arr,popt[0],popt[1])

            t_model_arr = t_model_arr/(24*3600)+np.min(t_fits)


            model_curve_set.append([t_model_arr,f_model_arr])
            t_range_burst.append( [rt.freq_drift_t_f(np.min(f_set_arr),*popt)[0]/(24*3600)+np.min(t_fits) ,
                                   rt.freq_drift_t_f(np.max(f_set_arr),*popt)[0]/(24*3600)+np.min(t_fits) ] )
            f_range_burst.append([np.min(f_set_arr),np.max(f_set_arr)])
            v_beam.append(popt[0])
            t_set_arr_set.append(t_set_arr)
            f_set_arr_set.append(f_set_arr)
        except:
            pass

    
    return (v_beam, f_range_burst, t_range_burst, model_curve_set,
            t_set_arr_set,f_set_arr_set,
           t_model_arr,f_model_arr)

def append_into_json(old_json, v_beam, f_range_burst, t_range_burst):
    
    event_detail = []
    for idx,v_cur in enumerate(v_beam):
        event_detail.append({
            'v_beam':v_cur,
            'freq_range':((f_range_burst[idx])),
            'time_range':((t_range_burst[idx])),
            'str_time':mdates.num2date(t_range_burst[idx][0]).strftime("%H:%M:%S")})
    
    old_json['event']={
        'detection': True,
        'type':'III',
        'detail': event_detail}
    
    return old_json


def averaging(arr_query, n_point, axis=0, start_idx=-1, end_idx=-1):
    start_idx = 0 if start_idx < 0 else start_idx
    end_idx = arr_query.shape[axis] if end_idx < 0 else end_idx
    out_size = ((end_idx-start_idx) // n_point)

    res = 0
    if axis == 1:
        res = np.mean(np.array(([arr_query[:, (start_idx+idx):(start_idx+(out_size)*n_point+idx):n_point]
                                 for idx in range(n_point)])), axis=0)
    else:
        res = np.mean(np.array(([arr_query[(start_idx+idx):(start_idx+(out_size)*n_point+idx):n_point, :]
                                 for idx in range(n_point)])), axis=0)
    return res

def model_based_shift_per_ch(data, t_ax, f_ax, vbeam=0.15, model=rt.parkerfit, ops_forward = True):
    """ 
    shift the data based on the model curve of type III burst.
    given a speed of the electron
    
    Parameters:
    -----------
    data: np.array
        the data to be shifted, shape: (time, frequency)
    t_ax: np.array
        time axis of the data, in mjd
    f_ax: np.array
        frequency axis of the data in MHz
    vbeam: float
        the speed of the electron beam, in c
    model: function
        the model function of the type III burst, default is parkerfit
    ops_forward: bool
        if True, shift the data forward, else shift the data backward
        (forward means from raw spec to shifted spec)
        (backward means from shifted spec to raw spec)
    
    Returns:
    --------
    data_aligned: np.array
        the aligned data, shape: (time, frequency)
    """

    # assuming 0.2c of beam speed for the channel-wised time shifting alignment
    delay_ch = np.array([rt.freq_drift_t_f(f_this, vbeam, 0) for f_this in f_ax])
    dt_raw =  np.mean(np.diff(t_ax))*24*3600 
    idx_shift = np.round(delay_ch /dt_raw ).astype(int)
    # align the data
    data_aligned = np.zeros_like(data)

    for i in range(data.shape[1]):
        if ops_forward:
            data_aligned[:,i] = np.roll(data[:,i], -idx_shift[i])
        else:
            data_aligned[:,i] = np.roll(data[:,i], idx_shift[i])

    return data_aligned


def pix_to_val(pix,arr):
    """ 
    Input a number in pixel space, output the value in the original space
    arr should be monotonically increasing or decreasing
    """

    if pix<0:
        return arr[0]
    elif pix>(arr.shape[0])-1:
        return arr[-1]
    else:
        boundary = [np.floor(pix),np.floor(pix)+1]
        return (pix-boundary[0])*arr[int(boundary[1])]+(boundary[1]-pix)*arr[int(boundary[0])]
    
def val_to_pix(val,arr):
    """ 
    Input a number in original space, output the value in the pixel space
    arr should be monotonically increasing or decreasing
    """
    if val<arr[0]:
        return 0
    elif val>arr[-1]:
        return (arr.shape[0])-1
    else:
        pix_lower = np.abs(arr-val).argmin()
        pix_upper = pix_lower+1
        return pix_lower+(val-arr[pix_lower])/(arr[pix_upper]-arr[pix_lower])
