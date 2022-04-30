# a script to run for a folder of solar radio burst events

import sys

import matplotlib.dates as mdates
import datetime
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
import scipy

import matplotlib as mpl
# try to use the precise epoch
mpl.rcParams['date.epoch']='1970-01-01T00:00:00'
try:
    mdates.set_epoch('1970-01-01T00:00:00')
except:
    pass
import detectRadioburst as drb

from skimage.transform import probabilistic_hough_line

import detectRadioburst as drb


folder_of_interest = './L857852_SAP000_B000_S0_P000_bf/'
out_folder = './detect/'
plot_fig = True
dump_info_to_json = True
write_csv = True


import radioTools as rt
import glob
import json
import os

fnames = glob.glob(folder_of_interest+'/*.fits')

import os
csv_fname = 'event.csv'
os.system('rm '+csv_fname)
id_event  = 0
with open(csv_fname,'w') as fp:
    fp.write('''ID, t, t0_num, t1_num,f_0,f_1, dfdt(MHz/s), v_b(c)
             ''')
fp.close()

from tqdm import tqdm

for fname in tqdm(fnames):
    # read in 
    (dyspec,t_fits,f_fits,hdu)  = drb.read_fits(fname)
    (dyspec,f_fits) =  drb.cut_low(dyspec,f_fits,f_low_cut_val=25)
    (data_fits_new_tmp,data_fits_new) = drb.preproc(
        dyspec,gauss_sigma=1.5)
    bmap = drb.binarization(data_fits_new,N_order=6,peak_r=1)
    lines = drb.hough_detect(bmap,dyspec,threshold=50,
                             line_gap=10,line_length=25,
            theta=np.linspace(np.pi/2-np.pi/8,np.pi/2-1/180*np.pi,300))
    if len(lines)<=1:
        continue
    line_sets = drb.line_grouping(lines)
    (v_beam, f_range_burst, t_range_burst, model_curve_set,
         t_set_arr_set,f_set_arr_set,t_model_arr,f_model_arr
        )= drb.get_info_from_linegroup(line_sets,t_fits,f_fits)

    fname_json  = fname.replace('.fits','.json')
    with open(fname_json, 'r') as fp:
        old_json = json.load(fp)
    fp.close()
    drb.append_into_json(old_json, v_beam, f_range_burst, t_range_burst)
    with open(out_folder+os.path.basename(fname_json), 'w') as fp:
        json.dump(old_json,fp)
        fp.close()
        
    if plot_fig==True:
        fig,ax = plt.subplots(1,1,figsize=[6,3],dpi=200)
        lines = sorted(lines, key=lambda i: i[0][1])
        ax.imshow(data_fits_new.T,aspect='auto',origin='lower', 
                           vmin=(np.mean(data_fits_new)-2*np.std(data_fits_new)),
                           vmax=(np.mean(data_fits_new)+3*np.std(data_fits_new)),cmap='gray',
                           extent=[t_fits[0],t_fits[-1],f_fits[0],f_fits[-1]])
        for idx,model in enumerate(model_curve_set):
            plt.plot(model[0],model[1],ls='--')
            plt.plot(t_range_burst[idx],f_range_burst[idx],'k+')



        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.set_xlabel('Time (UT)')
        ax.set_ylabel('Frequency (MHz)')
        ax.set_title(hdu[0].header['CONTENT'])
        
        fig.savefig(out_folder+os.path.basename(fname)+'.jpg')
        
    if write_csv==True:
        with open(csv_fname,'a') as fp:

            for idx,v_cur in enumerate(v_beam):
                fp.write(str(id_event)+','+mdates.num2date(t_range_burst[idx][0]).strftime("%H:%M:%S")+','
                     +str(t_range_burst[idx][0])+','+str(t_range_burst[idx][1])+','
                     +str(f_range_burst[idx][0])+','+str(f_range_burst[idx][1])+','
                     +str((np.max(f_set_arr_set[idx])-np.min(f_set_arr_set[idx]))/
                     (np.max(t_set_arr_set[idx])-np.min(t_set_arr_set[idx])))+','
                     +str(v_beam[idx])
                     +'''
                     ''')
                id_event+=1
        fp.close()
        
        