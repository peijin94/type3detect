# Type III Radio Burst Automatic Recognition Algorithm

The automatic recognition of the Type III radio burst using Hough Transform.

## Installation

### Python

install from pypi
```bash
pip install type3detect
```

install from git
```bash
git clone https://github.com/peijin94/type3detect.git
cd src/python
python pip install .
```

### Matlab:
```bash
git clone https://github.com/peijin94/type3detect.git
cd src/matlab
```

## Example 

Read fits file and detect type III radio bursts in dynamic spectrum.

```python
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from type3detect import detectRadioburst as drb
from type3detect import radioTools as rt

fname  = './LOFAR_20220701_070000_LBA_OUTER_S0.fits'

# preprocess
(dyspec,t_fits,f_fits,hdu) = drb.read_fits(fname) # read LOFAR dynspec fits
(dyspec,f_fits) = drb.cut_low(dyspec,f_fits,f_low_cut_val=30) # remove below freq (RFI)
(data_fits_new,data_fits_new_smooth) = drb.preproc(dyspec,gauss_sigma=1.5)

# binarization
bmap = drb.binarization(data_fits_new,N_order=6,peak_r=1.002)

# detect verticle features
lines = drb.hough_detect(bmap,dyspec,threshold=40,line_gap=10,line_length=30,
            theta=np.linspace(np.pi/2-np.pi/8,np.pi/2-1/180*np.pi,300))
line_sets = drb.line_grouping(lines)

# get electron beam information from radio bursts
(v_beam, f_range_burst, t_range_burst, model_curve_set,
     t_set_arr_set,f_set_arr_set,t_model_arr,f_model_arr
    )= drb.get_info_from_linegroup(line_sets,t_fits,f_fits)

# detailed demo in ./type3detect_demo.ipynb
```

Example: Implementation with for LOFAR dynamic spectrum.

Demo Notebook: [type3detect_demo.ipynb](https://github.com/peijin94/type3detect/blob/master/type3detect_demo.ipynb)

![img](https://github.com/peijin94/type3detect/raw/master/img/LOFAR_20220413_135000_LBA_OUTER.fits.jpg)


### Binarization

Use several ways to transform the flux intensity data into binary. ([LocalMax(matlab)](src/matlab/get_local_max_map.m)) ([binarization(Python)](src/python/src/detectRadioBurst.py))


### Active-Contour method for Backbone

The Active Contour Method for the backbone of the radio burst (**ACBone**)

Take the result of the Hough transform as a initail position and iteratively move the line to find the backbone position. The demo code : ([ACBone-Matlab](src/matlab/active_contour.m), [ACBone-Python](src/python/type3detect/ACBone.py))

![img](https://github.com/peijin94/type3detect/raw/master/img/activecontour.GIF)

Eventually, we can obtain the centerline of a Type III radio burst.


## Citation

Make sure to cite the paper if you use the idea or code in this repo: [A type III radio burst automatic analysis system and statistic results for a half solar cycle with Nançay Decameter Array data](https://www.aanda.org/component/article?access=doi&doi=10.1051/0004-6361/201833260#R16) Peijin Zhang. A&A 2018.10

bibtex:
```
@article{zhang2018type,
  title={A type III radio burst automatic analysis system and statistic results for a half solar cycle with Nan{\c{c}}ay Decameter Array data},
  author={Zhang, PJ and Wang, Chuan Bing and Ye, Lin},
  journal={Astronomy \& Astrophysics},
  volume={618},
  pages={A165},
  year={2018},
  publisher={EDP Sciences}
}
```
