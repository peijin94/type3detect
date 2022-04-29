# Type III Radio Burst Recognition

The automatic recognition of the Type III radio burst using Hough Transform and idea of active contour.

## Method

### Read and preprocess the data

The nancay data is in the form of binary, we developed a Matlab data driver to read in the data. And use the code [constback_sub](src/matlab/constback_sub.m) to remove the long-term background.

### Binarization

Use several ways to transform the flux intensity data into binary. ([LocalMax](src/matlab/get_local_max_map.m))

### Recognition

Use the Hough transform to recognize the time and frequency range of the Type III radio burst ([ProcEvent](src/matlab/proc_event.m))

### Identify the backbone

The Active Contour Method for the backbone of the radio burst (**ACBone**)

Take the result of the Hough transform as a initail position and iteratively move the line to find the backbone position. The demo code : ([ACBone-Matlab](src/matlab/active_contour.m), [ACBone-Python](src/python/ACBone.py))

![img](img/activecontour.GIF)

Eventually, we can obtain the centerline of a Type III radio burst.

Update: Python implementation with for LOFAR, (src/python)

![img](img/LOFAR_20220413_135000_LBA_OUTER.fits.jpg)

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
