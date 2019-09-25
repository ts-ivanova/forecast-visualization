#### Python WRF portprocessing

This piece of code is written by Stoyan Pisov and can be used 
if you ask him nicely: pisov<at>phys.uni-sofia.bg

Thee python script read the generated NetCDFv4 data files from WRF model
and visualise the following fields:

Temp  -> Temperature [T2]
Rain  -> Rainfall difference between two data files RAINC[t]+RAINNC[t]+RAINSH[t]-RAINC[t-1]-RAINNC[t-1]-RAINSH[t-1]
Hail  -> Hail difference between two data files GRAUPELNC[t]+HAILNC[t]-GRAUPELNC[t-1]-HAILNC[t-1]
Humid -> Calculated Humidity above 95% (Fog prediction) 


The script is tuned and tested for Bulgaria but ... in theory should render any domain.
You may need to install additional python packages such as:

```
sudo apt-get install python-netcdf4 python-pip python-gdal

pip install cartopy
pip install matplotlib
pip install scipy
pip install tzlocal
pip install pillow
```

(already installed on nestum cluster)

The script also require Natural Earth raster images archive NE2_HR_LC_SR_W.zip 
availabale at http://kelsocartography.com/downloads/gis/2012/natural_earth/ne_2.0.0_rc2/natural_earth_2.0.0_rc2_raster/

Basic usage:

Help:

```
python ./field2plot-0.1.py -h
```

by deafult script will lok in current workinf directory "./" for data files with prefix "wrfout_d02" and will render the 
specified field by -f [Temp, Rain, Hail, Humid]

```
python ./field2plot-0.1.py -f Temp

python ./field2plot-0.1.py -f Rain
```

!!!Important!!!

Script is expecting at least two data files in order to make the sequence
