#### Postprocesing WRF module output with RIP on Physon cluster

The commands below describe the basic procedure for postprocessing the ``d02`` 
data from WRF model output stored in ``FC_ROOT`` using installed RIP module:


```
#load rip module
module load rip

#clean old data
rm -f WRF* *.cgm *.rgb
#convert NetCDF output to ARW
ripdp_wrfarw WRF all $FC_ROOT/wrfprd/wrfout_d02*

#generate CGM file for temperature field
rip WRF bg_temp.in
#generate RGB files
ctrans -d sun -res 800x800 bg_temp.cgm > bg_temp.rgb

#generate animated GIF sequence using ImageMagick ``convert`` tool
convert -delay 50 bg_temp_$DATE.rgb bg_temp_current_hr.gif
#clean temporary ``frames`` directory
rm -f frames/*
#convert RGB files to PNG
convert bg_temp.rgb -depth 24  frames/%03d.png
#generate MP4 and OGV require several steps and use ``pn2yuv`` and ``ffmpeg`` tools
png2yuv -I p -f 5 -b 0 -n 144 -j frames/%03d.png > temp.yuv
vpxenc --good --cpu-used=0 --auto-alt-ref=1 --end-usage=vbr --passes=2 --threads=2 --fps=5/1 --target-bitrate=3000 -o temp.webm temp.yuv
/opt/physon/bin/ffmpeg -y -r 5 -i temp.webm -c:v libx264 -preset slow -pix_fmt yuv420p -refs:v 4 -qp 4 temp.mp4
/opt/physon/bin/ffmpeg -y -r 5 -i temp.webm -c:v libtheora  -qscale:v 10 -preset veryslow -qp 0 temp.ogv


#generate CGM file for pressure and precipitation fields
rip WRF bg_precip.in
#generate RGB files
ctrans -d sun -res 800x800 bg_precip.cgm > bg_precip_$DATE.rgb

#generate animated GIF sequence using ImageMagick ``convert`` tool
convert -delay 50 bg_precip_$DATE.rgb bg_precip_current_hr.gif
#clean temporary ``frames`` directory
rm -f frames/*
#convert RGB files to PNG
convert bg_precip_$DATE.rgb -depth 24  frames/%03d.png
#generate MP4 and OGV require several steps and use ``pn2yuv`` and ``ffmpeg`` tools
png2yuv -I p -f 5 -b 0 -n 144 -j frames/%03d.png > precip.yuv
vpxenc --good --cpu-used=0 --auto-alt-ref=1 --end-usage=vbr --passes=2 --threads=2 --fps=5/1 --target-bitrate=3000 -o precip.webm precip.yuv
/opt/physon/bin/ffmpeg -y -r 5 -i precip.webm -c:v libx264 -preset slow -pix_fmt yuv420p -refs:v 4 -qp 4 precip.mp4
/opt/physon/bin/ffmpeg -y -r 5 -i precip.webm -c:v libtheora  -qscale:v 10 -preset veryslow -qp 0 precip.ogv
```
