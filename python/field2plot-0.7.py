import sys, getopt
import glob
import os
from netCDF4 import Dataset as netcdf
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from osgeo import gdal, osr
from matplotlib.image import imread
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from tzlocal import get_localzone
from dateutil import parser
import datetime
import matplotlib
from matplotlib.patches import Circle
import matplotlib.patheffects as path_effects
from matplotlib import cm
import math
from matplotlib.colors import ListedColormap
import subprocess
import wrf
import sites

def listfiles(basedir, prefix):
  files = []
  try:
#    os.chdir(basedir)
    for file in sorted(glob.glob(basedir+'/'+prefix+"*")):
      files.append(file)
  except Exception as e:
    print('Exception reading basefolder {} {}'.format(basedir,e))
  return files

def buildBaseLayer(DPI, alpha, YSIZE, extent):
 
  plt.axes().set_aspect('equal', 'datalim')

  fname = 'tmp.tif'
  img = imread(fname)
  img = img[::-1]

  #ax = plt.axes(projection=ccrs.PlateCarree())
  ax = plt.axes(projection=ccrs.Mercator())
  ax.set_extent(extent)

  plt.imshow(img, origin='lower', transform=ccrs.PlateCarree(), extent=extent, zorder=0, cmap=plt.cm.binary)


  # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
  physical_map = cfeature.NaturalEarthFeature(
    category='physical',
    name='land',
    scale='10m',
    facecolor=cfeature.COLORS['land'])

  ocean_map = cfeature.NaturalEarthFeature(
   category='physical',
   name='ocean',
   scale='10m',
   facecolor='none')

  main_rivers_map = cfeature.NaturalEarthFeature(
    category='physical',
    name='rivers_lake_centerlines',
    scale='10m',
    facecolor='none')

  rivers_map = cfeature.NaturalEarthFeature(
    category='physical',
    name='rivers_europe',
    scale='10m',
    facecolor='none')

  countries_map = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_0_countries',
    scale='10m',
    facecolor='none')

  shpfilename = shpreader.natural_earth(resolution='10m',
                                      category='cultural',
                                      name='admin_0_countries')

  reader = shpreader.Reader(shpfilename)
  countries = reader.records()

  for country in countries:
    if country.attributes['ADM0_A3'] == 'BGR':
      ax.add_geometries(country.geometry, ccrs.PlateCarree(),
               	facecolor=(0.90, 0.90, 0.80),
                label=country.attributes['ADM0_A3'],
                alpha = alpha)
    else:
      ax.add_geometries(country.geometry, ccrs.PlateCarree(),
                facecolor=(0.90, 0.90, 0.90),
                label=country.attributes['ADM0_A3'],
                alpha = alpha)

  urban_map = cfeature.NaturalEarthFeature(
    category='cultural',
    name='urban_areas',
    scale='10m',
    facecolor='gray')

  SOURCE = 'Natural Earth'
  LICENSE = 'public domain'

  ax.add_feature(countries_map, edgecolor='black', alpha=alpha)
  ax.add_feature(urban_map, edgecolor='none', alpha=0.5)
  ax.add_feature(main_rivers_map, edgecolor=cfeature.COLORS['water'])
  ax.add_feature(rivers_map, edgecolor=cfeature.COLORS['water'])

  return ax

def main(argv):
  basedir='./'
  prefix='wrfout_d02'
  field = None
  outprefix = 'out'
  ZeroK = 273.16
  levels = 31
  DPI = 100
  YSIZE = 720.
  alpha = 0.1
# Define some constants
  pq0 = 379.90516
  a2  = 17.2693882
  a3  = 273.16
  a4  = 35.86

  try:
    opts, args = getopt.getopt(argv,"h:b:p:f:",["basedir=","prefix=","field=","outprefix"])
  except getopt.GetoptError:
    print 'field2plot.py -b <basedir> ['+basedir+'] -p <prefix> ['+prefix+'] -f <field> [] -o <outprefix>['+outprefix+']'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print 'field2plot.py -b <basedir> ['+basedir+'] -p <prefix> ['+prefix+'] -f <field> [] -o <outprefix>['+outprefix+']'
      sys.exit()
    elif opt in ("-b", "--basedir"):
      basedir = arg
    elif opt in ("-f", "--field"):
      field = arg
    elif opt in ("-p", "--prefix"):
      prefix = arg
    elif opt in ("-o", "--outprefix"):
      outprefix = arg
      

  if (field == None):
    print "2D field is not supplied. Please provide WRF model field (Temp, Rain, Hail, Humid, Cloud or Solar)"
    sys.exit(1)

  flist = listfiles(basedir, prefix)

  field2D = []
  T2_field2D = []
  timestep = 0

  if not(len(flist)):
    print 'No candidates for impot files found ...'
    sys.exit(1)

  xi = []
  yi = []
  x  = []
  y  = []
  Timestamps = []
  Times = []
  minVal = 0.
  maxVal = 0.
  for i in range(0,len(flist)):
    field2D.append([])
    T2_field2D.append([])
  for file in flist:
    print 'Processing: ', file

    ncfile = netcdf(file)
    strDateTime = ncfile.variables['Times'][0].tostring().replace('_', ' ')
    local_tz = get_localzone()
    date = parser.parse(strDateTime)
    strDateTimeLocal = local_tz.localize(date);
    Times.append(strDateTimeLocal)

    strDatePrefix = local_tz.localize(date).strftime("%Y%m%d%H%M%S")
    Timestamps.append(strDatePrefix)

    if (timestep == 0):
      xlong = ncfile.variables['XLONG'][0]
      xlat = ncfile.variables['XLAT'][0]
      znu = ncfile.variables['ZNU'][0]
      west_east = len(xlong[0])
      south_north = len(xlat)
      nc_west_east = ncfile.getncattr('WEST-EAST_GRID_DIMENSION')
      nc_south_north = ncfile.getncattr('SOUTH-NORTH_GRID_DIMENSION')
      bottom_top = len(znu)
      # load projection information from the ncfile
      truelat1 = ncfile.TRUELAT1
      truelat2 = ncfile.TRUELAT2
      ref_lat  = ncfile.CEN_LAT
      ref_lon  = ncfile.CEN_LON
      stand_lon= ncfile.STAND_LON
      nc_dx = ncfile.DX
      nc_dy = ncfile.DY

      xmax = xlong[0][west_east - 1]
      xmin = xlong[0][0]
      ymin = xlat[0][0]
      ymax = xlat[south_north - 1][0]
      dx = (xmax - xmin) / (west_east-1)
      dy = (ymax - ymin) / (south_north-1)
      xi = np.arange(xmin, xmax, dx)
      yi = np.arange(ymin, ymax,dy)
      for i in range(0, south_north - 1):
        for j in range(0, west_east - 1):
          x.append(xlong[i][j])
          y.append(xlat[i][j])
 
    fieldGrid = []
    T2Grid = []
    if (field == 'Temp'):
      fieldGrid = ncfile.variables['T2'][0]
      for i in range(0, south_north - 1):
        for j in range(0, west_east - 1):          
          field2D[timestep].append(fieldGrid[i][j]-ZeroK)
    elif (field == 'Solar'):
      fieldGrid = ncfile.variables['SWDOWN'][0]
      for i in range(0, south_north - 1):
        for j in range(0, west_east - 1):
          field2D[timestep].append(fieldGrid[i][j])
    elif (field == 'Humid'):
      fieldGrid = ncfile.variables['Q2'][0]
      Q2Grid    = ncfile.variables['Q2'][0]
      T2Grid    = ncfile.variables['T2'][0]
      PSFCGrid  = ncfile.variables['PSFC'][0]
      for i in range(0, south_north - 1):
        for j in range(0, west_east - 1):
          SatVapPress = 611* math.exp(a2 * (T2Grid[i][j]-a3)/(T2Grid[i][j]-a4))
          fieldGrid[i][j] = Q2Grid[i][j] * (PSFCGrid[i][j] - SatVapPress) / (0.622 * SatVapPress)
          field2D[timestep].append(fieldGrid[i][j])
    elif (field == 'Rain'):
      T2Grid    = ncfile.variables['T2'][0]
      if (timestep == 0):
        fieldGrid = ncfile.variables['RAINNC'][0]+ncfile.variables['SNOWNC'][0]+ncfile.variables['GRAUPELNC'][0]+ncfile.variables['HAILNC'][0]
      else:
        fieldGrid = ncfile.variables['RAINNC'][0]+ncfile.variables['SNOWNC'][0]+ncfile.variables['GRAUPELNC'][0]+ncfile.variables['HAILNC'][0]
        for i in range(0, south_north - 1):
          for j in range(0, west_east - 1):
            fieldGrid[i][j] -= Buffer[i][j]
      Buffer = ncfile.variables['RAINNC'][0]+ncfile.variables['SNOWNC'][0]+ncfile.variables['GRAUPELNC'][0]+ncfile.variables['HAILNC'][0] 
      for i in range(0, south_north - 1):
        for j in range(0, west_east - 1):
          field2D[timestep].append(fieldGrid[i][j])
          T2_field2D[timestep].append(T2Grid[i][j]-ZeroK)
    elif (field == 'Snow'):
      if (timestep == 0):
        fieldGrid = ncfile.variables['SNOWNC'][0]+ncfile.variables['GRAUPELNC'][0]
      else:
        fieldGrid = ncfile.variables['SNOWNC'][0]+ncfile.variables['GRAUPELNC'][0]
        for i in range(0, south_north - 1):
          for j in range(0, west_east - 1):
            fieldGrid[i][j] -= Buffer[i][j]
      Buffer = ncfile.variables['SNOWNC'][0]+ncfile.variables['GRAUPELNC'][0]
      for i in range(0, south_north - 1):
        for j in range(0, west_east - 1):
          field2D[timestep].append(fieldGrid[i][j])
    elif (field == 'Hail'):
      if (timestep == 0):
        fieldGrid = ncfile.variables['HAILNC'][0]+ncfile.variables['GRAUPELNC'][0]
      else:
        fieldGrid = ncfile.variables['HAILNC'][0]+ncfile.variables['GRAUPELNC'][0]
        for i in range(0, south_north - 1):
          for j in range(0, west_east - 1):
            fieldGrid[i][j] -= Buffer[i][j]
      Buffer = ncfile.variables['HAILNC'][0]+ncfile.variables['GRAUPELNC'][0] 
      for i in range(0, south_north - 1):
        for j in range(0, west_east - 1):
          field2D[timestep].append(fieldGrid[i][j])
    elif (field == 'Cloud'):
      fieldGrid = ncfile.variables['CLDFRA'][0][0]
      for i in range(1, bottom_top - 1):
         fieldGrid += ncfile.variables['CLDFRA'][0][i]
      #fieldGrid *= 0
      for i in range(0, south_north - 1):
        for j in range(0, west_east - 1):
          field2D[timestep].append(fieldGrid[i][j])
      #print('ts: {} min: {} max: {}'.format(timestep, min(field2D[timestep]), max(field2D[timestep])))
    else:
      print ("2D field {} is not supplied. Please provide WRF model field (Temp, Rain, Cloud or Solar)".format(field))
      sys.exit(1)

    timestep +=1

  extent = [xlong[0][0],  xlong[0][west_east-1], xlat[0][west_east-1], xlat[south_north-1][0]]
  minVal = min(min(field2D[0:timestep-1]))
  maxVal = max(max(field2D[0:timestep-1]))
  #minVal = min(min(field2D[0:]))
  #maxVal = max(max(field2D[0:]))

  if (field == 'Rain' or field == 'Hail' or field == 'Snow'):
    levels = 10
    clevs = np.logspace(math.log10(0.005*maxVal), math.log10(maxVal), levels, endpoint=True)
    if (field == 'Rain'):
      cmapBase = plt.cm.Blues
    else:
      cmapBase = plt.cm.BuPu
    new_cmap = cmapBase(np.arange(cmapBase.N))
    new_cmap[0] = 0.0
    cmap = ListedColormap(new_cmap)

    alpha = 0.55
  elif (field == 'Humid'):
    levels = 10
    clevs = np.linspace(0.0, 1.0, levels, endpoint=True)
    cmapBase = plt.cm.gray
    new_cmap = cmapBase(np.arange(cmapBase.N))
    for i in range(0, len(new_cmap)):
      if (i < int(round(0.90*len(new_cmap)))):
        new_cmap[i][3] = 0.0

    cmap = ListedColormap(new_cmap)

    alpha = 0.55
  elif (field == 'Cloud'):
    levels = 10
    clevs = np.linspace(0.0, maxVal, levels, endpoint=True)
    cmapBase = plt.cm.gray
    new_cmap = cmapBase(np.arange(cmapBase.N))
    for i in range(0, len(new_cmap)):
      if (i < int(round(0.2*len(new_cmap)))):
        new_cmap[i][3] = 0.0
        new_cmap[i][0] = 1.0 
        new_cmap[i][1] = 1.0 
        new_cmap[i][2] = 1.0 
      else:
        new_cmap[i][3] = new_cmap[i][0]
        new_cmap[i][0] = 1.0 
        new_cmap[i][1] = 1.0
        new_cmap[i][2] = 1.0

    cmap = ListedColormap(new_cmap)

    alpha = 0.55
  elif (field == 'Temp'):
    clevs = np.linspace(-40, 40, levels, endpoint=True)
    cmap = plt.cm.jet
    alpha = 0.1
  elif (field == 'Solar'):
    clevs = np.linspace(0.0, 1000, levels, endpoint=True)
    cmap = plt.cm.jet
    alpha = 0.1
  else:
    clevs = np.linspace(minVal, maxVal, levels, endpoint=True)
    cmap = plt.cm.jet
  
# crop the image

  print('Domain boundaries xmin: {} ymin: {} xmax: {} ymax: {}'.format(xmin, ymin, xmax, ymax))
  cmd = ['gdalwarp','-overwrite','-t_srs','EPSG:4326', '-te', str(xmin), str(ymin), str(xmax), str(ymax),'NE2_HR_LC_SR_W/NE2_HR_LC_SR_W.tif','tmp.tif']
  proc = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  stdout,stderr=proc.communicate()
  exit_code=proc.wait()

  if exit_code: #Oops, something went wrong!
    raise RuntimeError(stderr)
  else:print stdout
 
  for ts in range(1, timestep):
    fname = field+'_'+str(ts-1).zfill(4)+'.png'

    print fname

    GD = interpolate.griddata((np.ravel(x), np.ravel(y)),  np.ravel(field2D[ts]), (xi[None,:], yi[:,None]), method='cubic', rescale=True)

    fig = plt.figure(dpi=DPI, frameon=False)
    figSize = fig.get_size_inches()

    ax = buildBaseLayer(DPI, alpha, YSIZE, extent)
    
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi

    if (field == 'Rain' or field == 'Hail' or field == 'Snow' or field == 'Humid' or field == 'Cloud'):
      #for contour in GD:
      #  for tpl in contour:
      #    print('size: {}'.format(len(tpl)))
      plt.contourf(xi, yi, GD, clevs,transform=ccrs.PlateCarree(),cmap=cmap,extend='both', antialiased=True, linecolor='none',zorder=2)
    else:
      plt.contourf(xi, yi, GD, clevs,transform=ccrs.PlateCarree(),cmap=cmap,extend='both', antialiased=True, linecolor='none',alpha=1.)

    #Add temperatures for Rain plot
    if (field == 'Rain'):
      color_norm = cm.colors.Normalize(vmin=-35, vmax=45)  
      for site in sites.sites:
        site_name = site['name']
        site_lat = site['lat']
        site_lon = site['lon']
        indx = wrf.ll_to_ij(
          1,
          truelat1,
          truelat2,
          stand_lon,
          nc_dx,
          nc_dy,
          ref_lat,
          ref_lon,
          site_lat,
          site_lon
        )
	j0 = nc_west_east / 2 + indx[0] - 1
        i0 = nc_south_north / 2 + indx[1] - 1
	x0 = ( j0 * dx ) / ( xmax - xmin )
        y0 = ( i0 * dy ) / ( ymax - ymin )
        xx0 = x0 * 1.2
        yy0 = y0 * 1.1042
        temperature = T2_field2D[ts][j0 + i0 * (west_east - 1)]

        circle = Circle(
          xy=[xx0  * bbox.width, yy0 * bbox.height ],
          radius=0.18,
          color=cm.colors.rgb2hex(
            cm.jet(color_norm(temperature),bytes=False)
          ),
          alpha=0.4,
          transform=fig.dpi_scale_trans,
          zorder=3
        )
        ax.add_patch(circle)
        ax.text(x0, y0, "{0:.0f}".format(temperature),
          #horizontalalignment=site['halign'],
          #verticalalignment=site['valign'],
          horizontalalignment='center',
          verticalalignment='center',
          transform=ax.transAxes,
          weight = 'bold',
          color = '#393a3b',
          fontsize = 'x-large',
          zorder=4)
        ax.text(x0 + 0.2 / bbox.width, y0, "{}".format(site['cname']),
          #horizontalalignment=site['halign'],
          #verticalalignment=site['valign'],
          horizontalalignment='left',
          verticalalignment='center',
          transform=ax.transAxes,
          weight = 'bold',
          color = '#efefef',
          fontsize = 'large',
          zorder=5).set_path_effects([
            path_effects.Stroke(linewidth=1, foreground='black'),
            path_effects.Normal()
          ])
          

    ax.text(0.95, 0.95, Times[ts],
        horizontalalignment='right',
        verticalalignment='top',
        transform=ax.transAxes,
        weight = 'bold',
        color = '#393a3b',
        fontsize = 'x-large',
        zorder=5)
    ax.text(0.95, 0.05, 'Nestum@SofiaTech',
        horizontalalignment='right',
        verticalalignment='bottom',
        transform=ax.transAxes,
        weight = 'bold',
        color = '#393a3b',
        fontsize = 'medium',
        zorder=5)

    fig.set_size_inches(figSize[0]*YSIZE/(figSize[1]*DPI), YSIZE/DPI)
    fig.savefig(fname, bbox_inches='tight', pad_inches=0)
    plt.close(fig)

if __name__ == "__main__":
   main(sys.argv[1:])
