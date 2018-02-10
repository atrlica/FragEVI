# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/"
arcpy.env.overwriteOutput = True
try:
    snapR = arcpy.Raster("E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif")
    print("found 1m raster for grid guide")
except Exception as e:
    print(e)

arcpy.env.snapRaster = snapR

# resample NDVI to canopy grid
arcpy.Resample_management(in_raster="E:/FragEVI/data/NDVI/NDVI.img", out_raster="E:/FragEVI/data/NDVI_1m_res_cangrid.tif/", cell_size=1, resampling_type="BILINEAR")
