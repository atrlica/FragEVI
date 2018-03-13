# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/"
arcpy.env.overwriteOutput = True
try:
    snapR = arcpy.Raster("E:/FragEVI/processed/boston/bos.barr.tif")
    print("found 1m raster for grid guide")
except Exception as e:
    print(e)

arcpy.env.snapRaster = snapR

arcpy.Resample_management(in_raster="E:/FragEVI/processed/boston/bos.isa.rereg.tif", out_raster="E:/FragEVI/processed/boston/bos.isa.RRR.tif", cell_size=1, resampling_type="MAJORITY")
