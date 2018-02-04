# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/"
arcpy.env.overwriteOutput = True
EVI = arcpy.Raster("E:/FragEVI/processed/EVI/030005-6_2010-2012_EVI_NAD83.tif")
arcpy.env.snapRaster = EVI

# aggregate 1m data to EVI grid
grass30m = Aggregate(in_raster="E:/FragEVI/processed/bos.grass_only.tif", cell_factor="30", aggregation_type="SUM", extent_handling="EXPAND", ignore_nodata="DATA")
grass30m.save("E:/FragEVI/processed/bos.grass30m.tif")
print("snapped 1m grass cover to Landsat grid")

barr30m = Aggregate(in_raster="E:/FragEVI/processed/bos.barr_only.tif", cell_factor="30", aggregation_type="SUM", extent_handling="EXPAND", ignore_nodata="DATA")
barr30m.save("E:/FragEVI/processed/bos.barr30m.tif")
print("snapped 1m barren cover to Landsat grid")

can30m = Aggregate(in_raster="E:/FragEVI/processed/bos.can_only.tif", cell_factor="30", aggregation_type="SUM", extent_handling="EXPAND", ignore_nodata="DATA")
can30m.save("E:/FragEVI/processed/bos.can30m.tif")
print("snapped 1m canopy cover to Landsat grid")




