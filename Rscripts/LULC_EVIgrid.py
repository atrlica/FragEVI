# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/"
arcpy.env.overwriteOutput = True
try:
    snapR = arcpy.Raster("E:/FragEVI/processed/EVI/030005-6_2010-2012_EVI_NAD83.tif")
    print("found 30m raster for grid guide")
except Exception as e:
    print(e)

arcpy.env.snapRaster = snapR

# rasterize LULC polygons using maximum combined area, 30m EVI grid
arcpy.PolygonToRaster_conversion(in_features="E:/FragEVI/data/LULC/LU_polys_NAD83UTM19N/LU_polys_NAD83UTM19N.shp", value_field="LUCODE", out_rasterdataset="E:/FragEVI/processed/LULC30m.tif", cell_assignment="MAXIMUM_COMBINED_AREA", priority_field="NONE", cellsize="30")
print("rasterized LULC polygons to EVI layer")
