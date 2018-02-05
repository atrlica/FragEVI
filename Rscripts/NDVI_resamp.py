# import system modules 
import arcpy, arcinfo
from arcpy import env

# Set environment settings
arcpy.env.workspace = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/"
arcpy.env.overwriteOutput = True
canR = arcpy.Raster("E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif")
arcpy.env.snapRaster = canR

# resample 2.4m NDVI raster to 1m canopy grid
print("resampling NDVI raster to 1m grid")
ndviRes = arcpy.Resample_management(in_raster="E:/FragEVI/data/NDVI/NDVI.img", out_raster="/NDVI_Resample", cell_size="1", resampling_type="BILINEAR")
ndviRes.save("E:/FragEVI/data/NDVI/NDVI_1m_res_cangrid.tif")
print("wrote resampled NDVI to disk")
