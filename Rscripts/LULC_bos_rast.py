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
csRes = arcpy.GetRasterProperties_management(in_raster="E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif", property_type="CELLSIZEX", band_index="")
cs = csRes.getOutput(0)
print("cellsizex="+str(cs))

arcpy.Select_analysis(in_features="E:/FragEVI/data/LULC/LU_polys_NAD83UTM19N/LU_polys_NAD83UTM19N.shp", out_feature_class="E:/FragEVI/processed/boston/LU_polys_bos.shp", where_clause=""""TOWN" = 'BOSTON'""")
print("selected Boston LULC polys")
arcpy.PolygonToRaster_conversion(in_features="E:/FragEVI/processed/boston/LU_polys_bos.shp", value_field="LUCODE", out_rasterdataset="E:/FragEVI/processed/boston/LU_bos_r1m.tif", cell_assignment="CELL_CENTER", priority_field="NONE", cellsize=str(cs))
print("rasterized to grid "+str(cs)+"m")
