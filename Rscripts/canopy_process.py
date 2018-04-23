# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "E:/FragEVI/FragEVIworking.gdb"
arcpy.env.overwriteOutput = True

try:
    can01Filt = arcpy.Raster("E:/FragEVI/processed/boston/bos_can01_filt.tif")
    print("found filtered gap raster, performing Expand")
except:
    print("filtering canopy gaps for area")
    # convert 1m 0/1 canopy raster to poly
    canRastIn = "E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif"
    arcpy.RasterToPolygon_conversion(in_raster=canRastIn, out_polygon_features="bos_can_P", simplify="NO_SIMPLIFY", raster_field="VALUE")
    print("converted canopy raster to polygon")

    # select polygons that are gaps(0) + larger than threshold size, or canopy(1)
    gapSizeThresh = 50
    arcpy.Select_analysis(in_features="bos_can_P", out_feature_class="bos_can_01_filt", where_clause="gridcode = 1 OR (gridcode=0 AND Shape_Area>"+str(gapSizeThresh)+")")
    print("eliminated small canopy gaps")

    # convert filtered polygon to Raster
    arcpy.PolygonToRaster_conversion(in_features="bos_can_01_filt", value_field="gridcode", out_rasterdataset="bos_can_01_R", cell_assignment="CELL_CENTER", priority_field="gridcode", cellsize="1")
    arcpy.CopyRaster_management(in_raster = "bos_can_01_R", out_rasterdataset = "E:/FragEVI/processed/boston/bos_can01_filt.tif", format = "TIFF")
    print("converted filtered gap polys to raster")
    can01Filt = arcpy.Raster("E:/FragEVI/processed/boston/bos_can01_filt.tif")


# Create expand map of 10, 20, 30 m buffers around gaps
# get expand map of all distance classes by 1m increments 0-100m
buffs = range(1,101)
for b in range(0,len(buffs)):
    try:
        buffR = Expand(in_raster=can01Filt, number_cells=buffs[b], zone_values=0)
        buffR.save("E:/FragEVI/processed/boston/bos.nocan_"+str(buffs[b])+"mbuff.tif")
        print("did " + str(buffs[b]) + "m buffer and wrote raster .tif to /processed")
    except Exception as e:
        print(e)

