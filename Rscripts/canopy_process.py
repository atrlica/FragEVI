# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "E:/FragEVI/FragEVIworking.gdb"
arcpy.env.overwriteOutput = True

try:
##    can01Filt = arcpy.Raster("E:/FragEVI/processed/boston/bos_can01_filt.tif")
    can01Filt = arcpy.Raster("E:/FragEVI/processed/boston/bos_can01_filt_redux.tif")
    print("found filtered gap raster, performing Expand")
except:
    print("filtering canopy gaps for area")
    canRastIn = "E:/FragEVI/processed/boston/bos.can.redux.tif"
    
    ## convert the SOB to an integer raster god Arc why are you like this? :-/
    print("have raw raster file, commensing integer conversion")
##    arcpy.gp.Int_sa(canRastIn, "E:/FragEVI/FragEVIworking.gdb/canR_int")
    print("finished tedious reprocessing to integer raster")
    
    # convert 1m 0/1 canopy raster to poly
    print("integer conversion complete, converting to (simplified?) raster")
####    canRastIn = "E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif"
##    canRastIn = "E:/FragEVI/FragEVIworking.gdb/canR_int"
####    arcpy.RasterToPolygon_conversion(in_raster=canRastIn, out_polygon_features="bos_can_P", simplify="NO_SIMPLIFY", raster_field="VALUE")
##    arcpy.RasterToPolygon_conversion(in_raster="canR_int", out_polygon_features="E:/FragEVI/FragEVIworking.gdb/can_red_P", simplify="SIMPLIFY", raster_field="Value")
    print("converted canopy raster to (simplified?) polygon")

    # select polygons that are gaps(0) + larger than threshold size, or canopy(1)
    print("commensing filtering of small canopy gaps")
    gapSizeThresh = 50
##    arcpy.Select_analysis(in_features="bos_can_P", out_feature_class="bos_can_01_filt", where_clause="gridcode = 1 OR (gridcode=0 AND Shape_Area>"+str(gapSizeThresh)+")")
    arcpy.Select_analysis(in_features="can_red_P", out_feature_class="can_red_filtP", where_clause="gridcode = 1 OR (gridcode=0 AND Shape_Area>"+str(gapSizeThresh)+")")
    print("eliminated small canopy gaps")

    # convert filtered polygon to Raster
    print("converting filtered canopy polygon to raster")
##    arcpy.PolygonToRaster_conversion(in_features="bos_can_01_filt", value_field="gridcode", out_rasterdataset="bos_can_01_R", cell_assignment="CELL_CENTER", priority_field="gridcode", cellsize="1")
    arcpy.PolygonToRaster_conversion(in_features="can_red_filtP", value_field="gridcode", out_rasterdataset="can_red_filtR", cell_assignment="CELL_CENTER", priority_field="gridcode", cellsize="1")
##    arcpy.CopyRaster_management(in_raster = "bos_can_01_R", out_rasterdataset = "E:/FragEVI/processed/boston/bos_can01_filt.tif", format = "TIFF")
    arcpy.CopyRaster_management(in_raster = "can_red_filtR", out_rasterdataset = "E:/FragEVI/processed/boston/bos_can01_filt_redux.tif", format = "TIFF")
    print("converted filtered canopy polygon to raster")
##    can01Filt = arcpy.Raster("E:/FragEVI/processed/boston/bos_can01_filt.tif")
    can01Filt = arcpy.Raster("E:/FragEVI/processed/boston/bos_can01_filt_redux.tif")


# Create expand map of canopy buffers by 1 m out to 100 m
# get expand map of all distance classes by 1m increments 0-100m
print("calculating 1m buffers for filtered canopy raster")
arcpy.env.snapRaster = "E:/FragEVI/processed/boston/bos.can.redux.tif"
# buffs = range(1,101)
buffs = [3,30]
for b in range(0,len(buffs)):
    try:
        buffR = Expand(in_raster=can01Filt, number_cells=buffs[b], zone_values=0)
        buffR.save("E:/FragEVI/processed/boston/bos.nocan_"+str(buffs[b])+"mbuff.tif")
        print("did " + str(buffs[b]) + "m buffer and wrote raster .tif to /processed")
    except Exception as e:
        print(e)

