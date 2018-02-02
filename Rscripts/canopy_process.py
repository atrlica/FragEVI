# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")

# Set environment settings
arcpy.env.workspace = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/"
arcpy.env.overwriteOutput = True

dbpath = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/"

# convert 1m 0/1 canopy raster to poly
canRastIn = "E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif"
arcpy.RasterToPolygon_conversion(in_raster=canRastIn, out_polygon_features=dbpath + "bos_can_P", simplify="NO_SIMPLIFY", raster_field="VALUE")
print("converted canopy raster to polygon")

# select polygons that are 0 canopy
arcpy.Select_analysis(in_features=dbpath + "bos_can_P", out_feature_class=dbpath+"bos_can_0", where_clause="gridcode = 0")
print("isolated canopy gaps")

# select gaps that are larger than gap size threshold
gapSizeThresh = 50
arcpy.Select_analysis(in_features=dbpath+"bos_can_0", out_feature_class=dbpath+"bos_can_0_filt", where_clause="Shape_Area >= "+str(gapSizeThresh))
print("removed canopy gaps below " + str(gapSizeThresh) + " m2 in area")

# convert filtered gap polygon to Raster
arcpy.PolygonToRaster_conversion(in_features=dbpath+"bos_can_0_filt", value_field="gridcode", out_rasterdataset=dbpath+"bos_can_0_R", cell_assignment="CELL_CENTER", priority_field="gridcode", cellsize="1")
print("converted filtered gap polys to raster")

# Create expand map of 10, 20, 30 m buffers around gaps
buffs = [10, 20, 30]
rasterDump = "F:/FragEVI/processed/"
for b in range(0,2):
    try:
        buffR = arcpy.gp.Expand_sa(dbpath+"bos_can_0_R", dbpath+"nocan_"+str(buffs[b])+"mbuff", str(buffs[b]), "0")
        buffR.save("E:/FragEVI/processed/nocan_"+str(buffs[b])+"mbuff.tif")
        # arcpy.RasterToOtherFormat(dbpath+"nocan_"+str(buffs[b])+"mbuff", rasterDump+"nocan_"+str(buffs[b])+"mbuff", "TIFF")
        # arcpy.CopyRaster_management(dbpath+"nocan_"+str(buffs[b])+"mbuff", "E:/FragEVI/processed/nocan_"+str(buffs[b])+"mbuff.tif")
        print("did " + str(buffs[b]) + "m buffer and wrote raster .tif to /processed")
    except Exception as e:
        print(e)




