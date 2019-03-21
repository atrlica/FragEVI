# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "E:/FragEVI/FragEVIworking.gdb"
arcpy.env.overwriteOutput = True

### prepare an outer buffer for roads of class 2 3 4 5 
##arcpy.Select_analysis(in_features="EOTROADS_ARC", out_feature_class="E:/FragEVI/FragEVIworking.gdb/bos_roads5", where_clause=""""MGIS_TOWN" = 'BOSTON' AND "CLASS" = 5""")
##arcpy.Select_analysis(in_features="EOTROADS_ARC", out_feature_class="E:/FragEVI/FragEVIworking.gdb/bos_roads234", where_clause=""""CLASS" IN ( 2, 3, 4 ) AND "MGIS_TOWN" = 'BOSTON'""")
##arcpy.Buffer_analysis(in_features="bos_roads234", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road234_buff10", buffer_distance_or_field="10 Meters", line_side="FULL", line_end_type="ROUND", dissolve_option="ALL", dissolve_field="", method="PLANAR")
##arcpy.Buffer_analysis(in_features="bos_roads5", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road5_buff6", buffer_distance_or_field="6 Meters", line_side="FULL", line_end_type="ROUND", dissolve_option="ALL", dissolve_field="", method="PLANAR")
##arcpy.Union_analysis(in_features="road5_buff6 #;road234_buff10 #", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road2345_buffs", join_attributes="ALL", cluster_tolerance="", gaps="GAPS")
##arcpy.Dissolve_management(in_features="road2345_buffs", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road2345_buffsD", dissolve_field="", statistics_fields="", multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES")
print("have outer buffers for roads of class 234 and 5, joined")

# prepare inner buffer width = 2m
##arcpy.Buffer_analysis(in_features="bos_roads234", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road234_buff8", buffer_distance_or_field="8 Meters", line_side="FULL", line_end_type="ROUND", dissolve_option="ALL", dissolve_field="", method="PLANAR")
##arcpy.Buffer_analysis(in_features="bos_roads5", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road5_buff4", buffer_distance_or_field="4 Meters", line_side="FULL", line_end_type="ROUND", dissolve_option="ALL", dissolve_field="", method="PLANAR")
##arcpy.Union_analysis(in_features="road5_buff4 #;road234_buff8 #", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road2345_Ibuffs", join_attributes="ALL", cluster_tolerance="", gaps="GAPS")
##arcpy.Dissolve_management(in_features="road2345_Ibuffs", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road2345_IbuffsD", dissolve_field="", statistics_fields="", multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES")
print("have inner buffers for roads of class 234 and 5, joined")

## project the buffer polys to NAD83 to match canopy raster
##arcpy.Project_management(in_dataset="E:/FragEVI/FragEVIworking.gdb/road2345_buffsD", out_dataset="E:/FragEVI/FragEVIworking.gdb/road2345_buffsD_NAD83", out_coor_system="PROJCS['NAD_1983_UTM_Zone_19N',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-69.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]", transform_method="", in_coor_system="PROJCS['NAD_1983_StatePlane_Massachusetts_Mainland_FIPS_2001',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',200000.0],PARAMETER['False_Northing',750000.0],PARAMETER['Central_Meridian',-71.5],PARAMETER['Standard_Parallel_1',41.71666666666667],PARAMETER['Standard_Parallel_2',42.68333333333333],PARAMETER['Latitude_Of_Origin',41.0],UNIT['Meter',1.0]]", preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")
##arcpy.Project_management(in_dataset="E:/FragEVI/FragEVIworking.gdb/road2345_IbuffsD", out_dataset="E:/FragEVI/FragEVIworking.gdb/road2345_IbuffsD_NAD83", out_coor_system="PROJCS['NAD_1983_UTM_Zone_19N',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-69.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]", transform_method="", in_coor_system="PROJCS['NAD_1983_StatePlane_Massachusetts_Mainland_FIPS_2001',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',200000.0],PARAMETER['False_Northing',750000.0],PARAMETER['Central_Meridian',-71.5],PARAMETER['Standard_Parallel_1',41.71666666666667],PARAMETER['Standard_Parallel_2',42.68333333333333],PARAMETER['Latitude_Of_Origin',41.0],UNIT['Meter',1.0]]", preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")
print("projected buffer polys to NAD83")

## convert both polygons to rasters snapped to the canopy raster grid
##arcpy.env.snapRaster = "E:/FragEVI/processed/boston/bos.can.redux.tif"
##arcpy.PolygonToRaster_conversion(in_features="road2345_buffsD_NAD83", value_field="OBJECTID", out_rasterdataset="E:/FragEVI/processed/boston/newplanting/road2345_buffD_NAD83_R", cell_assignment="MAXIMUM_COMBINED_AREA", priority_field="NONE", cellsize="1")
##arcpy.PolygonToRaster_conversion(in_features="road2345_IbuffsD_NAD83", value_field="OBJECTID", out_rasterdataset="E:/FragEVI/processed/boston/newplanting/road2345_IbuffsD_NAD83_R", cell_assignment="MAXIMUM_COMBINED_AREA", priority_field="NONE", cellsize="1")
print("converted polygon buffers to raster at canopy grid")

## (in R) reprocess the overlaid buffer rasters to produce a single 1/0 "donut" 1m raster broken up by the 4m buffer canopy buffer and eliminating non dev/residential buffers
## raster is "processed/boston/newplanting/road2345_donut_R.tif"
print("skipping over R script that converts buffers to suitable planting donut strips")

## convert donut raster to polygons, convert to line object (linear planting distance derivative)
arcpy.Int_3d(in_raster_or_constant="E:/FragEVI/processed/boston/newplanting/road2345_donut_R.tif", out_raster="E:/FragEVI/FragEVIworking.gdb/road2345_donut_Rint")
arcpy.RasterToPolygon_conversion(in_raster="E:/FragEVI/FragEVIworking.gdb/road2345_donut_Rint", out_polygon_features="E:/FragEVI/FragEVIworking.gdb/road2345_donut_polySimp", simplify="SIMPLIFY", raster_field="Value")
print("converted donut raster to poly")

## filter polygons for too small areas, convert to line features, and export data attribute table
arcpy.Select_analysis(in_features="E:/FragEVI/FragEVIworking.gdb/road2345_donut_polySimp", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road2345_donut_polySimp_filt", where_clause="Shape_Area > 2")
arcpy.PolygonToLine_management(in_features="road2345_donut_polySimp_filt", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road2345_donut_line", neighbor_option="IGNORE_NEIGHBORS")
print("donut polys filtered for size and converted to lines") 
