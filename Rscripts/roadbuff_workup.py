# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "E:/FragEVI/FragEVIworking.gdb"
arcpy.env.overwriteOutput = True

# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "EOTROADS_ARC"
arcpy.Select_analysis(in_features="EOTROADS_ARC", out_feature_class="E:/FragEVI/FragEVIworking.gdb/bos_roads5", where_clause=""""MGIS_TOWN" = 'BOSTON' AND "CLASS" = 5""")

# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "EOTROADS_ARC"
arcpy.Select_analysis(in_features="EOTROADS_ARC", out_feature_class="E:/FragEVI/FragEVIworking.gdb/bos_roads234", where_clause=""""CLASS" IN ( 2, 3, 4 ) AND "MGIS_TOWN" = 'BOSTON'""")

# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "bos_roads234"
arcpy.Buffer_analysis(in_features="bos_roads234", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road234_buff10", buffer_distance_or_field="10 Meters", line_side="FULL", line_end_type="ROUND", dissolve_option="ALL", dissolve_field="", method="PLANAR")

# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "bos_roads5"
arcpy.Buffer_analysis(in_features="bos_roads5", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road5_buff6", buffer_distance_or_field="6 Meters", line_side="FULL", line_end_type="ROUND", dissolve_option="ALL", dissolve_field="", method="PLANAR")

arcpy.Union_analysis(in_features="road5_buff6 #;road234_buff10 #", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road2345_buffs", join_attributes="ALL", cluster_tolerance="", gaps="GAPS")

# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "road2345_buffs"
arcpy.Dissolve_management(in_features="road2345_buffs", out_feature_class="E:/FragEVI/FragEVIworking.gdb/road2345_buffsD", dissolve_field="", statistics_fields="", multi_part="MULTI_PART", unsplit_lines="DISSOLVE_LINES")





