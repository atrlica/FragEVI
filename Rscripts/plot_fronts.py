# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "E:/FragEVI/FragEVIworking.gdb/"
arcpy.env.overwriteOutput = True

# can set up a loop to create all these files for all the Rope and all the Andy plots
#plot R18 -- RIGHT
try:
    front="E:/FragEVI/processed/zones/R18_front.shp"
    print("found front of plot R18")
except Exception as e:
    print(e)


buffs = [10, 20, 30, 40, 50]
for b in range(0,5):
    try:
        distField = str(buffs[b])+" Meters"
        if arcpy.Exists("E:/FragEVI/processed/zones/R18_"+str(buffs)+"m.shp"):
           arcpy.Delete_management("E:/FragEVI/processed/zones/R18_"+str(buffs)+"m.shp")
        frontBuff = arcpy.Buffer_analysis(in_features=front, out_feature_class="E:/FragEVI/processed/zones/R18_"+str(buffs[b])+"m.shp", buffer_distance_or_field=distField, line_side="RIGHT", line_end_type="FLAT", dissolve_option="NONE", dissolve_field="", method="PLANAR")
        print("did " + str(buffs[b]) + "m front buffer and wrote .shp to /processed/zones")
    except Exception as e:
        print(e)

# plot R19 -- LEFT
try:
    front="E:/FragEVI/processed/zones/R19_front.shp"
    print("found front of plot R19")
except Exception as e:
    print(e)


buffs = [10, 20, 30, 40, 50]
for b in range(0,5):
    try:
        distField = str(buffs[b])+" Meters"
        if arcpy.Exists("E:/FragEVI/processed/zones/R19_"+str(buffs)+"m.shp"):
           arcpy.Delete_management("E:/FragEVI/processed/zones/R19_"+str(buffs)+"m.shp")
        frontBuff = arcpy.Buffer_analysis(in_features=front, out_feature_class="E:/FragEVI/processed/zones/R19_"+str(buffs[b])+"m.shp", buffer_distance_or_field=distField, line_side="LEFT", line_end_type="FLAT", dissolve_option="NONE", dissolve_field="", method="PLANAR")
        print("did " + str(buffs[b]) + "m front buffer and wrote .shp to /processed/zones")
    except Exception as e:
        print(e)


# plot R20 -- LEFT
try:
    front="E:/FragEVI/processed/zones/R20_front.shp"
    print("found front of plot R20")
except Exception as e:
    print(e)


buffs = [10, 20, 30, 40, 50]
for b in range(0,5):
    try:
        distField = str(buffs[b])+" Meters"
        if arcpy.Exists("E:/FragEVI/processed/zones/R20_"+str(buffs)+"m.shp"):
           arcpy.Delete_management("E:/FragEVI/processed/zones/R20_"+str(buffs)+"m.shp")
        frontBuff = arcpy.Buffer_analysis(in_features=front, out_feature_class="E:/FragEVI/processed/zones/R20_"+str(buffs[b])+"m.shp", buffer_distance_or_field=distField, line_side="LEFT", line_end_type="FLAT", dissolve_option="NONE", dissolve_field="", method="PLANAR")
        print("did " + str(buffs[b]) + "m front buffer and wrote .shp to /processed/zones")
    except Exception as e:
        print(e)

