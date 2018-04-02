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

plots = ["HF1", "HF2", "HF3", "HF4", "HF5", "HF6"]
side = ["LEFT", "LEFT", "LEFT", "RIGHT", "LEFT", "LEFT"]
buffs = [10, 20, 30, 40, 50]
nfronts = len(plots)
for a in range(0,nfronts):
    try:
        front="E:/FragEVI/processed/zones/"+str(plots[a])+"_front.shp"
        print("found front of plot "+str(plots[a]))
        for b in range(0,5):
              distField = str(buffs[b])+" Meters"
              if arcpy.Exists("E:/FragEVI/processed/zones/"+str(plots[a])+"_"+str(buffs[b])+"m.shp"):
                  arcpy.Delete_management("E:/FragEVI/processed/zones/"+str(plots[a])+"_"+str(buffs[b])+"m.shp")
              frontBuff = arcpy.Buffer_analysis(in_features=front, out_feature_class="E:/FragEVI/processed/zones/"+str(plots[a])+"_"+str(buffs[b])+"m.shp", buffer_distance_or_field=distField, line_side=str(side[a]), line_end_type="FLAT", dissolve_option="NONE", dissolve_field="", method="PLANAR")
              print("did " + str(buffs[b]) + "m front buffer and wrote .shp to /processed/zones")
    except Exception as e:
        print(e)

plots = ["R1", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11", "R13", "R14", "R15", "R16", "R17", "R18", "R19", "R20", "R21", "R22"]
side = ["LEFT", "RIGHT", "RIGHT", "LEFT", "LEFT", "RIGHT", "RIGHT", "LEFT", "LEFT", "LEFT", "RIGHT", "RIGHT", "RIGHT", "RIGHT", "LEFT", "RIGHT", "LEFT", "LEFT", "LEFT", "LEFT",]
buffs = [10, 20, 30, 40, 50]
nfronts = len(plots)
for a in range(0,nfronts):
    try:
        front="E:/FragEVI/processed/zones/"+str(plots[a])+"_front.shp"
        print("found front of plot "+str(plots[a]))
        for b in range(0,5):
              distField = str(buffs[b])+" Meters"
              if arcpy.Exists("E:/FragEVI/processed/zones/"+str(plots[a])+"_"+str(buffs[b])+"m.shp"):
                  arcpy.Delete_management("E:/FragEVI/processed/zones/"+str(plots[a])+"_"+str(buffs[b])+"m.shp")
              frontBuff = arcpy.Buffer_analysis(in_features=front, out_feature_class="E:/FragEVI/processed/zones/"+str(plots[a])+"_"+str(buffs[b])+"m.shp", buffer_distance_or_field=distField, line_side=str(side[a]), line_end_type="FLAT", dissolve_option="NONE", dissolve_field="", method="PLANAR")
              print("did " + str(buffs[b]) + "m front buffer and wrote .shp to /processed/zones")
    except Exception as e:
        print(e)


##
### plot R19 -- LEFT
##try:
##    front="E:/FragEVI/processed/zones/R19_front.shp"
##    print("found front of plot R19")
##except Exception as e:
##    print(e)
##
##
##buffs = [10, 20, 30, 40, 50]
##for b in range(0,5):
##    try:
##        distField = str(buffs[b])+" Meters"
##        if arcpy.Exists("E:/FragEVI/processed/zones/R19_"+str(buffs)+"m.shp"):
##           arcpy.Delete_management("E:/FragEVI/processed/zones/R19_"+str(buffs)+"m.shp")
##        frontBuff = arcpy.Buffer_analysis(in_features=front, out_feature_class="E:/FragEVI/processed/zones/R19_"+str(buffs[b])+"m.shp", buffer_distance_or_field=distField, line_side="LEFT", line_end_type="FLAT", dissolve_option="NONE", dissolve_field="", method="PLANAR")
##        print("did " + str(buffs[b]) + "m front buffer and wrote .shp to /processed/zones")
##    except Exception as e:
##        print(e)
##
##
### plot R20 -- LEFT
##try:
##    front="E:/FragEVI/processed/zones/R20_front.shp"
##    print("found front of plot R20")
##except Exception as e:
##    print(e)
##
##
##buffs = [10, 20, 30, 40, 50]
##for b in range(0,5):
##    try:
##        distField = str(buffs[b])+" Meters"
##        if arcpy.Exists("E:/FragEVI/processed/zones/R20_"+str(buffs)+"m.shp"):
##           arcpy.Delete_management("E:/FragEVI/processed/zones/R20_"+str(buffs)+"m.shp")
##        frontBuff = arcpy.Buffer_analysis(in_features=front, out_feature_class="E:/FragEVI/processed/zones/R20_"+str(buffs[b])+"m.shp", buffer_distance_or_field=distField, line_side="LEFT", line_end_type="FLAT", dissolve_option="NONE", dissolve_field="", method="PLANAR")
##        print("did " + str(buffs[b]) + "m front buffer and wrote .shp to /processed/zones")
##    except Exception as e:
##        print(e)
##
