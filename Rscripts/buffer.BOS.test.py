# import system modules 
import arcpy
from arcpy import env

# Set environment settings
env.workspace = "C:\Users\atrlica\Documents\ArcGIS\Default.gdb"
env.overwriteOutput = True

# select poly -- individual polygon approach
# in_feat = "E:/FragEVI/processed/bos_can.shp"
# out_feat = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/canChunk"
# whereClause = "FID = " + str(N)
# arcpy.Select_analysis(in_feat, out_feat, whereClause)

# loop select+MRBuffer on all following polys, union into main file after each (this may never work)

buffLength = -20 # select buffer width to use
hoods = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
for N in range(0, 24):
    try:
        # select neighborhood
        in_neigh = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/Bos_canNeighNAD83_int"
        out_feat = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/N_select"
        whereClause = "FID_Boston_Neighborhoods_NAD83 = " + str(hoods[N])
        arcpy.Select_analysis(in_neigh, out_feat, whereClause)
        print("working on neighborhood " + str(hoods[N]))
        # peform buffer
        in_feat = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/N_select"
        buffOut = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/Bos_buff" + str(buffLength*-1) + "_N" + str(hoods[N])
        distanceField = buffLength  
        sideType = "OUTSIDE_ONLY"
        endType = "ROUND"
        dissolveType = "NONE"
        dissolveField = ""
        method = "PLANAR"
        arcpy.Buffer_analysis(in_feat, buffOut, distanceField, sideType, endType, dissolveType, dissolveField, method)
        print("finished " + str(buffLength*-1) + "m buffer in neighborhood " + str(hoods[N]))
    except Exception as e:
        print e
        
# end for loop
