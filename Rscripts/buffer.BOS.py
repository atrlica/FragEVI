# import system modules 
import arcpy
from arcpy import env

# Set environment settings
env.workspace = "C:\Users\atrlica\Documents\ArcGIS\Default.gdb"

# FID_Boston_Neighborhoods 0:25
# hoods = range(26)
# for N in hoods:
# try:
# Select neighborhood
in_neigh = "FID_Boston_Neighborhoods"
out_feat = "C:\Users\atrlica\Documents\ArcGIS\Default.gdb\pytest_N0"
whereClause = "FID_Boston_Neighborhoods = 0" 
arcpy.Select_analysis(in_neigh, out_feat, whereClause)

# Buffer areas of impact around major roads
can = "bos_can"
buffOut = "c:\Users\atrlica\Documents\ArcGIS\Default.gdb"
distanceField = -10  # for 10 m buffer
sideType = "OUTSIDE_ONLY"
endType = "ROUND"
dissolveType = "NONE"
arcpy.Buffer_analysis(can, buffOut, distanceField, sideType, endType, dissolveType)

# end for loop
