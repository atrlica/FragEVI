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

# reproject and resample ISA 1m grid to boston canopy 1 m grid
arcpy.ProjectRaster_management(in_raster="E:/FragEVI/processed/bos_isa_1m.tif", out_raster="E:/FragEVI/processed/isa_cangrid.tif", out_coor_system="PROJCS['NAD_1983_UTM_Zone_19N',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',500000.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-69.0],PARAMETER['Scale_Factor',0.9996],PARAMETER['Latitude_Of_Origin',0.0],UNIT['Meter',1.0]]", resampling_type="NEAREST", cell_size=1, geographic_transform="", Registration_Point="", in_coor_system="PROJCS['Lambert_Conformal_Conic_2SP',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['false_easting',200000.0],PARAMETER['false_northing',750000.0],PARAMETER['central_meridian',-71.5],PARAMETER['standard_parallel_1',41.71666666666667],PARAMETER['standard_parallel_2',42.68333333333333],PARAMETER['latitude_of_origin',41.0],UNIT['Meter',1.0]]")
