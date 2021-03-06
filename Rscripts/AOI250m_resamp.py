# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "C:/Users/atrlica/Documents/ArcGIS/Default.gdb/"
arcpy.env.overwriteOutput = True
try:
    snapR = arcpy.Raster("E:/FragEVI/processed/EVI/MOD_2010-2012_EVI_NAD83.tif")
    print("found 250m EVI for grid guide")
except Exception as e:
    print(e)

arcpy.env.snapRaster = snapR
csRes = arcpy.GetRasterProperties_management(in_raster="E:/FragEVI/processed/EVI/MOD_2010-2012_EVI_NAD83.tif", property_type="CELLSIZEX", band_index="")
cs = csRes.getOutput(0)
print("cellsizex="+str(cs))

# reproject and resample ISA 250m to the EVI 250m grid (already reprojected to the AOI CRS)
arcpy.Resample_management(in_raster="E:/FragEVI/processed/isa.250m.tif", out_raster="E:/FragEVI/processed/isa.250m.res.tif", cell_size=str(cs), resampling_type="BILINEAR")
print("reprojected/resampled ISA 250m to EVI 250m grid")

### the below won't work for some reason, odd errors in Arc. Did all LULC fraction processing in R directly.
### resample/aggregate the 30m flagged 0/1 LULC rasters to the 250m grid
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.aoi.tif", "E:/FragEVI/processed/boston/bos.aoi30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished aoi 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.ndvi.tif", "E:/FragEVI/processed/boston/bos.ndvi30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished ndvi 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.can.tif", "E:/FragEVI/processed/boston/bos.can30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished can 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.grass.tif", "E:/FragEVI/processed/boston/bos.grass30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished grass 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.barr.tif", "E:/FragEVI/processed/boston/bos.barr30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished barr 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.isa.tif", "E:/FragEVI/processed/boston/bos.isa30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished isa 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.nonimpbarr.tif", "E:/FragEVI/processed/boston/bos.nonimpbarr30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished nonimpbarr 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.vegisa.tif", "E:/FragEVI/processed/boston/bos.vegisa30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished vegisa 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.ed10.tif", "E:/FragEVI/processed/boston/bos.ed1030m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished ed10 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.ed20.tif", "E:/FragEVI/processed/boston/bos.ed2030m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished ed20 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.ed30.tif", "E:/FragEVI/processed/boston/bos.ed3030m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished ed30 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.buff.20only.tif", "E:/FragEVI/processed/boston/bos.buff.20only30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished buff.20only 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.buff.30only.tif", "E:/FragEVI/processed/boston/bos.buff.30only30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished buff.30only 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.buff.Intonly.tif", "E:/FragEVI/processed/boston/bos.buff.Intonly30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished buff.Intonly 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.forest.tif", "E:/FragEVI/processed/boston/bos.forest30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished forest 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.dev.tif", "E:/FragEVI/processed/boston/bos.dev30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished dev 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.hd.res.tif", "E:/FragEVI/processed/boston/bos.hd.res30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished hd.res 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.med.res.tif", "E:/FragEVI/processed/boston/bos.med.res30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished med.res 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.low.res.tif", "E:/FragEVI/processed/boston/bos.low.res30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished low.res 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.lowveg.tif", "E:/FragEVI/processed/boston/bos.lowveg30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished lowveg 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.other.tif", "E:/FragEVI/processed/boston/bos.other30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished other 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.water.tif", "E:/FragEVI/processed/boston/bos.water30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished water 30m")
