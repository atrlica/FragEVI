# import system modules 
import arcpy, arcinfo
from arcpy import env
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

# Set environment settings
arcpy.env.workspace = "E:/FragEVI/FragEVIworking.gdb/"
arcpy.env.overwriteOutput = True
try:
    snapR = arcpy.Raster("E:/FragEVI/processed/EVI/030005-6_2010-2012_EVI_NAD83.tif")
    print("found 30m EVI for grid guide")
except Exception as e:
    print(e)

arcpy.env.snapRaster = snapR
csRes = arcpy.GetRasterProperties_management(in_raster="E:/FragEVI/data/dataverse_files/bostoncanopy_1m.tif", property_type="CELLSIZEX", band_index="")
cs = csRes.getOutput(0)
print("cellsizex="+str(cs))

# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.aoi.tif", "E:/FragEVI/processed/boston/bos.aoi30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished aoi 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.ndvi.tif", "E:/FragEVI/processed/boston/bos.ndvi30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished ndvi 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.can.tif", "E:/FragEVI/processed/boston/bos.can30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished can 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.can.redux.tif", "E:/FragEVI/processed/boston/bos.can.redux30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished can.redux 30m")
arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.grass.tif", "E:/FragEVI/processed/boston/bos.grass30m.tif", "30", "MEAN", "EXPAND", "DATA")
print("finished grass 30m")
arcpy.gp.Aggregate_sa("E:/BosBiog/processed/bos.nogolfturf1m.tif", "E:/BosBiog/processed/bos.nogolfturf30m.tif", "30", "MEAN", "EXPAND", "DATA")
print("finished nogolfturf 30m, wrote to E:/BosBiog")
arcpy.gp.Aggregate_sa("E:/BosBiog/processed/bos.golfturf1m.tif", "E:/BosBiog/processed/bos.golfturf30m.tif", "30", "MEAN", "EXPAND", "DATA")
print("finished golfturf 30m, wrote to E:/BosBiog")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.barr.tif", "E:/FragEVI/processed/boston/bos.barr30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished barr 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.isa.RR2.tif", "E:/FragEVI/processed/boston/bos.isa30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished isa 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.nonimpbarr.tif", "E:/FragEVI/processed/boston/bos.nonimpbarr30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished nonimpbarr 30m")
arcpy.gp.Aggregate_sa("E:/BosBiog/processed/bos.barr1m.tif", "E:/BosBiog/processed/bos.barr30m.tif", "30", "MEAN", "EXPAND", "DATA")
print("finished pervious barren 30m, wrote to E:BosBiog")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.vegisa.tif", "E:/FragEVI/processed/boston/bos.vegisa30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished vegisa 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.ed10m.redux.tif", "E:/FragEVI/processed/boston/bos.ed10m.redux30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished ed10.redux 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.ed20m.redux.tif", "E:/FragEVI/processed/boston/bos.ed20m.redux30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished ed20.redux 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.ed30m.redux.tif", "E:/FragEVI/processed/boston/bos.ed30m.redux30m.tif", "30", "MEAN", "EXPAND", "DATA")
##print("finished ed30.redux 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.buff.20only.tif", "E:/FragEVI/processed/boston/bos.buff.20only30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished buff.20only 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.buff.30only.tif", "E:/FragEVI/processed/boston/bos.buff.30only30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished buff.30only 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.buff.Intonly.tif", "E:/FragEVI/processed/boston/bos.buff.Intonly30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished buff.Intonly 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.forest.tif", "E:/FragEVI/processed/boston/bos.forest30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished forest 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.dev.tif", "E:/FragEVI/processed/boston/bos.dev30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished dev 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.hd.res.tif", "E:/FragEVI/processed/boston/bos.hd.res30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished hd.res 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.med.res.tif", "E:/FragEVI/processed/boston/bos.med.res30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished med.res 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.low.res.tif", "E:/FragEVI/processed/boston/bos.low.res30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished low.res 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.lowveg.tif", "E:/FragEVI/processed/boston/bos.lowveg30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished lowveg 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.other.tif", "E:/FragEVI/processed/boston/bos.other30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished other 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.water.tif", "E:/FragEVI/processed/boston/bos.water30m.tif", "30", "MEAN", "EXPAND", "DATA")
# print("finished water 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/data/dataverse_files/bostonbiomass_1m.tif", "E:/FragEVI/processed/boston/bos.biom30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished biomass 30m")
arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.ed10m.biom.redux.tif", "E:/FragEVI/processed/boston/bos.biomass.ed10only30m.tif", "30", "SUM", "EXPAND", "DATA")
print("finished biomass.edgeonly 30m")
# arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.biomass.forestonly.tif", "E:/FragEVI/processed/boston/bos.biomass.forestonly30m.tif", "30", "SUM", "EXPAND", "DATA")
# print("finished biomass.forestonly 30m")
##arcpy.gp.Aggregate_sa("E:/FragEVI/processed/boston/bos.biomass.ed10only_allLULC.tif", "E:/FragEVI/processed/boston/bos.biomass.ed10only_allLULC30m.tif", "30", "SUM", "EXPAND", "DATA")
##print("finished biomass.ed10only_allLULC 30m")



