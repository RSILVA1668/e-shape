
"""
This scrips contains functions needed in the e-shape pilot monitoring fishing activities.
To use those functions: from e_shape_functions import *
This will import functions and databases as variables
"""
######### Import libraries needed
import os , gdal , osr , datetime
os.environ ['PROJ_LIB' ] = r'C:\Users\macl\AppData\Local\Continuum\anaconda3\pkgs\proj4-5.2.0-h6538335_1006\Library\share'
from matplotlib import colors as mcolors
import random
import numpy as np
import pandas as pd
import warnings
import shapely
import sys , time
import geopandas as gpd
from netCDF4 import Dataset as dt
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import gdal
import math
from math import radians , cos , sin , asin , sqrt
from shapely.geometry import *
from matplotlib.widgets import RectangleSelector
from pandas.plotting import register_matplotlib_converters
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)



######### Import of databases (Ports of PT, PT EEZ shapefile, sales from DGRM...)
main_path = r'G:\Mon Drive\E-SHAPE_to_share\0_Data'
df_allmmsi = pd.read_csv ( main_path + r'\11_vessel_data_match\all_mmsi_and_parameters_0-2684_gear_ALL_EU.csv' )
vessel_database_path = main_path + r'\11_vessel_data_match\combineOnlyPortugal_all_mmsi_and_parameters_DGRM_vesselCode.csv'
vessel_database = pd.read_csv ( vessel_database_path , sep = ';' )
eez_shp_path = r'C:\Users\macl\Documents\E-shape\1_EEZ_PT_all\eez.shp'
eez_shp = gpd.read_file ( eez_shp_path )

### 3 first digits of MMSI Id number correspond to a country
df_flags = pd.read_excel ( main_path + r'\9_vessels_registeries\MMSI_id.xlsx' )

### Ports of PT as Polygons with their name. This allow to spot when emissions are inside or outside the ports
shp_file = pd.read_csv(main_path + r'\15_Portugal-Ports\ports_polygon.csv')
polygons = [ shapely.wkt.loads ( polygon ) for polygon in shp_file [ 'geometry' ] ]
ports = gpd.GeoDataFrame ( shp_file.drop ( columns = [ 'geometry' ] ) , geometry = polygons )
a = [ ]
for i in range ( 0 , len ( ports) ) :
    a.append ( ports.iloc [ i , -1 ] )
p = MultiPolygon ( a )

# This part of the code is supposed to estimate the depth for each emissions but it is highly uncertain
# filestr1 = main_path.split ( '.' ) [ 0 ] + r'\17_Emodnet\13_Bathymetry'
#
# tiles = [ '\F1_2018' , '\F2_2018' , '\F3_2018' , '\G1_2018' , '\G2_2018' , '\G3_2018' ]
#
# filestr = filestr1.split ( '.' ) [ 0 ] + '%s.dtm' % (tiles [ 0 ])
# fh = dt ( filestr , mode = 'r' )  ## Read the file
# xF1 = np.arange ( fh.Longitude_TL , fh.Longitude_BR , fh.Element_x_size )
# yF1 = np.arange ( fh.Latitude_BL , fh.Latitude_TL , fh.Element_x_size )
# depthF1 = fh.variables [ 'DEPTH_SMOOTH' ] [ : ]
# depthF1 = pd.DataFrame ( depthF1 )
#
# filestr = filestr1.split ( '.' ) [ 0 ] + '%s.dtm' % (tiles [ 1 ])
# fh = dt ( filestr , mode = 'r' )  ## Read the file
# xF2 = np.arange ( fh.Longitude_TL , fh.Longitude_BR , fh.Element_x_size )
# yF2 = np.arange ( fh.Latitude_BL , fh.Latitude_TL , fh.Element_x_size )
# depthF2 = fh.variables [ 'DEPTH_SMOOTH' ] [ : ]
# depthF2 = pd.DataFrame ( depthF2 )
#
# filestr = filestr1.split ( '.' ) [ 0 ] + '%s.dtm' % (tiles [ 2 ])
# fh = dt ( filestr , mode = 'r' )  ## Read the file
# xF3 = np.arange ( fh.Longitude_TL , fh.Longitude_BR , fh.Element_x_size )
# yF3 = np.arange ( fh.Latitude_BL , fh.Latitude_TL , fh.Element_x_size )
# depthF3 = fh.variables [ 'DEPTH_SMOOTH' ] [ : ]
# depthF3 = pd.DataFrame ( depthF3 )
#
# filestr = filestr1.split ( '.' ) [ 0 ] + '%s.dtm' % (tiles [ 3 ])
# fh = dt ( filestr , mode = 'r' )  ## Read the file
# xG1 = np.arange ( fh.Longitude_TL , fh.Longitude_BR , fh.Element_x_size )
# yG1 = np.arange ( fh.Latitude_BL , fh.Latitude_TL , fh.Element_x_size )
# depthG1 = fh.variables [ 'DEPTH_SMOOTH' ] [ : ]
# depthG1 = pd.DataFrame ( depthG1 )
#
# filestr = filestr1.split ( '.' ) [ 0 ] + '%s.dtm' % (tiles [ 4 ])
# fh = dt ( filestr , mode = 'r' )  ## Read the file
# xG2 = np.arange ( fh.Longitude_TL , fh.Longitude_BR , fh.Element_x_size )
# yG2 = np.arange ( fh.Latitude_BL , fh.Latitude_TL , fh.Element_x_size )
# depthG2 = fh.variables [ 'DEPTH_SMOOTH' ] [ : ]
# depthG2 = pd.DataFrame ( depthG2 )
#
# filestr = filestr1.split ( '.' ) [ 0 ] + '%s.dtm' % (tiles [ 5 ])
# fh = dt ( filestr , mode = 'r' )  ## Read the file
# xG3 = np.arange ( fh.Longitude_TL , fh.Longitude_BR , fh.Element_x_size )
# yG3 = np.arange ( fh.Latitude_BL , fh.Latitude_TL , fh.Element_x_size )
# depthG3 = fh.variables [ 'DEPTH_SMOOTH' ] [ : ]
# depthG3 = pd.DataFrame ( depthG3 )
# Functions to find the depth
# def find_nearest ( yy , xx , value1 , value2 ) :
#     yy = np.asarray ( yy )
#     xx = np.asarray ( xx )
#     idx1 = (np.abs ( yy - value1 )).argmin ( )
#     idx2 = (np.abs ( xx - value2 )).argmin ( )
#     return idx1 , idx2
#

# def find_depth ( lat , lon ) :
#     global depth
#     if np.abs ( xF1 - lon ).argmin ( ) > 0 :
#         if np.abs ( yF1 - lat ).argmin ( ) > 0 :
#             idx1 , idx2 = find_nearest ( yF1 , xF1 , lat , lon )
#             depth.append ( depthF1.iloc [ idx1 , idx2 ] )
#         elif np.abs ( yG1 - lat ).argmin ( ) > 0 :
#             idx1 , idx2 = find_nearest ( yG1 , xF1 , lat , lon )
#             depth.append ( depthG2.iloc [ idx1 , idx2 ] )
#         else :
#             depth.append ( None )
#     elif np.abs ( xF2 - lon ).argmin ( ) > 0 :
#         if np.abs ( yF2 - lat ).argmin ( ) > 0 :
#             idx1 , idx2 = find_nearest ( yF2 , xF2 , lat , lon )
#             depth.append ( depthF2.iloc [ idx1 , idx2 ] )
#         elif np.abs ( yG2 - lat ).argmin ( ) > 0 :
#             idx1 , idx2 = find_nearest ( yG2 , xF2 , lat , lon )
#             depth.append ( depthG3.iloc [ idx1 , idx2 ] )
#         else :
#             depth.append ( None )
#     else :
#         depth.append ( None )
#     return depth


def landings_geo(year):
# This is a function that can sort the sales from landing declaration depending the year
# The output of it is 3 dataframes for Azores, Madeiras and PT Continental to do stats on the sales
    global fishin_pt_cont,fishin_az,fishin_mad
    dfport = pd.read_csv ( r'C:\Users\macl\Documents\E-shape\11_MC_Python_Codes\ports_PT_GPS.csv' )
    dfport = dfport.iloc [ : , 1 : ]
    gdfport = gpd.GeoDataFrame ( dfport , geometry = gpd.points_from_xy ( dfport.LON , dfport.LAT ) )
    del gdfport [ 'LON' ]
    del gdfport [ 'LAT' ]
    fishin = pd.read_excel (
        main_path + r'\14_IPMA_LandingDeclarations\DADOS_ENVIADOS_DEIMOS\DESEMBARQUES_%s.xlsx' % (str ( year )) )
    madeiras = Polygon ( [ (-17.83161 , 33.355298498597385) , (-15.810131068255945 , 33.355298498597385) ,
                           (-15.810131068255945 , 32.26575824509192) , (-17.831615443255945 , 32.26575824509192) ,
                           (-17.831615443255945 , 33.355298498597385) ] )
    PT_continental = Polygon ( [ (-10.5529112 , 42.491914) , (-6.4220518 , 42.491914) , (-6.4220518 , 36.5312314) ,
                                 (-10.5529112 , 36.5312314) , (-10.5529112 , 42.491914) ] )
    azores = Polygon ( [ (-29.84233865556249 , 39.95774470243776) , (-22.59136209306249 , 39.95774470243776) ,
                         (-22.59136209306249 , 36.75383990140275) , (-29.84233865556249 , 36.75383990140275) ,
                         (-29.84233865556249 , 39.95774470243776) ] )
    port_az = [ ]
    port_mad = [ ]
    port_pt_cont = [ ]
    for i in range ( len ( gdfport ) ) :
        if gdfport.iloc [ i , 1 ].within ( azores ) :
            port_az.append ( i )
        if gdfport.iloc [ i , 1 ].within ( madeiras ) :
            port_mad.append ( i )
        if gdfport.iloc [ i , 1 ].within ( PT_continental ) :
            port_pt_cont.append ( i )
    ports_az = dfport.loc [ port_az ]
    ports_mad = dfport.loc [ port_mad ]
    ports_pt_cont = dfport.loc [ port_pt_cont ]
    liste_az = [ ]
    liste_mad = [ ]
    liste_pt_cont = [ ]
    for i in range ( len ( fishin ) ) :
        if fishin.loc [ i , 'PORTO_DESEMBARQUE' ] in list(ports_az['0']) :
            liste_az.append ( i )
        if fishin.loc [ i , 'PORTO_DESEMBARQUE' ] in list(ports_mad['0']) :
            liste_mad.append ( i )
        if fishin.loc [ i , 'PORTO_DESEMBARQUE' ] in list(ports_pt_cont['0']) :
            liste_pt_cont.append ( i )
    fishin_az = fishin.loc [ liste_az ]
    fishin_mad = fishin.loc [ liste_mad ]
    fishin_pt_cont = fishin.loc [ liste_pt_cont ]


def calculate_initial_compass_bearing ( pointA , pointB ) :
    """
    Calculates the bearing between two points.
    The formulae used is the following:
        θ = atan2(sin(Δlong).cos(lat2),
                  cos(lat1).sin(lat2) − sin(lat1).cos(lat2).cos(Δlong))
                  """
    lat1 = math.radians ( pointA.y )
    lat2 = math.radians ( pointB.y )
    diffLong = math.radians ( pointB.x - pointA.x )
    x = math.sin ( diffLong ) * math.cos ( lat2 )
    y = math.cos ( lat1 ) * math.sin ( lat2 ) - (math.sin ( lat1 )
                                                 * math.cos ( lat2 ) * math.cos ( diffLong ))
    initial_bearing = math.atan2 ( x , y )
    # Now we have the initial bearing but math.atan2 return values
    # from -180° to + 180° which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = math.degrees ( initial_bearing )
    compass_bearing = (initial_bearing + 360) % 360
    return compass_bearing


def distance ( lat1 , lat2 , lon1 , lon2 ) :
    # The math module contains a function named
    # radians which converts from degrees to radians.
    lon1 = radians ( lon1 )
    lon2 = radians ( lon2 )
    lat1 = radians ( lat1 )
    lat2 = radians ( lat2 )
    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin ( dlat / 2 ) ** 2 + cos ( lat1 ) * cos ( lat2 ) * sin ( dlon / 2 ) ** 2
    c = 2 * asin ( sqrt ( a ) )
    # Radius of earth in kilometers. Use 3956 for miles
    r = 6371
    # calculate the result
    return (c * r)


def distance_gps ( point1 , point2 ) :
    # The math module contains a function named
    # radians which converts from degrees to radians.

    lon1 = radians ( point1.x )
    lon2 = radians ( point2.x )
    lat1 = radians ( point1.y )
    lat2 = radians ( point2.y )
    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin ( dlat / 2 ) ** 2 + cos ( lat1 ) * cos ( lat2 ) * sin ( dlon / 2 ) ** 2
    c = 2 * asin ( sqrt ( a ) )
    # Radius of earth in kilometers. Use 3956 for miles
    r = 6371
    # calculate the result
    return (c * r)


def plotpath1 ( paths , dataframe , int , relvel_threshold=1.5 , timelaps_threshold=3600 ) :
# This function plots a path with a green point at the first emission and a red cross at the last
# there is a blue star for emission when there is more than 1h between 2 emissions
# there are red points when the emission has a relative velocity under 2.5 knots
# The int number is the index of the paths in the "paths" dataframe, "dataframe" is used to search the
# first and last emissions
    global current_ax
    current_ax = eez_shp.plot ( )
    beg_col1 = paths.columns.get_loc ( "Begining Date" )
    end_col1 = paths.columns.get_loc ( "Ending Date" )
    for i in range ( 0 , len ( ports ) ) :
        x , y = ports.iloc [ i , -1 ].exterior.xy
        plt.plot ( x , y , 'r' )
    geom_col = paths.columns.get_loc ( "geometry" )
    x , y = paths.iloc [ int , geom_col ].xy
    plt.plot ( x , y , 'k.-' )
    test = dataframe [
        (dataframe [ 'DATE' ] <= paths.iloc [ int , end_col1 ]) & (
                dataframe [ 'DATE' ] >= paths.iloc [ int , beg_col1 ]) & (
                dataframe [ 'MMSI' ] == paths.iloc [ int , 0 ]) & (dataframe [ 'RelVel' ] <= relvel_threshold) ]
    geom_col1 = dataframe.columns.get_loc ( "geometry" )

    for i in range ( len ( test ) ) :
        x , y = test.iloc [ i , geom_col1 ].xy
        plt.plot ( x , y , 'ro' , LineWidth = 0.3 )
    test = dataframe [
        (dataframe [ 'DATE' ] <= paths.iloc [ int , end_col1 ]) & (
                dataframe [ 'DATE' ] >= paths.iloc [ int , beg_col1 ]) & (
                dataframe [ 'MMSI' ] == paths.iloc [ int , 0 ]) & (dataframe [ 'TimeLaps' ] >= timelaps_threshold) ]
    geom_col1 = dataframe.columns.get_loc ( "geometry" )
    for i in range ( len ( test ) ) :
        x , y = test.iloc [ i , geom_col1 ].xy
        plt.plot ( x , y , 'c*' , LineWidth = 0.3 )
    a = Point ( paths.iloc [ int , geom_col ].coords [ 0 ] )
    x1 , y1 = a.xy
    plt.plot ( x1 , y1 , 'go' )
    b = Point ( paths.iloc [ int , geom_col ].coords [ -1 ] )
    x1 , y1 = b.xy
    plt.plot ( x1 , y1 , 'r+' )
    print ( 'MMSI num : ' + str ( paths.iloc [ int , 0 ] ) + ' \nFrom ' + str ( paths.iloc [ int , beg_col1 ] )
            + ' to ' + str ( paths.iloc [ int , end_col1 ] ) + ' \nGear ' +
            str ( df_allmmsi [ df_allmmsi [ 'MMSI' ] == paths.iloc [ int , 0 ] ].iloc [ 0 , 10 ] ) )


# Functions that lead to plotting a path, and in another figure, the relative velocity profile and compass
# def onclick ( event ) :
#     global ix , iy , coords , current_ax
#     ix , iy = int ( event.xdata ) , event.ydata
#     print ( ix , iy )
#     coords = [ ix , iy ]
#     x1 , y1 = test.iloc [ coords [ 0 ] , -6 ].xy
#     current_ax.plot ( x1 , y1 , 'dc' )
#     return coords
#
#
# def clickndraw ( test , current_ax ) :
#     test = res ( test )
#     fig = plt.figure ( )
#     plt.plot ( test [ 'RelVel' ] , 'b' , label = 'Compass' , )
#     fig.canvas.mpl_connect ( 'button_press_event' , onclick )
#
#
# def double_plot ( dfsure , i ) :
#     global test
#     test = look_into ( dataframe , dfsure , i )
#     plotpath1 ( dfsure , dataframe , i )
#     plot_relvel_compass ( test )
#     clickndraw ( test , current_ax )

### I created from the data of Marine Traffic, csv files of areas 1 & 2 and cut them considering years (2012 - 2018)
### This function import this csv file turning it into a geodataframe


def import_dataframe ( year ) :
# This function takes the subsets that I made from areas 1 to 4 and add columns for the relative velocity, the compass
# It also remodels the dataframe removing the wrong data.
    global gdf , dataframe , fishin

    fishin = pd.read_excel (
        main_path + r'\14_IPMA_LandingDeclarations\DADOS_ENVIADOS_DEIMOS\DESEMBARQUES_%s.xlsx' % (str ( year )) )
    main_path_az = main_path.split ( '.' ) [ 0 ] + r'\5_AIS_IPMA\Subset_Area_1_2\areas%s.csv' % (str ( year ))
    # main_path_mad = main_path.split ( '.' ) [ 0 ] + r'\5_AIS_IPMA\Subset_Area_3\AREA3_%s.csv' % (str ( year ))
    # main_path_cont = main_path.split ( '.' ) [ 0 ] + r'\5_AIS_IPMA\Subset_Area_4\AREA4_%s.csv' % (str ( year ))
    df = pd.read_csv ( main_path_az )
    df.rename ( columns = lambda x : x.replace ( 'TIMESTAMP(UTC)' , 'DATE' ) ,
                inplace = True )  ## Rename the column into DATE
    df [ 'DATE' ] = pd.to_datetime ( df.DATE )  ## Translate every date into a datetime type
    df = df.sort_values ( by = [ 'MMSI' , 'DATE' ] )  ## Sort by MMSI number and date
    df = df.iloc [ : , 1 : ]  ## The first column is the indexes of the previous dataframe, no use.
    ## This DataFrame has "," instead of "." in LON/LAT columns. We have to change that and change the columns
    ## from "str" to float type
    if type ( df.loc [ 0 , 'LON' ] ) == str :
        df [ "LON" ] = [ x.replace ( "," , "." ) for x in df [ "LON" ] ]
        df [ "LON" ] = df [ "LON" ].astype ( float )
        df [ "LAT" ] = [ x.replace ( "," , "." ) for x in df [ "LAT" ] ]
        df [ "LAT" ] = df [ "LAT" ].astype ( float )
    ## Creating geodataframe
    gdf = gpd.GeoDataFrame ( df , geometry = gpd.points_from_xy ( df.LON , df.LAT ) )
    date_col = gdf.columns.get_loc ( "DATE" )
    geom_col = gdf.columns.get_loc ( "geometry" )
    ## Removing the same date twice
    diff = [ ]
    diff.append ( True )
    for i in range ( 1 , len ( gdf ) ) :
        diff.append ( gdf.iloc [ i , date_col ] != gdf.iloc [ i - 1 , date_col ] )
    gdf = gdf [ (diff) ]
    ## Creating a Relative Velocity Column
    gdf [ 'RelVel' ] = 0
    relvel = [ ]
    # for i in range ( len ( gdf ) ) :
    #     find_depth ( gdf.iloc [ i , 5 ] , gdf.iloc [ i , 4 ] )
    TimeLaps = [ ]
    Dist = [ ]
    compass = [ ]
    relvel.append ( 0 )
    compass.append ( 0 )
    TimeLaps.append ( 0 )
    Dist.append ( 0 )
    for i in range ( 0 , gdf.shape [ 0 ] - 1 ) :
        bearing = calculate_initial_compass_bearing ( gdf.iloc [ i , geom_col ] , gdf.iloc [ i + 1 , geom_col ] )
        dist = distance ( gdf.iloc [ i , geom_col ].y , gdf.iloc [ i + 1 , geom_col ].y , gdf.iloc [ i , geom_col ].x ,
                          gdf.iloc [ i + 1 , geom_col ].x )

        duration = gdf.iloc [ i + 1 , date_col ] - gdf.iloc [ i , date_col ]
        compass.append ( bearing )
        duration_in_h = (duration.total_seconds ( ) / 3600)
        vitrel = dist / duration_in_h
        relvel.append ( vitrel )
        TimeLaps.append ( duration.total_seconds ( ) )
        Dist.append ( dist )
    # gdf [ 'Depth' ] = depth
    gdf [ 'RelVel' ] = relvel
    gdf [ 'TimeLaps' ] = TimeLaps
    gdf [ 'Compass' ] = compass
    gdf [ 'Dist' ] = Dist
    gdf [ 'RelVel' ] *= 0.539957  ## into knots!!!!
    gdf [ 'SPEED(KNOTS x10)' ] /= 10
    gdf = gdf [ (gdf [ 'MMSI' ] >= 100000000) & (gdf [ 'MMSI' ] <= 800000000) ]
    gdf = gdf [ (gdf [ 'MMSI' ] != 111111111) & (gdf [ 'MMSI' ] != 222222222) & (gdf [ 'MMSI' ] != 333333333) & (
            gdf [ 'MMSI' ] != 444444444) ]
    gdf = gdf [ (gdf [ 'MMSI' ] != 555555555) & (gdf [ 'MMSI' ] != 666666666) & (gdf [ 'MMSI' ] != 777777777) & (
            gdf [ 'MMSI' ] != 888888888) & (gdf [ 'MMSI' ] != 999999999) ]
    myList = list ( gdf.columns )
    myList [ 3 ] = 'SPEED' ;
    gdf.columns = myList

    ### Some of the paths created later are wrong because of positions misplaced, this part is there to remove those
    ### wrong emissions that will make the paths whole again

    gdf = gdf.reset_index ( )
    gdf = gdf.iloc [ : , 1 : ]
    gdf30 = gdf [ gdf [ 'RelVel' ] > 30 ]
    gdf30 = gdf30.reset_index ( )
    todrop = [ ]
    for i in gdf30 [ 'index' ] :
        test = gdf.iloc [ i - 1 :i + 2 , : ]
        dist = distance ( test.iloc [ 0 , 9 ].y , test.iloc [ -1 , 9 ].y , test.iloc [ 0 , 9 ].x ,
                          test.iloc [ -1 , 9 ].x )
        duration = test.iloc [ -1 , 8 ] - test.iloc [ 0 , 8 ]
        duration_in_h = (duration.total_seconds ( ) / 3600)
        vitrel = (dist / duration_in_h) * 0.539957
        if vitrel < 30 :
            todrop.append ( i )
        else :
            test = gdf.iloc [ i - 2 :i + 1 , : ]
            dist = distance ( test.iloc [ 0 , 9 ].y , test.iloc [ -1 , 9 ].y , test.iloc [ 0 , 9 ].x ,
                              test.iloc [ -1 , 9 ].x )
            duration = test.iloc [ -1 , 8 ] - test.iloc [ 0 , 8 ]
            duration_in_h = (duration.total_seconds ( ) / 3600)
            vitrel = (dist / duration_in_h) * 0.539957
            if vitrel < 30 :
                todrop.append ( i - 1 )
    gdf.drop ( todrop , axis = 0 , inplace = True )
    relvel = [ ]
    relvel.append ( 0 )
    for i in range ( gdf.shape [ 0 ] - 1 ) :
        dist = distance ( gdf.iloc [ i , 9 ].y , gdf.iloc [ i + 1 , 9 ].y , gdf.iloc [ i , 9 ].x ,
                          gdf.iloc [ i + 1 , 9 ].x )
        duration = gdf.iloc [ i + 1 , 8 ] - gdf.iloc [ i , 8 ]
        duration_in_h = (duration.total_seconds ( ) / 3600)
        vitrel = (dist / duration_in_h) * 0.539957
        relvel.append ( vitrel )
    gdf [ 'RelVel' ] = relvel
    dataframe = pd.DataFrame ( gdf )
    return gdf


def make_paths ( gdf , hour_delay=24 , relvel=30 ) :
# We have a dataframe and we want to create paths taking into account a delay between two emissions and relative velocity
# Thresholds have been selected at 24 hours for the maximum time delay accepted between two emissions
# and 30 Knots for the maximum relative velocity (which is a large speed for a boat but considering the time delay
# and the overall data, it prevents from cutting trips in the middle.
    global paths
    date_col = gdf.columns.get_loc ( "DATE" )
    geom_col = gdf.columns.get_loc ( "geometry" )
    dist_col = gdf.columns.get_loc ( "Dist" )
    TimeLaps_col = gdf.columns.get_loc ( "TimeLaps" )
    RelVel_col = gdf.columns.get_loc ( "RelVel" )
    mmsi_col = gdf.columns.get_loc ( "MMSI" )
    j = 0
    lines = [ ]
    aa = 0
    paths = gpd.GeoDataFrame (
        { 'MMSI' : [ ] , 'Begining Date' : [ ] , 'Ending Date' : [ ] , 'Trip_Length' : [ ] , 'geometry' : [ ] } )
    # Start of the file
    for i in range ( 0 , len ( gdf ) - 1 ) :
        if gdf.iloc [ i , mmsi_col ] == gdf.iloc [ i + 1 , mmsi_col ] :  ## If MMSI number stays the same
            if gdf.iloc [ i , geom_col ].within ( p ) :  # if the GPS point inside the port
                if gdf.iloc [ i + 1 , geom_col ].within ( p ) :  # If the boat stays in a port between 2 emissions
                    i += 1  # the trip has not started
                elif gdf.iloc [ i , geom_col ].within ( p ) and not gdf.iloc [ i + 1 , geom_col ].within ( p ) and \
                        gdf.iloc [ i + 1 , TimeLaps_col ] < hour_delay * 3600 :  # Boat is leaving the port
                    if len ( lines ) == 0 :
                        if gdf.iloc [ i + 1 , RelVel_col ] < relvel and not gdf.iloc [ i + 1 , geom_col ].within (
                                p ) and i + 1 <= len ( gdf ) - 2 and gdf.iloc [
                            i + 1 , TimeLaps_col ] < hour_delay * 3600 :  ##
                            lines.append ( gdf.iloc [ i , geom_col ] )
                            paths.loc [ j , 'Begining Date' ] = gdf.iloc [ i , date_col ]
                        else :
                            i += 1
                    else :
                        # While relative velocity is correct and boat doesn't go back into a port. The timelaps between two emissions has to be small
                        if gdf.iloc [ i + 1 , RelVel_col ] < relvel and not gdf.iloc [ i + 1 , geom_col ].within (
                                p ) and i + 1 <= len ( gdf ) - 2 and gdf.iloc [
                            i + 1 , TimeLaps_col ] < hour_delay * 3600 :
                            lines.append ( gdf.iloc [ i , geom_col ] )
                            aa = aa + gdf.iloc [ i , dist_col ]
                        else :  # Once any of those condition is not met
                            paths.loc [ j , 'MMSI' ] = gdf.iloc [ i , mmsi_col ]
                            aa = aa + gdf.iloc [ i , dist_col ]
                            lines.append ( gdf.iloc [ i , geom_col ] )
                            paths.loc [ j , 'Ending Date' ] = gdf.iloc [ i , date_col ]  # Ending date
                            if gdf.iloc [ i + 1 , RelVel_col ] < relvel and gdf.iloc [ i + 1 , geom_col ].within (
                                    p ) and i + 1 <= len ( gdf ) - 2 and gdf.iloc [
                                i + 1 , TimeLaps_col ] < hour_delay * 3600 :
                                aa = aa + gdf.iloc [ i + 1 , dist_col ]
                                lines.append ( gdf.iloc [ i + 1 , geom_col ] )
                                paths.loc [ j , 'Ending Date' ] = gdf.iloc [ i + 1 , date_col ]
                            if len ( lines ) > 1 :
                                a = LineString ( list ( lines ) )  # Linestring
                                paths.loc [ j , 'Trip_Length' ] = aa
                                paths.loc [ j , 'geometry' ] = a
                            else :
                                paths.loc [ j , 'geometry' ] = gdf.iloc [ i , geom_col ]
                                paths.loc [ j , 'Trip_Length' ] = aa
                            lines [ : ] = [ ]  # Empty the list
                            aa = 0
                            j += 1  ## new line for the next vessels paths
            else :
                ### If the GPS point is not in a port, it is the first of the new path and the same creation of path
                if len ( lines ) == 0 :
                    if gdf.iloc [ i + 1 , RelVel_col ] < relvel and not gdf.iloc [ i + 1 , geom_col ].within (
                            p ) and i + 1 <= len ( gdf ) - 2 and gdf.iloc [ i + 1 , TimeLaps_col ] < hour_delay * 3600 :
                        lines.append ( gdf.iloc [ i , geom_col ] )
                        paths.loc [ j , 'Begining Date' ] = gdf.iloc [ i , date_col ]
                    else :
                        i += 1
                else :
                    if gdf.iloc [ i + 1 , RelVel_col ] < relvel and not gdf.iloc [ i + 1 , geom_col ].within (
                            p ) and i + 1 <= len ( gdf ) - 2 and gdf.iloc [ i + 1 , TimeLaps_col ] < hour_delay * 3600 :
                        lines.append ( gdf.iloc [ i , geom_col ] )
                        aa = aa + gdf.iloc [ i , dist_col ]
                    else :
                        lines.append ( gdf.iloc [ i , geom_col ] )
                        aa = aa + gdf.iloc [ i , dist_col ]
                        paths.loc [ j , 'MMSI' ] = gdf.iloc [ i , mmsi_col ]
                        paths.loc [ j , 'Ending Date' ] = gdf.iloc [ i , date_col ]
                        if gdf.iloc [ i + 1 , RelVel_col ] < relvel and gdf.iloc [ i + 1 , geom_col ].within (
                                p ) and i + 1 <= len ( gdf ) - 2 and gdf.iloc [
                            i + 1 , TimeLaps_col ] < hour_delay * 3600 :
                            aa = aa + gdf.iloc [ i + 1 , dist_col ]
                            lines.append ( gdf.iloc [ i + 1 , geom_col ] )
                            paths.loc [ j , 'Ending Date' ] = gdf.iloc [ i + 1 , date_col ]
                        if len ( lines ) > 1 :
                            a = LineString ( list ( lines ) )
                            paths.loc [ j , 'Trip_Length' ] = aa
                            paths.loc [ j , 'geometry' ] = a
                        else :
                            paths.loc [ j , 'geometry' ] = gdf.iloc [ i , geom_col ]
                            paths.loc [ j , 'Trip_Length' ] = aa
                        lines [ : ] = [ ]
                        aa = 0
                        j += 1
        else :
            paths.loc [ j , 'MMSI' ] = gdf.iloc [ i , mmsi_col ]
            paths.loc [ j , 'Ending Date' ] = gdf.iloc [ i , date_col ]
            if len ( lines ) > 1 :
                a = LineString ( list ( lines ) )
                paths.loc [ j , 'Trip_Length' ] = aa
                paths.loc [ j , 'geometry' ] = a
            else :
                paths.loc [ j , 'Begining Date' ] = gdf.iloc [ i , date_col ]
                paths.loc [ j , 'Trip_Length' ] = aa
                paths.loc [ j , 'geometry' ] = gdf.iloc [ i , geom_col ]
            lines [ : ] = [ ]
            aa = 0
            j += 1
        if i == len ( gdf ) - 1 :
            print ( 'fin' )
            lines.append ( gdf.iloc [ i , geom_col ] )
            aa = aa + gdf.iloc [ i , dist_col ]
            paths.loc [ j , 'MMSI' ] = gdf.iloc [ i , mmsi_col ]
            paths.loc [ j , 'Ending Date' ] = gdf.iloc [ i , date_col ]
            if len ( lines ) > 1 :
                a = LineString ( list ( lines ) )
                paths.loc [ j , 'Trip_Length' ] = aa
                paths.loc [ j , 'geometry' ] = a
            else :
                paths.loc [ j , 'Begining Date' ] = gdf.iloc [ i , date_col ]
                paths.loc [ j , 'Trip_Length' ] = aa
                paths.loc [ j , 'geometry' ] = gdf.iloc [ i , geom_col ]
    ### In the paths there are isolated points and Linestrings with the same points twice -> Erase
    paths = paths [ (paths [ 'geometry' ].geom_type == 'LineString') ]
    linepoint = [ ]
    for i in range ( 0 , len ( paths ) ) :
        if paths.geometry.iloc [ i ].coords [ 0 ] == paths.geometry.iloc [ i ].coords [
            -1 ] :  # if the length of the LineString is only 2 points
            linepoint.append ( False )
        else :
            linepoint.append ( True )
    paths = paths [ (linepoint) ]
    lengthlinestring = [ ]
    geom_col_paths = paths.columns.get_loc ( "geometry" )
    for i in range ( 0 , len ( paths ) ) :
        a = len ( paths.iloc [ i , geom_col_paths ].coords )
        lengthlinestring.append ( a )
    paths [ 'lengthlinestring' ] = lengthlinestring
    return paths


# paths_ports for statistics

### From the paths dataframe we want to select the ones related to ports of Azores
# def paths_ports ( gdf, paths ) :
#     geom_col = paths.columns.get_loc ( "geometry" )
#     anyport = [ ]
#     for i in range ( 0 , len ( paths ) ) :
#         if paths.iloc [ i , geom_col ].boundary [ 0 ].within ( p ) or paths.iloc [ i , geom_col ].boundary [
#             -1 ].within ( p ) :
#             anyport.append ( True )
#         else :
#             anyport.append ( False )
#     noports = [ not c for c in anyport ]
#     paths_ports = paths [ (anyport) ]
#     paths_noports = paths [ (noports) ]
#     lengthlinestring = [ ]
#     for i in range ( 0 , len ( paths ) ) :
#         a = len ( paths.iloc [ i , geom_col ].coords )
#         lengthlinestring.append ( a )
#     paths [ 'lengthlinestring' ] = lengthlinestring
#     mean_paths = [ ]  # Empty list to store a mean
#     median_paths = [ ]  # Empty list to store a median
#     av_relvel = [ ]
#     av_speed = [ ]
#     av_length = [ ]
#     av_time = [ ]
#     for i in range ( 0 , len ( paths ) ) :
#         time_col = gdf.columns.get_loc ( "TimeLaps" )
#         dist_col = gdf.columns.get_loc ( "Dist" )
#         # if the date is between the starting and ending date of the path and has the same MMSI
#         test = gdf [ (gdf [ 'DATE' ] <= paths.iloc [ i , 2 ]) & (gdf [ 'DATE' ] >= paths.iloc [ i , 1 ]) & (
#                 gdf [ 'MMSI' ] == paths.iloc [ i , 0 ]) ]
#         test.iloc [ 0 , time_col ] = 0
#         test.iloc [ 0 , dist_col ] = 0
#         mean_t = np.mean ( abs ( test [ 'TimeLaps' ] ) )
#         mean_l = np.mean ( abs ( test [ 'Dist' ] ) )
#         mean_relvel = np.mean ( abs ( test [ 'RelVel' ] ) )
#         mean_sp = np.mean ( abs ( test [ 'SPEED' ] ) )
#         av_relvel.append ( mean_relvel )
#         av_speed.append ( mean_sp )
#         av_length.append ( mean_l )
#         av_time.append ( mean_t )
#
#     port_departure = [ ]
#     port_arrival = [ ]
#     geom_col = paths.columns.get_loc ( "geometry" )
#     for i in range ( 0 , len ( paths ) ) :
#         if paths.iloc [ i , geom_col ].boundary [ 0 ].within ( p ) :
#             for j in range ( 0 , len ( ports_azores ) ) :
#                 if paths.iloc [ i , geom_col ].boundary [ 0 ].within ( ports_azores.iloc [ j , -1 ] ) :
#                     port_departure.append ( ports_azores [ 'ports' ] [ j ] )
#         else :
#             port_departure.append ( np.nan )
#     for i in range ( 0 , len ( paths ) ) :
#         if paths.iloc [ i , geom_col ].boundary [ 1 ].within ( p ) :
#             for j in range ( 0 , len ( ports_azores ) ) :
#                 if paths.iloc [ i , geom_col ].boundary [ 1 ].within ( ports_azores.iloc [ j , -1 ] ) :
#                     port_arrival.append ( ports_azores [ 'ports' ] [ j ] )
#         else :
#             port_arrival.append ( np.nan )
#
#     paths [ 'port_arrival' ] = port_arrival
#     paths [ 'port_departure' ] = port_departure
#     paths [ 'average_speed' ] = av_speed
#     paths [ 'average_relvel' ] = av_relvel
#     paths [ 'average_length' ] = av_length
#     paths [ 'average_time' ] = av_time
#
#     portsbothends = paths [ (paths [ 'port_arrival' ] != 'nan') & (paths [ 'port_departure' ] != 'nan') ]
#     sameport = portsbothends [ (portsbothends [ 'port_arrival' ] == portsbothends [ 'port_departure' ]) ]
#     diffport = portsbothends [ (portsbothends [ 'port_departure' ] != portsbothends [ 'port_arrival' ]) ]
#     leavingport = paths [ (paths [ 'port_departure' ] != 'nan') & (paths [ 'port_arrival' ] == 'nan') ]
#     enteringports = paths [ (paths [ 'port_arrival' ] != 'nan') & (paths [ 'port_departure' ] == 'nan') ]
#     noport = paths [ (paths [ 'port_arrival' ] == 'nan') & (paths [ 'port_departure' ] == 'nan') ]
#
#     statsameport = (len ( sameport ) / len ( paths )) * 100
#     statleavingport = (len ( leavingport ) / len ( paths )) * 100
#     statenteringport = (len ( enteringports ) / len ( paths )) * 100
#     statdiffports = (len ( diffport ) / len ( paths )) * 100
#     statnoport = (len ( noport ) / len ( paths )) * 100
#
#     return statsameport , statleavingport , statenteringport , statdiffports , portsbothends , diffport , sameport , leavingport , enteringports , paths , statnoport , noport


def paths_ports ( paths ) :
# This function adds ports of departure/arrival in the paths dataframe as extra columns
# ports come from a dataframe in the drive folder
    port_departure = [ ]
    port_arrival = [ ]
    geom_col = paths.columns.get_loc ( "geometry" )
    for i in range ( 0 , len ( paths ) ) :
        if paths.iloc [ i , geom_col ].boundary [ 0 ].within ( p ) :
            for j in range ( 0 , len ( ports ) ) :
                if paths.iloc [ i , geom_col ].boundary [ 0 ].within ( ports.iloc [ j , -1 ] ) :
                    port_departure.append ( ports [ 'PORT' ] [ j ] )
        else :
            port_departure.append ( np.nan )
    for i in range ( 0 , len ( paths ) ) :
        if paths.iloc [ i , geom_col ].boundary [ 1 ].within ( p ) :
            for j in range ( 0 , len ( ports ) ) :
                if paths.iloc [ i , geom_col ].boundary [ 1 ].within ( ports.iloc [ j , -1 ] ) :
                    port_arrival.append ( ports [ 'PORT' ] [ j ] )
        else :
            port_arrival.append ( np.nan )
    paths [ 'port_arrival' ] = port_arrival
    paths [ 'port_departure' ] = port_departure
    return paths


def fisheries ( dfsame , dataframe ) :
# This function is looking into trips and when there are multiple emissions below 2.5 knots it creates a ramp
#if the ramp is long enough, it means that there is probably a fishing event
    global dffishing , dfnotfishing
    begdate_col = dfsame.columns.get_loc ( "Begining Date" )
    enddate_col = dfsame.columns.get_loc ( "Ending Date" )
    mmsi_col = dfsame.columns.get_loc('MMSI')
    fishing = [ ]
    notfish = [ ]
    len_ind = [ ]
    len_test = [ ]
    len_time = [ ]
    for i in range ( len ( dfsame ) ) :
        test = dataframe [
            (dataframe [ 'DATE' ] <= dfsame.iloc [ i , enddate_col ]) & (dataframe [ 'DATE' ] >= dfsame.iloc [ i , begdate_col ]) & (
                    dataframe [ 'MMSI' ] == dfsame.iloc [ i , mmsi_col ]) ]
        ind = [ ]
        RelVel_col = test.columns.get_loc ( "RelVel" )

        ind.append ( 0 )
        for j in range ( 2 , len ( test ) - 3 ) :
            if (test.iloc [ j , RelVel_col ] > 3.01 and test.iloc [ j + 1 , RelVel_col ] < 3.01) :
                ind.append ( j )
            if (test.iloc [ j , RelVel_col ] < 3.01 and test.iloc [ j + 1 , RelVel_col ] > 3.01) :
                ind.append ( j )
        ind.append ( len ( test ) )
        if len ( ind ) >= 5 :
            len_ind.append ( len ( ind ) )
            len_test.append ( len ( test ) )
            len_time.append ( test [ 'TimeLaps' ].sum ( ) )
            test = test.reset_index ( )
            val_ramp = [ ]
            len_ramp = [ ]
            pos_ramp = [ ]
            d = 0
            c = 0
            # plt.figure ( )
            # plt.plot ( test.loc [ : , 'RelVel' ] )
            # plt.plot ( test.loc [ : , 'SPEED' ] )
            for k in range ( 0 , len ( ind ) - 1 ) :
                parts = test.iloc [ ind [ k ] + 1 :ind [ k + 1 ] , : ]
                ramp = parts [ 'RelVel' ].mean ( )
                val_ramp.append ( ramp )
                c = c + parts [ 'TimeLaps' ].sum ( )
                d = len ( parts )
                len_ramp.append ( c )
                pos_ramp.append ( d )
                b = np.linspace ( ramp , ramp , ind [ k + 1 ] - ind [ k ] )
                a = np.linspace ( ind [ k ] + 1 , ind [ k + 1 ] , ind [ k + 1 ] - ind [ k ] )
                # plt.plot ( a , b , '-r' )
                moomin = pd.DataFrame ( { 'positions_ramp' : pos_ramp , 'val_ramp' : val_ramp } )
                diff = [ ]
                diff = [ len_ramp [ 0 ] / 60 ]
                for l in range ( 0 , len ( len_ramp ) - 1 ) :
                    diff.append ( (len_ramp [ l + 1 ] - len_ramp [ l ]) / 60 )
                moomin [ 'duration_ramp (min)' ] = diff
                low = [ ]
                long = [ ]
            if len ( moomin.loc [ moomin [ 'val_ramp' ] < 2.5 ] ) >= 1 :
                fishing.append ( i )
                # print ( i )
                # print ( moomin)
            else :
                notfish.append ( i )
        else :
            notfish.append ( i )
    dffishing = dfsame.loc [ fishing ]
    dffishing = pd.DataFrame ( dffishing )
    dfnotfishing = dfsame.loc [ notfish ]
    dfnotfishing = pd.DataFrame ( dfnotfishing )


def subsets ( list_of_MMSI ) :
# This function takes a list of MMSI (given by IPMA) and only select those in the dataframes
    global subset
    main_data_folder = r'G:\Mon Drive\E-SHAPE_to_share\0_Data'
    # vessel_database = main_data_folder + r'\11_vessel_data_match\all_mmsi_and_parameters_0-2684_gear_ALL_EU.csv'
    vessel_database = main_data_folder + r'\11_vessel_data_match\combineOnlyPortugal_all_mmsi_and_parameters_DGRM_vesselCode.csv'
    vessel_database = pd.read_csv ( vessel_database , sep = ';' )
    vessel_database = vessel_database.set_index ( 'MMSI' )
    subset = pd.DataFrame (
        columns = [ 'IMO' , 'Name' , 'Call_Sign' , 'Flag' , 'main_gear' , 'sec_gear' , 'loa' , 'ton_gt' ,
                    'year_of_const' , 'COD_EMBA' ] )
    for mmsi in list_of_MMSI :
        try :
            subset = subset.append (
                vessel_database.loc [ mmsi , [ 'IMO' , 'Name' , 'Call_Sign' , 'Flag' , 'main_gear' , 'sec_gear' ,
                                               'loa' , 'ton_gt' , 'year_of_const' , 'COD_EMBA' ] ] )
        except KeyError :
            # print('%d is not a Portuguese Vessel vessel' % mmsi)
            # print('')
            database = (
                r'G:\Mon Drive\E-SHAPE_to_share\0_Data\11_vessel_data_match\all_mmsi_and_parameters_0-2684_gear_ALL_EU.csv'
            )
            database = pd.read_csv ( database )
            database = database.set_index ( 'MMSI' )
            subset = subset.append (
                database.loc [ mmsi , [ 'IMO' , 'Name' , 'Call_Sign' , 'Flag' , 'main_gear' , 'sec_gear' ,
                                        'loa' , 'ton_gt' , 'year_of_const' , 'COD_EMBA' ] ] )
    return subset


def res ( df ) :
# function to reset the index
    df = df.reset_index ( )
    df = df.iloc [ : , 1 : ]
    return df


def look_into ( dataframe , path , i ) :
# This function take the path with the index "i" and looks into the dataframe for every emission point into the
# test dataframe
    begdate_col = path.columns.get_loc ( "Begining Date" )
    enddate_col = path.columns.get_loc ( "Ending Date" )
    mmsi_col = path.columns.get_loc ( "MMSI" )
    test = dataframe [
        (dataframe [ 'DATE' ] <= path.iloc [ i , enddate_col ]) & (
                dataframe [ 'DATE' ] >= path.iloc [ i , begdate_col ]) & (
                dataframe [ 'MMSI' ] == path.iloc [ i , mmsi_col ]) ]
    if type ( test ) == gpd.GeoDataFrame :
        test = pd.DataFrame ( test )
    return test


def PT_vessels ( paths ) :
# Takes only the portuguese vessels in the paths dataframe
    global PT_paths , PT
    paths = res ( paths )
    PT = [ 204 , 255 , 263 ]  ### MMSI numbers of PT fleet [Continental, Madeira, Azores]
    mmsi_col = paths.columns.get_loc ( "MMSI" )
    mask = [ ]
    for i in range ( len ( paths ) ) :
        mms = int ( str ( paths.iloc [ i , mmsi_col ] ) [ 0 :3 ] )
        if mms in PT :
            mask.append ( i )

    PT_paths = paths.loc [ mask ]
    return PT_paths


def plot_relvel_compass ( test ) :
# Plots the relative velocity and the compass in a subplot
    fig = plt.figure ( )
    ax1 = fig.add_subplot ( 211 )
    ax1.plot ( test [ 'DATE' ] , test [ 'Compass' ] , 'b' , label = 'Compass' , )
    ax1.legend ( loc = 'upper left' )
    ax1.grid ( True )
    ax1.set ( title = "Heading of the vessel" , xlabel = "Date" , ylabel = "Degrees" )

    ax2 = fig.add_subplot ( 212 )
    ax2.plot ( test [ 'DATE' ] , test [ 'RelVel' ] , 'r' , label = 'Relative Velocity (Knots)' )
    ax2.legend ( loc = 'upper left' )
    ax2.grid ( True )
    ax2.set ( title = "Relative Velocity of the vessel" , xlabel = "Date" , ylabel = "Knots" )
    plt.show ( )


def find_near ( array1 , value ) :
# function that find the index of the nearest value in the array
    array1 = np.asarray ( array1 )
    idx = (np.abs ( array1 - value )).argmin ( )
    return idx


def closest_date ( array , date ) :
# Finds the closest date in the array
    data_venda_col = array.columns.get_loc ( "DATA_VENDA" )
    min_date = [ ]
    array_date = array [ 'DATA_VENDA' ]
    for d in array_date :
        min_date.append ( (d - date).total_seconds ( ) )
    if len ( list ( filter ( lambda x : (x < 0) , min_date ) ) ) != 0 :
        idx = find_near ( min_date , 0 )
        date = array.iloc [ idx , data_venda_col ]
    return date


def assoc_mmsi_cod_emba(mmsi,cod_emba ,paths,vessel_database):
    list_of_MMSI = list ( set ( paths [mmsi] ) )    # Set the list of mmsi in paths
    # this subset dataframe is a concatenation of the association MMSI/COD_EMBA
    assoc = pd.DataFrame ( )
    for mmsi_num in list_of_MMSI :
        assoc = pd.concat ( [ assoc , vessel_database [ (vessel_database [mmsi] == mmsi_num) ][ [mmsi,cod_emba]] ] )
    return assoc


def landing_association ( paths , fishin , vessel_database ) :
    # This function makes the association between the vessel database where can be found the MMSI Identification number
    # Coupled with the COD_EMBA Id number that exists in other dataframes
    # "paths" is the paths, "fishin" is the sales dataframe
    global lands , landings , distance_ports
    enddate_col = paths.columns.get_loc ( "Ending Date" )
    endport_col = paths.columns.get_loc ( "port_arrival" )
    triplength_col = paths.columns.get_loc ( "Trip_Length" )
    geom_col = paths.columns.get_loc ( "geometry" )
    paths = res ( paths )  # Reset the index of the paths dataframe
    fishin = fishin.rename ( columns = { 'NUMERO' : 'COD_EMBA' } )
    # We only take the fishing events from the MMSIs in the paths
    lands = pd.merge ( fishin , vessel_database [ [ 'MMSI' , 'COD_EMBA' ] ] )
    lands [ 'DATA_VENDA' ] = pd.to_datetime ( lands.DATA_VENDA )
    lands = res ( lands )
    landings = pd.DataFrame ( )
    for ii in range ( len ( paths ) ) :  # ii is the index of the path in the dataframe.
        mm = paths.iloc [ ii , 0 ]
        triplength = paths.iloc [ ii , triplength_col ]
        ending_date = paths.iloc [ ii , enddate_col ]
        ending_port = paths.iloc [ ii , endport_col ]
        sales_mmsi = lands [ (lands [ 'MMSI' ] == mm) ]
        sales_mmsi = res ( sales_mmsi )
        closest = closest_date ( sales_mmsi , ending_date )
        if (closest - ending_date).total_seconds ( ) / 86400 < 2 and (
                closest - ending_date).total_seconds ( ) > -27200 :
            sales_closest = sales_mmsi.loc [ sales_mmsi [ 'DATA_VENDA' ] == closest ]
            sales_closest [ 'find_path' ] = ii
            sales_closest [ 'ending_date' ] = ending_date
            sales_closest [ 'port_arrival' ] = ending_port
            sales_closest [ 'Trip_Length' ] = triplength
            landings = pd.concat ( [ landings , sales_closest ] )
    findpath_col = landings.columns.get_loc ( "find_path" )
    port_name_col = landings.columns.get_loc ( "PORTO_DESEMBARQUE" )
    distance_ports = [ ]
    for i in range ( len ( landings ) ) :
        distance_ports.append (
            distance_gps ( Point ( paths.iloc [ landings.iloc [ i , findpath_col ] , geom_col ].coords [ -1 ] ) ,
                           list ( ports [ (ports [ 'PORT' ] == landings.iloc [ i , port_name_col ]) ] [
                                      'geometry' ].centroid ) [ 0 ] ) )
    landings [ 'distance_ports' ] = distance_ports
    aa = landings.columns.to_list ( )
    aa.remove ( 'port_arrival' )
    aa.insert ( 4 , 'port_arrival' )
    aa.remove ( 'ending_date' )
    aa.insert ( 0 , 'ending_date' )
    landings = landings [ aa ]
    return landings


def select_good_paths(paths, gdf):
# This function looks into the complete paths dataframe, with the overall length of the trip, and creates 2 columns
# the first one is the pourcentage of emissions with a relative velocity under 2.5 knots
# the second one is the pourcentage of emissions with a timelaps between 2 emissions over 1 hour
# We chose to select the paths with at least 30 emissions, less than 80% of emissions with more than 1 hour
# and more than 15%  of emissions "possibly fishing"
    global paths_good
    relvel_under_2_5 = [ ]
    timelaps_over_1_hour = [ ]
    for i in range ( len ( paths ) ) :
        test = look_into ( gdf , paths , i )
        test.iloc [ 0 , -3 ] = 0
        test_under = test [ (test [ 'RelVel' ] < 2.5) ]
        test_over = test [ (test [ 'TimeLaps' ] > 3600) ]
        relvel_under_2_5.append ( 100 * (len ( test_under ) / len ( paths.iloc [ i , -1 ].coords )) )
        timelaps_over_1_hour.append ( 100 * (len ( test_over ) / len ( test )) )

    paths [ 'timelaps_over_1_hour' ] = timelaps_over_1_hour
    paths [ 'relvel_under_2_5' ] = relvel_under_2_5
    paths_good = paths [ (paths [ 'timelaps_over_1_hour' ] < 80) & (paths [ 'relvel_under_2_5' ] > 15) & (
                paths [ 'lengthlinestring' ] > 30) ]


def fish_location(gdf, landings, fish_name=0):
# This function can plot the emission points that are supposely "fishing" when the sales/trips association is made
# This association links the fishes to a trip, and here we take the low relative velocity in the trip and plot it on the map
# Then if you select
    if fish_name == 0:
 #       fish_fao = list ( set ( landings [ 'COD_FAO' ] ) )
        fish_fao = ['BET','SKJ','SWO']
    else :
        if type(fish_name) == str:
            fish_fao = [fish_name]
        else:
            fish_fao = fish_name
    colors = dict ( mcolors.BASE_COLORS , **mcolors.CSS4_COLORS )
    # Sort colors by hue, saturation, value and name.
    by_hsv = sorted ( (tuple ( mcolors.rgb_to_hsv ( mcolors.to_rgba ( color ) [ :3 ] ) ) , name)
                      for name , color in colors.items ( ) )
    sorted_names = [ name for hsv , name in by_hsv ]
    for i in range ( len ( fish_fao ) ) :
        if i == 0 :
            eez_shp.plot ( )
            if len ( fish_fao ) == 1 :
                name_fish = fish_name
        name_fish = str ( fish_fao [ i ] )
        landing_select = landings [ (landings [ 'COD_FAO' ] == name_fish) ]
        landing_select = res ( landing_select )
        b = pd.DataFrame ( )
        if len ( landing_select ) > 1 :
            color = random.choice ( sorted_names )
            for i in range ( len ( landing_select ) ) :
                test = look_into ( gdf , paths , landing_select.loc [ i , 'find_path' ] )
                test_fishing = test [ (test [ 'RelVel' ] < 2.5) ]
                b = pd.concat ( [ b , test_fishing [ 'geometry' ] ] )
            b = list ( b [ 0 ] )
            xs = [ point.x for point in b ]
            ys = [ point.y for point in b ]
            plt.scatter ( xs , ys , c = color , alpha = 0.5 , label = name_fish )
            plt.legend ( )


def geom_col_end(df):
    aa = df.columns.to_list ( )
    aa.remove ( 'geometry' )
    aa.append ( 'geometry' )
    df = df [ aa ]
    return df


def landings_geo(year):
    global fishin_pt_cont,fishin_az,fishin_mad
    dfport = pd.read_csv ( 'ports_PT_GPS.csv' )
    dfport = dfport.iloc [ : , 1 : ]
    gdfport = gpd.GeoDataFrame ( dfport , geometry = gpd.points_from_xy ( dfport.LON , dfport.LAT ) )
    del gdfport [ 'LON' ]
    del gdfport [ 'LAT' ]
    fishin = pd.read_excel (
        main_path + r'\14_IPMA_LandingDeclarations\DADOS_ENVIADOS_DEIMOS\DESEMBARQUES_%s.xlsx' % (str ( year )) )
    madeiras = Polygon ( [ (-17.83161 , 33.355298498597385) , (-15.810131068255945 , 33.355298498597385) ,
                           (-15.810131068255945 , 32.26575824509192) , (-17.831615443255945 , 32.26575824509192) ,
                           (-17.831615443255945 , 33.355298498597385) ] )
    PT_continental = Polygon ( [ (-10.5529112 , 42.491914) , (-6.4220518 , 42.491914) , (-6.4220518 , 36.5312314) ,
                                 (-10.5529112 , 36.5312314) , (-10.5529112 , 42.491914) ] )
    azores = Polygon ( [ (-29.84233865556249 , 39.95774470243776) , (-22.59136209306249 , 39.95774470243776) ,
                         (-22.59136209306249 , 36.75383990140275) , (-29.84233865556249 , 36.75383990140275) ,
                         (-29.84233865556249 , 39.95774470243776) ] )
    port_az = [ ]
    port_mad = [ ]
    port_pt_cont = [ ]
    for i in range ( len ( gdfport ) ) :
        if gdfport.iloc [ i , 1 ].within ( azores ) :
            port_az.append ( i )
        if gdfport.iloc [ i , 1 ].within ( madeiras ) :
            port_mad.append ( i )
        if gdfport.iloc [ i , 1 ].within ( PT_continental ) :
            port_pt_cont.append ( i )
    ports_az = dfport.loc [ port_az ]
    ports_mad = dfport.loc [ port_mad ]
    ports_pt_cont = dfport.loc [ port_pt_cont ]
    liste_az = [ ]
    liste_mad = [ ]
    liste_pt_cont = [ ]
    for i in range ( len ( fishin ) ) :
        if fishin.loc [ i , 'PORTO_DESEMBARQUE' ] in list(ports_az['0']) :
            liste_az.append ( i )
        if fishin.loc [ i , 'PORTO_DESEMBARQUE' ] in list(ports_mad['0']) :
            liste_mad.append ( i )
        if fishin.loc [ i , 'PORTO_DESEMBARQUE' ] in list(ports_pt_cont['0']) :
            liste_pt_cont.append ( i )
    fishin_az = fishin.loc [ liste_az ]
    fishin_mad = fishin.loc [ liste_mad ]
    fishin_pt_cont = fishin.loc [ liste_pt_cont ]
    return fishin_pt_cont , fishin_az , fishin_mad


def gear_paths(dfpaths):
# associates the mmsi with their gears
    main_gear = [ ]
    sec_gear = [ ]
    main_gear_col = df_allmmsi.columns.get_loc ( "main_gear" )
    sec_gear_col = df_allmmsi.columns.get_loc ( "sec_gear" )
    mmsi_col = dfpaths.columns.get_loc ( "MMSI" )
    mmsi_col1 = df_allmmsi.columns.get_loc ( "MMSI" )
    for i in range ( len ( dfpaths ) ) :
        a = str ( df_allmmsi [ df_allmmsi [ 'MMSI' ] == dfpaths.iloc [ i , mmsi_col ] ].iloc [ mmsi_col1 , main_gear_col ] )
        b = str ( df_allmmsi [ df_allmmsi [ 'MMSI' ] == dfpaths.iloc [ i , mmsi_col ] ].iloc [ mmsi_col1 , sec_gear_col ] )
        if 'since' in a :
            a = a.replace ( ';' , ' ' )
            aa = a.split ( ' ' )
            diffdate = [ ]
            gear = [ ]
            for indice , word in enumerate ( aa ) :
                if word == 'since' :
                    diffdate.append (
                        ((pd.to_datetime ( aa [ indice + 1 ] ) - dfpaths.iloc [ i , 1 ])).total_seconds ( ) )
                    gear.append ( aa [ indice - 1 ] )
            list_of_tuples = list ( zip ( gear , diffdate ) )
            c = pd.DataFrame ( list_of_tuples , columns = [ 'gear' , 'ind' ] )
            c1 = c [ (c [ 'ind' ] < 0) ]
            main_gear.append ( c1.loc [ c1 [ 'ind' ] == c1 [ 'ind' ].max ( ) ] [ 'gear' ].iloc [ 0 ] )
        else :
            main_gear.append ( a )
        if 'since' in b :
            b = b.replace ( ';' , ' ' )
            bb = b.split ( ' ' )
            diffdate = [ ]
            gear = [ ]
            for indice , word in enumerate ( bb ) :
                if word == 'since' :
                    diffdate.append (
                        ((pd.to_datetime ( bb [ indice + 1 ] ) - dfpaths.iloc [ i , 1 ])).total_seconds ( ) )
                    gear.append ( bb [ indice - 1 ] )
            list_of_tuples = list ( zip ( gear , diffdate ) )
            d = pd.DataFrame ( list_of_tuples , columns = [ 'gear' , 'ind' ] )
            d1 = d [ (d [ 'ind' ] < 0) ]
            sec_gear.append ( d1.loc [ d1 [ 'ind' ] == d1 [ 'ind' ].max ( ) ] [ 'gear' ].iloc [ 0 ] )
        else :
            sec_gear.append ( b )
    dfpaths [ 'main_gear' ] = main_gear
    dfpaths [ 'sec_gear' ] = sec_gear

    #dfpaths = geom_col_end(dfpaths)
    return dfpaths


def subset_AIS_list_IPMA(dataframe):
    LHP_IPMA = pd.read_csv (
        r'G:\Mon Drive\E-SHAPE_to_share\0_Data\20_fishing_fleet_IPMA\LHP_28August2020_green_blue.csv' , sep = ';' )
    LLD_IPMA = pd.read_csv ( r'G:\Mon Drive\E-SHAPE_to_share\0_Data\20_fishing_fleet_IPMA\LLD_28August2020_green.csv' ,
                             sep = ';' )
    a = list ( set ( dataframe [ 'MMSI' ] ) )
    b = [element for element in a if element in  list ( LHP_IPMA.iloc [ : , -1 ] ) ]
    c = [ element for element in a if element in list ( LLD_IPMA.iloc [ : , -1 ] )  ]
    d = b + c
    dataframe = dataframe.set_index('MMSI')
    dataframe_subset = dataframe.loc[d]
    dataframe_subset = dataframe_subset.reset_index()
    return dataframe_subset

############ Command to run ##################
# from e_shape_functions import *
#
# t1 = time.time ( )
# year = 2017
# fishin = pd.read_excel (
#     main_path + r'\14_IPMA_LandingDeclarations\DADOS_ENVIADOS_DEIMOS\DESEMBARQUES_%s.xlsx' % (str ( year )) )
# path_path = os.path.join ( main_path , '21_Vessel_paths' , 'AREA_1_2', 'paths%s.csv' % (str ( year )) )
# gdf_path = os.path.join ( main_path , '21_Vessel_paths', 'AREA_1_2' , 'gdf%s.csv' % (str ( year )) )
# landings_geo ( year )
# gdf = pd.read_csv ( gdf_path )
# paths = pd.read_csv ( path_path )
# lines = [ shapely.wkt.loads ( line ) for line in paths [ 'geometry' ] ]
# paths = gpd.GeoDataFrame ( paths.drop ( columns = [ 'geometry' ] ) , geometry = lines )
# points = [ shapely.wkt.loads ( point ) for point in gdf [ 'geometry' ] ]
# gdf = gpd.GeoDataFrame ( gdf.drop ( columns = [ 'geometry' ] ) , geometry = points )
# paths = paths.iloc [ : , 1 : ]
# gdf = gdf.iloc [ : , 1 : ]
# gdf [ 'DATE' ] = pd.to_datetime ( gdf.DATE )  ## Translate every date into a datetime type
# paths [ 'Begining Date' ] = pd.to_datetime ( paths.iloc [ : , 1 ] )  ## Translate every date into a datetime type
# paths [ 'Ending Date' ] = pd.to_datetime ( paths.iloc [ : , 2 ] )  ## Translate every date into a datetime type
# dfpaths = pd.DataFrame ( paths )
# dataframe = pd.DataFrame ( gdf )
# landings = landing_association ( paths , fishin , vessel_database )
# relvel_under_2_5 = [ ]
# timelaps_over_1_hour = [ ]
# for i in range ( len ( paths ) ) :
#     test = look_into ( gdf , paths , i )
#     test.iloc [ 0 , -3 ] = 0
#     test_under = test [ (test [ 'RelVel' ] < 2.5) ]
#     test_over = test [ (test [ 'TimeLaps' ] > 3600) ]
#     relvel_under_2_5.append ( 100 * (len ( test_under ) / len ( paths.iloc [ i , -1 ].coords )) )
#     timelaps_over_1_hour.append ( 100 * (len ( test_over ) / len ( test )) )
#
# paths [ 'timelaps_over_1_hour' ] = timelaps_over_1_hour
# paths [ 'relvel_under_2_5' ] = relvel_under_2_5
# paths = geom_col_end ( paths )
# dfpaths = pd.DataFrame ( paths )
# dataframe = pd.DataFrame ( gdf )
# dfpaths_good = dfpaths [ (dfpaths [ 'timelaps_over_1_hour' ] < 80) & (dfpaths [ 'relvel_under_2_5' ] > 15) & (
#             dfpaths [ 'lengthlinestring' ] > 30) & (dfpaths [ 'Trip_Length' ] > 5) ]
# paths_good = paths [ (paths [ 'timelaps_over_1_hour' ] < 80) & (paths [ 'relvel_under_2_5' ] > 15) & (
#             paths [ 'lengthlinestring' ] > 30) & (dfpaths [ 'Trip_Length' ] > 5) ]
# landings_good = landing_association ( paths_good , fishin , vessel_database )
# t2 = time.time ( )
# print ( t2 - t1 )