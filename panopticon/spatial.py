import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_labels_and_polys(tif_file):
    from panopticon.utilities import import_check
    exit_code = import_check(
        "tifffile", 'pip install tifffile')
    if exit_code != 0:
        return
    exit_code = import_check(
        "stardist", 'pip install stardist')
    if exit_code != 0:
        return
    from tifffile import imread, imwrite, imshow
    from stardist.models import StarDist2D
    img = imread(tif_file)
    img = img.swapaxes(0,2) #Assumes to be CXY, we switch to YXC
    model = StarDist2D.from_pretrained('2D_versatile_he')
    
    img = img/255
    labels, polys = model.predict_instances_big(img, axes='YXC', block_size=4096, prob_thresh=0.01,nms_thresh=0.001, min_overlap=128, context=128, normalizer=None, n_tiles=(4,4,1))
    return labels, polys

def get_gdf(labels, polys):
    from panopticon.utilities import import_check
    exit_code = import_check(
        "geopandas", 'pip install geopandas')
    if exit_code != 0:
        return
    import geopandas as gpd

    # Creating a list to store Polygon geometries
    geometries = []
    
    # Iterating through each nuclei in the 'polys' DataFrame
    for nuclei in range(len(polys['coord'])):
    
        # Extracting coordinates for the current nuclei and converting them to (y, x) format
        coords = [(y, x) for x, y in zip(polys['coord'][nuclei][0], polys['coord'][nuclei][1])]
    
        # Creating a Polygon geometry from the coordinates
        geometries.append(Polygon(coords))
    
    # Creating a GeoDataFrame using the Polygon geometries
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf['id'] = [f"ID_{i+1}" for i, _ in enumerate(gdf.index)]
    return gdf
    
