import json
import geojson
import numpy as np
from matplotlib.patches import  Rectangle
import matplotlib.pyplot as plt
import openslide
import pandas as pd
import imageio
import os 
from tqdm.auto import trange

class Annotation():
    """Class for well organized geojson data generate from QuPath annotation. 
    """
    def __init__(self, label:geojson.feature.FeatureCollection):
        """Reorganize geojsondata

        Args:
            label (geojson.feature.FeatureCollection): Raw jeojson data
        """
        self.rowdata = label 
        self.features = label['features']

    def __len__(self):
        return len(self.features)
    
    def get_coordinates(self, index):
        """Get No.[index] coordinates of the rectangle

        Args:
            index (int): index of self.features

        Returns:
            np.Array: 4 vertexes of the rectangle
        """
        return np.array(self.features[index]['geometry']['coordinates']).squeeze()[1:]
    
    def get_tissue_type(self, index): 
        """Get No.[index] tissue type of the rectang;e

        Args:
            index (int): index of self.features

        Returns:
            str: tissure type
        """
        return self.features[index]['properties']['classification']['name']
    
    def get_position(self, index):
        """Summarize min and max x,y coordinates

        Args:
            index (int): index of self.features

        Returns:
            minx, miny, maxx, maxy: min and max x,y coordinates
        """
        coordinate = self.get_coordinates(index)
        minx, miny = coordinate[-1][0], coordinate[-1][1]
        maxx, maxy = coordinate[1][0], coordinate[1][1]
        return minx, miny, maxx, maxy
    
    def get_summarize_df(self):
        """Summarize all ROI in form of pd.Dataframe

        Returns:
            pd.Dataframe: with 5 columns 'name', 'x': minx, 'y': miny, 'w': maxx-minx, 'h': maxy-miny
        """
        x, y, w, h, name = [], [], [], [], []
        for i in range(len(self.features)):
            name.append(self.get_tissue_type(i))
            minx, miny, maxx, maxy = self.get_position(i)
            x.append(minx)
            y.append(miny)
            w.append(maxx - minx)
            h.append(maxy - miny)
        df = pd.DataFrame({
            'name': name, 
            'x': x, 
            'y': y, 
            'w': w,
            'h': h
        })
        return df 
    
    def generate_patches(self, slide_path:str, save_root:str, size:int=112, tissue_type=None):
        """Generate corresponding patches for each ROI from WSI. 

        Args:
            slide_path (str): Path to corresponding WSI file
            save_root (str): Path to save dir 
            size (int, optional): Size of each patch. Defaults to 112.
            tissue_type (list, optional): List of tissue types of interest. Defaults to None.

        Returns:
            pd.Dataframe: Summize result of patch generation. 
        """
        slide = openslide.OpenSlide(slide_path)
        slide_name = slide_path.split('/')[-1]
        if os.path.exists(os.path.join(save_root, slide_name)) is not True: 
                os.mkdir(os.path.join(save_root, slide_name))
        patch_num = []
        df = self.get_summarize_df()
        if tissue_type is not None: 
            df = df[df['name'].isin(tissue_type)]
        for i in trange(len(df)):
            name, x, y, w, h = tuple(df.iloc[i])
            save_path = os.path.join(save_root, slide_name, name)
            if os.path.exists(save_path) is not True: 
                os.mkdir(save_path)
            nx = int(np.floor(w/size))
            ny = int(np.floor(h/size))
            patch_num.append(nx * ny)
            for n in trange(x, x + nx * size, size, desc='x', leave=False): 
                for m in trange(y, y + ny * size, size, desc='y', leave=False):
                    im = np.array(slide.read_region((n, m), 0, (size, size)))
                    imageio.imwrite(os.path.join(save_path, slide_name + f'_({n},{m}).png'), im)
        df['patch_num'] = patch_num
        slide.close()
        return df