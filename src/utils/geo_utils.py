import folium
from folium.plugins import HeatMap 
import pandas as pd
import googlemaps 
import continents
import os
import xyzservices.providers as xyz
import pickle 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch
from matplotlib.gridspec import GridSpec

api_key = os.getenv('GMAPS_API_KEY')
stadia_key = os.getenv('STADIA_KEY')
gmaps = googlemaps.Client(key=api_key)

def get_map(located,tiles="Stadia.AlidadeSmooth"):
    
    stadia = xyz.Stadia.AlidadeSmooth(api_key=stadia_key)  

    stadia["url"] = stadia["url"] + "?api_key={api_key}"

    m= folium.Map(location=[located['Lat'].mean(), located['Long'].mean()], zoom_start=2,max_bounds=True, min_zoom=2,tiles=stadia.build_url(),
    attr=stadia["attribution"])
    m.fit_bounds([[located['Lat'].min(), located['Long'].min()], [located['Lat'].max(), located['Long'].max()]])


    heat_data = [[row['Lat'], row['Long'], row['Count']] for index, row in located.iterrows()]
    HeatMap(heat_data,radius=8,blur=5,min_opacity=0.35).add_to(m)

    return m

def get_cached_locations(db):
    location = pickle.load(open('../data/institution_locations.pkl', 'rb'))
    db['Lat'] = db['Institution'].map(lambda x: location.get(x, {'lat':0,'lng':0})['lat'])
    db['Long'] = db['Institution'].map(lambda x: location.get(x, {'lat':0,'lng':0})['lng'])


    located =  db.drop(db[(db['Lat'] == 0) & (db['Long'] == 0)].index)
    return located

def get_geocode(institution):
    geocode_result = gmaps.geocode(institution)
    return geocode_result



def load_data():
    if not os.path.exists('../../data/BindingDB_All_202411.pkl') and os.path.exists('../../data/BindingDB_All_202411_tsv.zip'):
        df = pd.read_csv('../../data/BindingDB_All_202411_tsv.zip',sep='\t',compression='zip',on_bad_lines='skip',low_memory=False)
        pickle.dump(df,open('../../data/BindingDB_All_202411.pkl','wb'))
    # if 
    else:
        df = pd.read_pickle('../../data/BindingDB_All_202411.pkl')
    return df