import folium
from folium.plugins import HeatMap 
import pandas as pd
import googlemaps 
#import continents
import os
import xyzservices.providers as xyz
import pickle 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch
from matplotlib.gridspec import GridSpec

api_key = 'AIzaSyA39Qg4fyA86_bg-tSF-Rr7nHsJfb5VTOc'
stadia_key = '138f93a3-0b14-4cdc-b13e-832b7c5eb7a1'
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
    db = db.groupby('Institution').size().reset_index(name='Count')
    location = pickle.load(open('src/data/institution_locations.pkl', 'rb'))
    db['Lat'] = db['Institution'].map(lambda x: location.get(x, {'lat':0,'lng':0})['lat'])
    db['Long'] = db['Institution'].map(lambda x: location.get(x, {'lat':0,'lng':0})['lng'])


    located =  db.drop(db[(db['Lat'] == 0) & (db['Long'] == 0)].index)
    return located

def get_geocode(institution):
    geocode_result = gmaps.geocode(institution)
    return geocode_result



def load_data():
    if not os.path.exists('src/data/BindingDB_All_202411.pkl') and os.path.exists('src/data/BindingDB_All_202411_tsv.zip'):
        df = pd.read_csv('src/data/BindingDB_All_202411_tsv.zip',sep='\t',compression='zip',on_bad_lines='skip',low_memory=False)
        pickle.dump(df,open('src/data/BindingDB_All_202411.pkl','wb'))
    # if 
    else:
        df = pd.read_pickle('src/data/BindingDB_All_202411.pkl')

    return df


def continents_pie_chart(dy) :
    continent_institutions = dy.groupby('Continent')['Institution'].count()
    total_institutions = continent_institutions.sum()
    continent_institutions = continent_institutions[continent_institutions >= 0.05 * total_institutions]
    continent_institutions['Others'] = total_institutions - continent_institutions.sum()

    plt.figure(figsize=(8, 8))
    plt.pie(continent_institutions, labels=continent_institutions.index, autopct='%1.1f%%', startangle=140)
    plt.title('Institutions by Continent')
    plt.show()

def top_countries_bar_chart(dy) :
    top_10_countries = dy.groupby('Country')['Count'].sum().nlargest(10)
    plt.figure(figsize=(12, 6))
    plt.bar(top_10_countries.index, top_10_countries.values)
    plt.title('Top 10 Countries by Number of Contributions')
    plt.xlabel('Country')
    plt.ylabel('Number of Contributions')
    plt.xticks(rotation=45)
    plt.show()

def contribs_by_continents_pie(dy) :
    # Pie chart for number of contributions by continent, grouping those with < 5% as "Others"
    continent_contribs = dy.groupby('Continent')['Count'].sum()
    total_contribs = continent_contribs.sum()
    continent_contribs = continent_contribs[continent_contribs >= 0.05 * total_contribs]
    continent_contribs['Others'] = total_contribs - continent_contribs.sum()

    plt.figure(figsize=(8, 8))
    plt.pie(continent_contribs, labels=continent_contribs.index, autopct='%1.1f%%', startangle=140)
    plt.title('Contributions by Continent')
    plt.show()
    
def countries_with_most_institutions(dy) :
    # Displaying a bar chart for number of institutions by country
    top_10_institutions = dy.groupby('Country')['Institution'].count().nlargest(10)
    plt.figure(figsize=(12, 6))
    plt.bar(top_10_institutions.index, top_10_institutions.values)
    plt.title('Top 10 Countries by Number of Institutions')
    plt.xlabel('Country')
    plt.ylabel('Number of Institutions')
    plt.xticks(rotation=45)
    plt.show()

def institutions_per_continent(dy) :
    continent_institutions_full = dy.groupby('Continent')['Institution'].count()

    plt.figure(figsize=(12, 6))
    plt.bar(continent_institutions_full.index, continent_institutions_full.values)
    plt.title('Number of Institutions per Continent')
    plt.xlabel('Continent')
    plt.ylabel('Number of Contributing Institutions')
    plt.xticks(rotation=45)
    plt.show()

def get_institution_locations(df):
    full_locations = pickle.load(open('src/data/full_institution_locations.pkl', 'rb'))

    df['Institution'].value_counts()
    db = df['Institution'].value_counts()
    db = db.reset_index()
    db.columns = ['Institution', 'Count']
    countries ={}
    import src.utils.continents as continents
    for k,v in full_locations.items():
            components = v[0]['address_components']
            country = next((component['long_name'] for component in components if 'country' in component['types']), 'Unknown')
            countries[k] = country
            
    db['Country'] = db['Institution'].map(countries)
    db['Continent'] = db['Country'].map(continents.get_continent)
    db = db[db['Country'] != 'Unknown']
    db = db[db['Continent'] != 'Unknown']
    return db