import pandas as pd
import pickle
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
import plotly.express as px
import plotly.io as pio
import plotly.graph_objs as go
import statsmodels.formula.api as smf

code_to_region_WHO ={
 'ABW': 'AMR',
 'AFG': 'EMR',
 'AGO': 'AFR',
 'AIA': 'AMR',
 'ALB': 'EUR',
 'AND': 'EUR',
 'ARE': 'EMR',
 'ARG': 'AMR',
 'ARM': 'EUR',
 'ASM': 'WPR',
 'ATG': 'AMR',
 'AUS': 'WPR',
 'AUT': 'EUR',
 'AZE': 'EUR',
 'BDI': 'AFR',
 'BEL': 'EUR',
 'BEN': 'AFR',
 'BES': 'AMR',
 'BFA': 'AFR',
 'BGD': 'SEAR',
 'BGR': 'EUR',
 'BHR': 'EMR',
 'BHS': 'AMR',
 'BIH': 'EUR',
 'BLR': 'EUR',
 'BLZ': 'AMR',
 'BMU': 'AMR',
 'BOL': 'AMR',
 'BRA': 'AMR',
 'BRB': 'AMR',
 'BRN': 'WPR',
 'BTN': 'SEAR',
 'BWA': 'AFR',
 'CAF': 'AFR',
 'CAN': 'AMR',
 'CHE': 'EUR',
 'CHL': 'AMR',
 'CHN': 'WPR',
 'CIV': 'AFR',
 'CMR': 'AFR',
 'COD': 'AFR',
 'COG': 'AFR',
 'COK': 'WPR',
 'COL': 'AMR',
 'COM': 'AFR',
 'CPV': 'AFR',
 'CRI': 'AMR',
 'CUB': 'AMR',
 'CUW': 'AMR',
 'CYM': 'AMR',
 'CYP': 'EUR',
 'CZE': 'EUR',
 'DEU': 'EUR',
 'DJI': 'EMR',
 'DMA': 'AMR',
 'DNK': 'EUR',
 'DOM': 'AMR',
 'DZA': 'AFR',
 'ECU': 'AMR',
 'EGY': 'EMR',
 'ERI': 'AFR',
 'ESP': 'EUR',
 'EST': 'EUR',
 'ETH': 'AFR',
 'FIN': 'EUR',
 'FJI': 'WPR',
 'FLK': 'AMR',
 'FRA': 'EUR',
 'FRO': 'EUR',
 'FSM': 'WPR',
 'GAB': 'AFR',
 'GBR': 'EUR',
 'GEO': 'EUR',
 'GHA': 'AFR',
 'GIB': 'EUR',
 'GIN': 'AFR',
 'GLP': 'AMR',
 'GMB': 'AFR',
 'GNB': 'AFR',
 'GNQ': 'AFR',
 'GRC': 'EUR',
 'GRD': 'AMR',
 'GRL': 'EUR',
 'GTM': 'AMR',
 'GUF': 'AMR',
 'GUM': 'WPR',
 'GUY': 'AMR',
 'HND': 'AMR',
 'HRV': 'EUR',
 'HTI': 'AMR',
 'HUN': 'EUR',
 'IDN': 'SEAR',
 'IMN': 'EUR',
 'IND': 'SEAR',
 'IRL': 'EUR',
 'IRN': 'EMR',
 'IRQ': 'EMR',
 'ISL': 'EUR',
 'ISR': 'EUR',
 'ITA': 'EUR',
 'JAM': 'AMR',
 'JOR': 'EMR',
 'JPN': 'WPR',
 'KAZ': 'EUR',
 'KEN': 'AFR',
 'KGZ': 'EUR',
 'KHM': 'WPR',
 'KIR': 'WPR',
 'KNA': 'AMR',
 'KOR': 'WPR',
 'KWT': 'EMR',
 'LAO': 'WPR',
 'LBN': 'EMR',
 'LBR': 'AFR',
 'LBY': 'EMR',
 'LCA': 'AMR',
 'LIE': 'EUR',
 'LKA': 'SEAR',
 'LSO': 'AFR',
 'LTU': 'EUR',
 'LUX': 'EUR',
 'LVA': 'EUR',
 'MAR': 'EMR',
 'MCO': 'EUR',
 'MDA': 'EUR',
 'MDG': 'AFR',
 'MDV': 'SEAR',
 'MEX': 'AMR',
 'MHL': 'WPR',
 'MKD': 'EUR',
 'MLI': 'AFR',
 'MLT': 'EUR',
 'MMR': 'SEAR',
 'MNE': 'EUR',
 'MNG': 'WPR',
 'MNP': 'WPR',
 'MOZ': 'AFR',
 'MRT': 'AFR',
 'MSR': 'AMR',
 'MTQ': 'AMR',
 'MUS': 'AFR',
 'MWI': 'AFR',
 'MYS': 'WPR',
 'MYT': 'AFR',
 'NAM': 'AFR',
 'NCL': 'WPR',
 'NER': 'AFR',
 'NGA': 'AFR',
 'NIC': 'AMR',
 'NIU': 'WPR',
 'NLD': 'EUR',
 'NOR': 'EUR',
 'NPL': 'SEAR',
 'NRU': 'WPR',
 'NZL': 'WPR',
 'OMN': 'EMR',
 'PAK': 'EMR',
 'PAN': 'AMR',
 'PER': 'AMR',
 'PHL': 'WPR',
 'PLW': 'WPR',
 'PNG': 'WPR',
 'POL': 'EUR',
 'PRI': 'AMR',
 'PRK': 'SEAR',
 'PRT': 'EUR',
 'PRY': 'AMR',
 'PSE': 'EMR',
 'PYF': 'WPR',
 'QAT': 'EMR',
 'REU': 'AFR',
 'ROU': 'EUR',
 'RUS': 'EUR',
 'RWA': 'AFR',
 'SAU': 'EMR',
 'SDN': 'EMR',
 'SDN736': 'EMR',
 'SEN': 'AFR',
 'SGP': 'WPR',
 'SHN': 'AFR',
 'SLB': 'WPR',
 'SLE': 'AFR',
 'SLV': 'AMR',
 'SMR': 'EUR',
 'SOM': 'EMR',
 'SPM': 'AMR',
 'SRB': 'EUR',
 'SSD': 'AFR',
 'STP': 'AFR',
 'SUR': 'AMR',
 'SVK': 'EUR',
 'SVN': 'EUR',
 'SWE': 'EUR',
 'SWZ': 'AFR',
 'SXM': 'AMR',
 'SYC': 'AFR',
 'SYR': 'EMR',
 'TCA': 'AMR',
 'TCD': 'AFR',
 'TGO': 'AFR',
 'THA': 'SEAR',
 'TJK': 'EUR',
 'TKL': 'WPR',
 'TKM': 'EUR',
 'TLS': 'SEAR',
 'TON': 'WPR',
 'TTO': 'AMR',
 'TUN': 'EMR',
 'TUR': 'EUR',
 'TUV': 'WPR',
 'TZA': 'AFR',
 'UGA': 'AFR',
 'UKR': 'EUR',
 'URY': 'AMR',
 'USA': 'AMR',
 'UZB': 'EUR',
 'VCT': 'AMR',
 'VEN': 'AMR',
 'VGB': 'AMR',
 'VIR': 'AMR',
 'VNM': 'WPR',
 'VUT': 'WPR',
 'WLF': 'WPR',
 'WSM': 'WPR',
 'YEM': 'EMR',
 'ZAF': 'AFR',
 'ZMB': 'AFR',
 'ZWE': 'AFR'
}

country_to_code =  {
    'United States': 'USA',
    'China': 'CHN',
    'Germany': 'DEU',
    'Italy': 'ITA',
    'United Kingdom': 'GBR',
    'France': 'FRA',
    'India': 'IND',
    'Japan': 'JPN',
    'South Korea': 'KOR',
    'Canada': 'CAN',
    'Switzerland': 'CHE',
    'Poland': 'POL',
    'Australia': 'AUS',
    'Belgium': 'BEL',
    'Netherlands': 'NLD',
    'Spain': 'ESP',
    'Denmark': 'DNK',
    'Sweden': 'SWE',
    'Taiwan': 'TWN',
    'Egypt': 'EGY',
    'New Zealand': 'NZL',
    'Türkiye': 'TUR',
    'Lithuania': 'LTU',
    'Czechia': 'CZE',
    'Israel': 'ISR',
    'Austria': 'AUT',
    'South Africa': 'ZAF',
    'Pakistan': 'PAK',
    'Singapore': 'SGP',
    'Slovenia': 'SVN',
    'Greece': 'GRC',
    'Finland': 'FIN',
    'Russia': 'RUS',
    'Brazil': 'BRA',
    'Cayman Islands': 'CYM',
    'Iran': 'IRN',
    'Saudi Arabia': 'SAU',
    'Portugal': 'PRT',
    'Malaysia': 'MYS',
    'Thailand': 'THA',
    'Norway': 'NOR',
    'Hungary': 'HUN',
    'Bangladesh': 'BGD',
    'Hong Kong': 'HKG',
    'Mexico': 'MEX',
    'Argentina': 'ARG',
    'Jordan': 'JOR',
    'United Arab Emirates': 'ARE',
    'Slovakia': 'SVK',
    'Ireland': 'IRL',
    'Bermuda': 'BMU',
    'Chile': 'CHL',
    'Romania': 'ROU',
    'Serbia': 'SRB',
    'Philippines': 'PHL',
    'Estonia': 'EST',
    'Ukraine': 'UKR',
    'Monaco': 'MCO',
    'Macao': 'MAC',
    'Colombia': 'COL',
    'Croatia': 'HRV',
    'Morocco': 'MAR',
    'Tunisia': 'TUN',
    'Indonesia': 'IDN',
    'Latvia': 'LVA',
    'Vietnam': 'VNM',
    'Bulgaria': 'BGR',
    'Ecuador': 'ECU',
    'Zambia': 'ZMB',
    'Venezuela': 'VEN',
    'Lebanon': 'LBN',
    'Rwanda': 'RWA',
    'Oman': 'OMN',
    'Luxembourg': 'LUX',
    'Cyprus': 'CYP',
    'Botswana': 'BWA',
    'Uruguay': 'URY',
    'Algeria': 'DZA',
    'Grenada': 'GRD',
    'Jamaica': 'JAM',
    'Kuwait': 'KWT',
    'Sudan': 'SDN',
    'Peru': 'PER',
    'Qatar': 'QAT',
    'Zimbabwe': 'ZWE',
    'Puerto Rico': 'PRI',
    'Bosnia and Herzegovina': 'BIH',
    'Nepal': 'NPL',
    'Sri Lanka': 'LKA',
    'Armenia': 'ARM',
    'Kenya': 'KEN',
    'Cuba': 'CUB',
    'Ethiopia': 'ETH',
    'Nigeria': 'NGA'
}

HIV_sources = ['Human immunodeficiency virus 1', 'Human immunodeficiency virus', 'Human immunodeficiency virus type 1 group M subtype B (isolate BRU/LAI)', 'Human immunodeficiency virus type 1 group M subtype B (isolate YU-2)', 'Human immunodeficiency virus type 1 group M subtype B (isolate HXB2)', 'Human immunodeficiency virus 2', 'Human immunodeficiency virus type 1 group M subtype B (isolate PCV12)', 'Human immunodeficiency virus type 1 group M subtype B (isolate MN)']

Plasmodium_falciparum_sources = ['Plasmodium falciparum', 'Plasmodium falciparum (isolate 3D7)', 'Plasmodium falciparum (isolate K1 / Thailand)', 'Plasmodium falciparum (isolate FcB1 / Columbia)']

Poliovirus_sources = ['Poliovirus type 1 (strain Mahoney)']

Plasmodium_vivax_sources = ['Plasmodium vivax']

Tuberculosis_sources = ['Mycobacterium tuberculosis', 'Mycobacterium tuberculosis H37Rv']

Hepatitis_C_sources = ['Hepatitis C virus', 'Hepatitis C virus genotype 1a (isolate H)', 'Hepatitis C virus genotype 1b (isolate BK)', 'Hepatitis C virus genotype 3a (isolate NZL1)', 'Hepatitis C virus genotype 2b (isolate HC-J8)', 'Hepatitis C virus genotype 1b (isolate Taiwan)', 'Hepatitis C virus genotype 2', 'Hepatitis C virus genotype 4a (isolate ED43)', 'Hepatitis C virus genotype 6a (isolate EUHK2)']

Escherichia_coli_sources = ['Escherichia coli', 'Escherichia coli str. K-12 substr. MG1655', 'Escherichia coli (strain UTI89 / UPEC)', 'Escherichia coli (strain K12)', 'Escherichia coli O157:H7', 'Escherichia coli O6']

Staphylococcus_aureus_sources = ['Staphylococcus aureus', 'Staphylococcus aureus (strain MRSA252)', 'Staphylococcus aureus (strain Mu50 / ATCC 700699)', 'Staphylococcus aureus (strain MW2)', 'Staphylococcus aureus (strain NCTC 8325)']

disease_names = ["Plasmodium Falciparum", "Human Immunodeficiency Virus",
                  "Poliovirus", "Plasmodium vivax", 
                  "Mycobacterium tuberculosis", "Escherichia coli", "Hepatitis C", "Staphylococcus aureus"]

disease_indicators_names = ["Number of indigenous P. falciparum malaria cases",
                            "Number of new HIV infections",
                            "Poliomyelitis - number of reported cases",
                            "Number of indigenous P. vivax malaria cases",
                            "Incidence of tuberculosis (per 100 000 population per year)",
                            "Proportion of bloodstream infection due to Escherichia coli resistant to third-generation cephalosporins (%)",
                            "Prevalence of chronic hepatitis C (HCV) ) in the general population (%)",
                            "Proportion of bloodstream infection due to methicillin-resistant Staphylococcus aureus (MRSA) (%)"
                            ]

disease_indicator_codes = ["MALARIA_PF_INDIG", "HIV_0000000026", "WHS3_49", "MALARIA_PV_INDIG",
                           "MDG_0000000020", "AMR_INFECT_ECOLI", "HEPATITIS_HCV_PREVALENCE_PER100", "AMR_INFECT_MRSA"]

country_dimensions = pd.read_pickle("src/data/who_country_region_codes.pkl")
country_populations = pd.read_pickle("src/data/who_country_populations.pkl")

country_codes_no_population_data = country_dimensions[~country_dimensions['Code'].isin(country_populations['SpatialDim'])]['Code']
region_population = country_populations.groupby('ParentLocationCode')['NumericValue'].sum()
region_codes = sorted(country_dimensions['ParentCode'].unique().tolist()[:-1])

diseases_data = pd.read_pickle("src/data/who_diseases_data.pkl")

pf_malaria_data = diseases_data['pf_malaria']
hiv_data = diseases_data['hiv']
polio_data = diseases_data['polio']
pv_malaria_data = diseases_data['pv_malaria']
tuberculosis_data = diseases_data['tuberculosis']
escherichia_coli_data = diseases_data['escherichia_coli']
hepatitis_c_data = diseases_data['hepatitis_c']
staphylococcus_aureus_data = diseases_data['staphylococcus_aureus']

def country_cases_to_cases_per_1000(data_df, country_codes_no_population_data, country_populations):
    data_df_country = data_df[data_df['SpatialDimType'] == 'COUNTRY']
    data_df_country = data_df_country.groupby('SpatialDim')['NumericValue'].mean()
    data_df_country = data_df_country.reset_index()
    data_df_country = data_df_country.rename(columns={'SpatialDim': 'Code', 'NumericValue': 'Cases'})

    data_df_country = data_df_country[~data_df_country['Code'].isin(country_codes_no_population_data)]
    data_df_country['Cases per 1000'] = data_df_country.apply(lambda row: row['Cases'] / country_populations[country_populations['SpatialDim'] == row['Code']]['NumericValue'].values[0] * 1000, axis=1)
    return data_df_country

def country_cases_to_region_cases_per_1000(data_df, region_population):
    data_df_region = data_df.groupby('ParentLocationCode')['NumericValue'].sum()
    data_df_region = data_df_region.reset_index()
    data_df_region = data_df_region.rename(columns={'ParentLocationCode': 'Code', 'NumericValue': 'Cases'})
    data_df_region['Cases per 1000'] = data_df_region.apply(lambda row: row['Cases'] / region_population[row['Code']] * 1000, axis=1)
    return data_df_region

def country_cases_per_x_to_cases_per_1000(data_df, x, country_codes_no_population_data, country_populations):
    data_df_country = data_df[data_df['SpatialDimType'] == 'COUNTRY']
    data_df_country = data_df_country.groupby('SpatialDim')['NumericValue'].mean()
    data_df_country = data_df_country.reset_index()
    data_df_country = data_df_country.rename(columns={'SpatialDim': 'Code', 'NumericValue': f'Cases per {x}'})
    data_df_country = data_df_country[~data_df_country['Code'].isin(country_codes_no_population_data)]
    data_df_country['Cases per 1000'] = data_df_country.apply(lambda row: row[f'Cases per {x}'] * (1000/x), axis=1)
    data_df_country['Cases'] = data_df_country.apply(lambda row: row['Cases per 1000'] / 1000 * country_populations[country_populations['SpatialDim'] == row['Code']]['NumericValue'].values[0], axis=1)
    data_df_country['NumericValue'] = data_df_country['Cases']
    data_df_country['ParentLocationCode'] = data_df_country['Code'].apply(country_code_to_region_code)
    return data_df_country

def country_code_to_region_code(country_code):
    return country_dimensions[country_dimensions['Code'] == country_code]['ParentCode'].values[0]

def plot_who_heatmap():
    pf_malaria_country_cases = country_cases_to_cases_per_1000(pf_malaria_data, country_codes_no_population_data, country_populations)
    pf_malaria_region_cases = country_cases_to_region_cases_per_1000(pf_malaria_data, region_population)

    hiv_country_cases = country_cases_to_cases_per_1000(hiv_data, country_codes_no_population_data, country_populations)
    hiv_region_cases = country_cases_to_region_cases_per_1000(hiv_data, region_population)

    polio_country_cases = country_cases_to_cases_per_1000(polio_data, country_codes_no_population_data, country_populations)
    polio_region_cases = country_cases_to_region_cases_per_1000(polio_data, region_population)

    pv_malaria_country_cases = country_cases_to_cases_per_1000(pv_malaria_data, country_codes_no_population_data, country_populations)
    pv_malaria_region_cases = country_cases_to_region_cases_per_1000(pv_malaria_data, region_population)

    tuberculosis_country_cases = country_cases_per_x_to_cases_per_1000(tuberculosis_data, 1000, country_codes_no_population_data, country_populations)
    tuberculosis_region_cases = country_cases_to_region_cases_per_1000(tuberculosis_data, region_population)

    escherichia_coli_country_cases = country_cases_per_x_to_cases_per_1000(escherichia_coli_data, 100, country_codes_no_population_data, country_populations)
    escherichia_coli_region_cases = country_cases_to_region_cases_per_1000(escherichia_coli_data, region_population)

    hepatitis_c_country_cases = country_cases_per_x_to_cases_per_1000(hepatitis_c_data, 100, country_codes_no_population_data, country_populations)
    hepatitis_c_region_cases = country_cases_to_region_cases_per_1000(hepatitis_c_data, region_population)

    staphylococcus_aureus_country_cases = country_cases_per_x_to_cases_per_1000(staphylococcus_aureus_data, 100, country_codes_no_population_data, country_populations)
    staphylococcus_aureus_region_cases = country_cases_to_region_cases_per_1000(staphylococcus_aureus_data, region_population)
    print('sa')
    all_region_cases = [pf_malaria_region_cases, hiv_region_cases, polio_region_cases, pv_malaria_region_cases, tuberculosis_region_cases, escherichia_coli_region_cases, hepatitis_c_region_cases, staphylococcus_aureus_region_cases]

    for i in range(len(all_region_cases)):
        all_region_cases[i] = all_region_cases[i].rename(columns={'Cases per 1000': f'{disease_names[i]}'})

    all_region_cases = pd.concat(all_region_cases, axis=1)
    all_region_cases.drop(columns=['Cases'], inplace=True)

    code_column = all_region_cases['Code']
    first_code_column = code_column.iloc[:, 0]
    all_region_cases.drop(columns=['Code'], inplace=True)

    # scaler = StandardScaler()
    # all_region_cases = pd.DataFrame(scaler.fit_transform(all_region_cases), columns=all_region_cases.columns)

    scaler = MinMaxScaler()
    all_region_cases = pd.DataFrame(scaler.fit_transform(all_region_cases), columns=all_region_cases.columns)

    all_region_cases.insert(0, 'Cases', first_code_column)
    print(all_region_cases)

    heatmap_data = all_region_cases.set_index('Cases').T
    plt.figure(figsize=(10, 7))
    sns.heatmap(heatmap_data, fmt=".2f", cmap='coolwarm')
    plt.show()

def df_institutions_target_country_region_WHO(df,locations) :
    dx = df[['Institution', 'Target Source Organism According to Curator or DataSource']].copy()
    countries ={}
    for k,v in locations.items():
            components = v[0]['address_components']
            country = next((component['long_name'] for component in components if 'country' in component['types']), 'Unknown')
            countries[k] = country
            
    dx['Country'] = dx['Institution'].map(countries)
    dx = dx[dx['Country'] != 'Unknown']
    dx['Code'] = dx['Country'].map(country_to_code)
    dx['Continent'] = dx['Code'].map(code_to_region_WHO)
    return dx.dropna()


def create_interactive_plot(country_data, disease):
    """
    Create an interactive scatter plot of prevalence vs research experiments for countries using Plotly
    
    Returns:
    --------
    plotly.graph_objs._figure.Figure
        The created interactive plot
    """
    # Create interactive scatter plot
    fig = px.scatter(
        country_data, 
        x='Cases per 1000', 
        y='count',
        hover_name='Code',
        hover_data={
            'Cases per 1000': ':.2f',
            'count': ':.2f'
        },
        title=f'Research Experiments vs {disease} Prevalence by Country',
        labels={
            'Cases per 1000': f'{disease} Prevalence (Normalized)', 
            'count': 'Number of Research Experiments (Normalized)'
        },
        color_discrete_sequence=['blue']
    )
    
    # Customize the layout
    fig.update_layout(
        plot_bgcolor='white',
        height=800,
        width=1200,
    )
    
    # Add grid lines
    fig.update_xaxes(
        showgrid=True, 
        gridwidth=1, 
        gridcolor='lightgray'
    )
    fig.update_yaxes(
        showgrid=True, 
        gridwidth=1, 
        gridcolor='lightgray'
    )
    
    # Adjust marker sizes
    fig.update_traces(
        marker=dict(
            size=10,
            line=dict(width=1, color='DarkSlateGrey')
        )
    )
    
    return fig


def plot_prevalence_wrt_research(df_research) :
    pf_malaria_country_cases = country_cases_to_cases_per_1000(pf_malaria_data, country_codes_no_population_data, country_populations)

    polio_country_cases = country_cases_to_cases_per_1000(polio_data, country_codes_no_population_data, country_populations)

    hiv_country_cases = country_cases_to_cases_per_1000(hiv_data, country_codes_no_population_data, country_populations)

    pv_malaria_country_cases = country_cases_to_cases_per_1000(pv_malaria_data, country_codes_no_population_data, country_populations)

    tuberculosis_country_cases = country_cases_per_x_to_cases_per_1000(tuberculosis_data, 1000, country_codes_no_population_data, country_populations)

    escherichia_coli_country_cases = country_cases_per_x_to_cases_per_1000(escherichia_coli_data, 100, country_codes_no_population_data, country_populations)

    hepatitis_c_country_cases = country_cases_per_x_to_cases_per_1000(hepatitis_c_data, 100, country_codes_no_population_data, country_populations)

    staphylococcus_aureus_country_cases = country_cases_per_x_to_cases_per_1000(staphylococcus_aureus_data, 100, country_codes_no_population_data, country_populations)

    #HIV

    df_HIV = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(HIV_sources)]
    research_counts_HIV_country = df_HIV["Code"].value_counts()
    country_data_HIV = pd.merge(
        hiv_country_cases, 
        research_counts_HIV_country, 
        on='Code', 
        how='outer'
    )
    country_data_HIV.fillna(0, inplace=True)

    country_data_HIV["Cases per 1000"] = (country_data_HIV["Cases per 1000"]-country_data_HIV["Cases per 1000"].min())/(country_data_HIV["Cases per 1000"].max() -country_data_HIV["Cases per 1000"].min())
    country_data_HIV["count"] = (country_data_HIV["count"]-country_data_HIV["count"].min())/(country_data_HIV["count"].max() -country_data_HIV["count"].min())

    fig = create_interactive_plot(country_data_HIV, 'HIV')
    fig.show()

    #Pf Malaria

    df_Plasmodium_falciparum = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Plasmodium_falciparum_sources)]
    research_counts_Plasmodium_falciparum_country = df_Plasmodium_falciparum["Code"].value_counts()
    country_data_Plasmodium_falciparum = pd.merge(
        pf_malaria_country_cases, 
        research_counts_Plasmodium_falciparum_country, 
        on='Code', 
        how='outer'
    )
    country_data_Plasmodium_falciparum.fillna(0, inplace=True)

    country_data_Plasmodium_falciparum["Cases per 1000"] = (country_data_Plasmodium_falciparum["Cases per 1000"]-country_data_Plasmodium_falciparum["Cases per 1000"].min())/(country_data_Plasmodium_falciparum["Cases per 1000"].max() -country_data_Plasmodium_falciparum["Cases per 1000"].min())
    country_data_Plasmodium_falciparum["count"] = (country_data_Plasmodium_falciparum["count"]-country_data_Plasmodium_falciparum["count"].min())/(country_data_Plasmodium_falciparum["count"].max() -country_data_Plasmodium_falciparum["count"].min())

    fig = create_interactive_plot(country_data_Plasmodium_falciparum, 'Plasmodium falciparum')
    fig.show()

    # Poliovirus

    df_Poliovirus = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Poliovirus_sources)]
    research_counts_Poliovirus_country = df_Poliovirus["Code"].value_counts()
    country_data_Poliovirus = pd.merge(
        polio_country_cases, 
        research_counts_Poliovirus_country, 
        on='Code', 
        how='outer'
    )
    country_data_Poliovirus.fillna(0, inplace=True)

    country_data_Poliovirus["Cases per 1000"] = (country_data_Poliovirus["Cases per 1000"]-country_data_Poliovirus["Cases per 1000"].min())/(country_data_Poliovirus["Cases per 1000"].max() -country_data_Poliovirus["Cases per 1000"].min())
    country_data_Poliovirus["count"] = (country_data_Poliovirus["count"]-country_data_Poliovirus["count"].min())/(country_data_Poliovirus["count"].max() -country_data_Poliovirus["count"].min())

    fig = create_interactive_plot(country_data_Poliovirus,'Poliovirus')
    fig.show()

    #Plasmodium vivax

    df_Plasmodium_vivax = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Plasmodium_vivax_sources)]
    research_counts_Plasmodium_vivax_country = df_Plasmodium_vivax["Code"].value_counts()
    country_data_Plasmodium_vivax = pd.merge(
        pv_malaria_country_cases, 
        research_counts_Plasmodium_vivax_country, 
        on='Code', 
        how='outer'
    )
    country_data_Plasmodium_vivax.fillna(0, inplace=True)

    country_data_Plasmodium_vivax["Cases per 1000"] = (country_data_Plasmodium_vivax["Cases per 1000"]-country_data_Plasmodium_vivax["Cases per 1000"].min())/(country_data_Plasmodium_vivax["Cases per 1000"].max() -country_data_Plasmodium_vivax["Cases per 1000"].min())
    country_data_Plasmodium_vivax["count"] = (country_data_Plasmodium_vivax["count"]-country_data_Plasmodium_vivax["count"].min())/(country_data_Plasmodium_vivax["count"].max() -country_data_Plasmodium_vivax["count"].min())

    fig = create_interactive_plot(country_data_Plasmodium_vivax, 'Plasmodium vivax')
    fig.show()


    #Tuberculosis

    df_Tuberculosis = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Tuberculosis_sources)]
    research_counts_Tuberculosis_country = df_Tuberculosis["Code"].value_counts()
    country_data_Tuberculosis = pd.merge(
        tuberculosis_country_cases, 
        research_counts_Tuberculosis_country, 
        on='Code', 
        how='outer'
    )
    country_data_Tuberculosis.fillna(0, inplace=True)

    country_data_Tuberculosis["Cases per 1000"] = (country_data_Tuberculosis["Cases per 1000"]-country_data_Tuberculosis["Cases per 1000"].min())/(country_data_Tuberculosis["Cases per 1000"].max() -country_data_Tuberculosis["Cases per 1000"].min())
    country_data_Tuberculosis["count"] = (country_data_Tuberculosis["count"]-country_data_Tuberculosis["count"].min())/(country_data_Tuberculosis["count"].max() -country_data_Tuberculosis["count"].min())

    fig = create_interactive_plot(country_data_Tuberculosis, 'Tuberculosis')
    fig.show()

    #Hepatitis C

    df_Hepatitis_C = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Hepatitis_C_sources)]
    research_counts_Hepatitis_C_country = df_Hepatitis_C["Code"].value_counts()
    country_data_Hepatitis_C = pd.merge(
        hepatitis_c_country_cases, 
        research_counts_Hepatitis_C_country, 
        on='Code', 
        how='outer'
    )
    country_data_Hepatitis_C.fillna(0, inplace=True)

    country_data_Hepatitis_C["Cases per 1000"] = (country_data_Hepatitis_C["Cases per 1000"]-country_data_Hepatitis_C["Cases per 1000"].min())/(country_data_Hepatitis_C["Cases per 1000"].max() -country_data_Hepatitis_C["Cases per 1000"].min())
    country_data_Hepatitis_C["count"] = (country_data_Hepatitis_C["count"]-country_data_Hepatitis_C["count"].min())/(country_data_Hepatitis_C["count"].max() -country_data_Hepatitis_C["count"].min())

    fig = create_interactive_plot(country_data_Hepatitis_C, 'Hepatitis C')
    fig.show()

    #Escherichia coli

    df_Escherichia_coli = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Escherichia_coli_sources)]
    research_counts_Escherichia_coli_country = df_Escherichia_coli["Code"].value_counts()
    country_data_Escherichia_coli = pd.merge(
        escherichia_coli_country_cases, 
        research_counts_Escherichia_coli_country, 
        on='Code', 
        how='outer'
    )
    country_data_Escherichia_coli.fillna(0, inplace=True)

    country_data_Escherichia_coli["Cases per 1000"] = (country_data_Escherichia_coli["Cases per 1000"]-country_data_Escherichia_coli["Cases per 1000"].min())/(country_data_Escherichia_coli["Cases per 1000"].max() -country_data_Escherichia_coli["Cases per 1000"].min())
    country_data_Escherichia_coli["count"] = (country_data_Escherichia_coli["count"]-country_data_Escherichia_coli["count"].min())/(country_data_Escherichia_coli["count"].max() -country_data_Escherichia_coli["count"].min())

    fig = create_interactive_plot(country_data_Escherichia_coli, 'Escherichia coli')
    fig.show()  

    #Staphylococcus aureus

    df_Staphylococcus_aureus = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Staphylococcus_aureus_sources)]
    research_counts_Staphylococcus_aureus_country = df_Staphylococcus_aureus["Code"].value_counts()
    country_data_Staphylococcus_aureus = pd.merge(
        staphylococcus_aureus_country_cases, 
        research_counts_Staphylococcus_aureus_country, 
        on='Code', 
        how='outer'
    )
    country_data_Staphylococcus_aureus.fillna(0, inplace=True)

    country_data_Staphylococcus_aureus["Cases per 1000"] = (country_data_Staphylococcus_aureus["Cases per 1000"]-country_data_Staphylococcus_aureus["Cases per 1000"].min())/(country_data_Staphylococcus_aureus["Cases per 1000"].max() -country_data_Staphylococcus_aureus["Cases per 1000"].min())
    country_data_Staphylococcus_aureus["count"] = (country_data_Staphylococcus_aureus["count"]-country_data_Staphylococcus_aureus["count"].min())/(country_data_Staphylococcus_aureus["count"].max() -country_data_Staphylococcus_aureus["count"].min())

    fig = create_interactive_plot(country_data_Staphylococcus_aureus, 'Staphylococcus_aureus')
    fig.show()  


def download_interactive_plots(df_research) :
    pf_malaria_country_cases = country_cases_to_cases_per_1000(pf_malaria_data, country_codes_no_population_data, country_populations)

    polio_country_cases = country_cases_to_cases_per_1000(polio_data, country_codes_no_population_data, country_populations)

    hiv_country_cases = country_cases_to_cases_per_1000(hiv_data, country_codes_no_population_data, country_populations)

    pv_malaria_country_cases = country_cases_to_cases_per_1000(pv_malaria_data, country_codes_no_population_data, country_populations)

    tuberculosis_country_cases = country_cases_per_x_to_cases_per_1000(tuberculosis_data, 1000, country_codes_no_population_data, country_populations)

    escherichia_coli_country_cases = country_cases_per_x_to_cases_per_1000(escherichia_coli_data, 100, country_codes_no_population_data, country_populations)

    hepatitis_c_country_cases = country_cases_per_x_to_cases_per_1000(hepatitis_c_data, 100, country_codes_no_population_data, country_populations)

    staphylococcus_aureus_country_cases = country_cases_per_x_to_cases_per_1000(staphylococcus_aureus_data, 100, country_codes_no_population_data, country_populations)

    #HIV
    print('HIV')

    df_HIV = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(HIV_sources)]
    research_counts_HIV_country = df_HIV["Code"].value_counts()
    country_data_HIV = pd.merge(
        hiv_country_cases, 
        research_counts_HIV_country, 
        on='Code', 
        how='outer'
    )
    country_data_HIV.fillna(0, inplace=True)

    country_data_HIV["Cases per 1000"] = (country_data_HIV["Cases per 1000"]-country_data_HIV["Cases per 1000"].min())/(country_data_HIV["Cases per 1000"].max() -country_data_HIV["Cases per 1000"].min())
    country_data_HIV["count"] = (country_data_HIV["count"]-country_data_HIV["count"].min())/(country_data_HIV["count"].max() -country_data_HIV["count"].min())

    fig = create_interactive_plot(country_data_HIV, 'HIV')
    pio.write_html(fig, "plot_HIV.html")

    #Pf Malaria

    print('Pf Malaria')

    df_Plasmodium_falciparum = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Plasmodium_falciparum_sources)]
    research_counts_Plasmodium_falciparum_country = df_Plasmodium_falciparum["Code"].value_counts()
    country_data_Plasmodium_falciparum = pd.merge(
        pf_malaria_country_cases, 
        research_counts_Plasmodium_falciparum_country, 
        on='Code', 
        how='outer'
    )
    country_data_Plasmodium_falciparum.fillna(0, inplace=True)

    country_data_Plasmodium_falciparum["Cases per 1000"] = (country_data_Plasmodium_falciparum["Cases per 1000"]-country_data_Plasmodium_falciparum["Cases per 1000"].min())/(country_data_Plasmodium_falciparum["Cases per 1000"].max() -country_data_Plasmodium_falciparum["Cases per 1000"].min())
    country_data_Plasmodium_falciparum["count"] = (country_data_Plasmodium_falciparum["count"]-country_data_Plasmodium_falciparum["count"].min())/(country_data_Plasmodium_falciparum["count"].max() -country_data_Plasmodium_falciparum["count"].min())

    fig = create_interactive_plot(country_data_Plasmodium_falciparum, 'Plasmodium falciparum')
    pio.write_html(fig, "plot_Plasmodium_falciparum.html")

    # Poliovirus

    print('Poliovirus')

    df_Poliovirus = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Poliovirus_sources)]
    research_counts_Poliovirus_country = df_Poliovirus["Code"].value_counts()
    country_data_Poliovirus = pd.merge(
        polio_country_cases, 
        research_counts_Poliovirus_country, 
        on='Code', 
        how='outer'
    )
    country_data_Poliovirus.fillna(0, inplace=True)

    country_data_Poliovirus["Cases per 1000"] = (country_data_Poliovirus["Cases per 1000"]-country_data_Poliovirus["Cases per 1000"].min())/(country_data_Poliovirus["Cases per 1000"].max() -country_data_Poliovirus["Cases per 1000"].min())
    country_data_Poliovirus["count"] = (country_data_Poliovirus["count"]-country_data_Poliovirus["count"].min())/(country_data_Poliovirus["count"].max() -country_data_Poliovirus["count"].min())

    fig = create_interactive_plot(country_data_Poliovirus,'Poliovirus')
    pio.write_html(fig, "plot_Poliovirus.html")

    #Plasmodium vivax

    print('Plasmodium vivax')

    df_Plasmodium_vivax = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Plasmodium_vivax_sources)]
    research_counts_Plasmodium_vivax_country = df_Plasmodium_vivax["Code"].value_counts()
    country_data_Plasmodium_vivax = pd.merge(
        pv_malaria_country_cases, 
        research_counts_Plasmodium_vivax_country, 
        on='Code', 
        how='outer'
    )
    country_data_Plasmodium_vivax.fillna(0, inplace=True)

    country_data_Plasmodium_vivax["Cases per 1000"] = (country_data_Plasmodium_vivax["Cases per 1000"]-country_data_Plasmodium_vivax["Cases per 1000"].min())/(country_data_Plasmodium_vivax["Cases per 1000"].max() -country_data_Plasmodium_vivax["Cases per 1000"].min())
    country_data_Plasmodium_vivax["count"] = (country_data_Plasmodium_vivax["count"]-country_data_Plasmodium_vivax["count"].min())/(country_data_Plasmodium_vivax["count"].max() -country_data_Plasmodium_vivax["count"].min())

    fig = create_interactive_plot(country_data_Plasmodium_vivax, 'Plasmodium vivax')
    pio.write_html(fig, "plot_Plasmodium_vivax.html")


    #Tuberculosis

    print('Tuberculosis')

    df_Tuberculosis = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Tuberculosis_sources)]
    research_counts_Tuberculosis_country = df_Tuberculosis["Code"].value_counts()
    country_data_Tuberculosis = pd.merge(
        tuberculosis_country_cases, 
        research_counts_Tuberculosis_country, 
        on='Code', 
        how='outer'
    )
    country_data_Tuberculosis.fillna(0, inplace=True)

    country_data_Tuberculosis["Cases per 1000"] = (country_data_Tuberculosis["Cases per 1000"]-country_data_Tuberculosis["Cases per 1000"].min())/(country_data_Tuberculosis["Cases per 1000"].max() -country_data_Tuberculosis["Cases per 1000"].min())
    country_data_Tuberculosis["count"] = (country_data_Tuberculosis["count"]-country_data_Tuberculosis["count"].min())/(country_data_Tuberculosis["count"].max() -country_data_Tuberculosis["count"].min())

    fig = create_interactive_plot(country_data_Tuberculosis, 'Tuberculosis')
    pio.write_html(fig, "plot_Tuberculosis.html")

    #Hepatitis C

    print('Hepatitis C')

    df_Hepatitis_C = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Hepatitis_C_sources)]
    research_counts_Hepatitis_C_country = df_Hepatitis_C["Code"].value_counts()
    country_data_Hepatitis_C = pd.merge(
        hepatitis_c_country_cases, 
        research_counts_Hepatitis_C_country, 
        on='Code', 
        how='outer'
    )
    country_data_Hepatitis_C.fillna(0, inplace=True)

    country_data_Hepatitis_C["Cases per 1000"] = (country_data_Hepatitis_C["Cases per 1000"]-country_data_Hepatitis_C["Cases per 1000"].min())/(country_data_Hepatitis_C["Cases per 1000"].max() -country_data_Hepatitis_C["Cases per 1000"].min())
    country_data_Hepatitis_C["count"] = (country_data_Hepatitis_C["count"]-country_data_Hepatitis_C["count"].min())/(country_data_Hepatitis_C["count"].max() -country_data_Hepatitis_C["count"].min())

    fig = create_interactive_plot(country_data_Hepatitis_C, 'Hepatitis C')
    pio.write_html(fig, "plot_Hepatitis_C.html")

    #Escherichia coli

    print('Escherichia coli')

    df_Escherichia_coli = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Escherichia_coli_sources)]
    research_counts_Escherichia_coli_country = df_Escherichia_coli["Code"].value_counts()
    country_data_Escherichia_coli = pd.merge(
        escherichia_coli_country_cases, 
        research_counts_Escherichia_coli_country, 
        on='Code', 
        how='outer'
    )
    country_data_Escherichia_coli.fillna(0, inplace=True)

    country_data_Escherichia_coli["Cases per 1000"] = (country_data_Escherichia_coli["Cases per 1000"]-country_data_Escherichia_coli["Cases per 1000"].min())/(country_data_Escherichia_coli["Cases per 1000"].max() -country_data_Escherichia_coli["Cases per 1000"].min())
    country_data_Escherichia_coli["count"] = (country_data_Escherichia_coli["count"]-country_data_Escherichia_coli["count"].min())/(country_data_Escherichia_coli["count"].max() -country_data_Escherichia_coli["count"].min())

    fig = create_interactive_plot(country_data_Escherichia_coli, 'Escherichia coli')
    pio.write_html(fig, "plot_Escherichia_coli.html") 

    #Staphylococcus aureus

    print('Staphylococcus aureus')

    df_Staphylococcus_aureus = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Staphylococcus_aureus_sources)]
    research_counts_Staphylococcus_aureus_country = df_Staphylococcus_aureus["Code"].value_counts()
    country_data_Staphylococcus_aureus = pd.merge(
        staphylococcus_aureus_country_cases, 
        research_counts_Staphylococcus_aureus_country, 
        on='Code', 
        how='outer'
    )
    country_data_Staphylococcus_aureus.fillna(0, inplace=True)

    country_data_Staphylococcus_aureus["Cases per 1000"] = (country_data_Staphylococcus_aureus["Cases per 1000"]-country_data_Staphylococcus_aureus["Cases per 1000"].min())/(country_data_Staphylococcus_aureus["Cases per 1000"].max() -country_data_Staphylococcus_aureus["Cases per 1000"].min())
    country_data_Staphylococcus_aureus["count"] = (country_data_Staphylococcus_aureus["count"]-country_data_Staphylococcus_aureus["count"].min())/(country_data_Staphylococcus_aureus["count"].max() -country_data_Staphylococcus_aureus["count"].min())

    fig = create_interactive_plot(country_data_Staphylococcus_aureus, 'Staphylococcus_aureus')
    pio.write_html(fig, "plot_Staphylococcus_aureus.html")  

def df_regional_prevalence_for_each_disease():
    pf_malaria_region_cases = country_cases_to_region_cases_per_1000(pf_malaria_data, region_population)
    hiv_region_cases = country_cases_to_region_cases_per_1000(hiv_data, region_population)
    polio_region_cases = country_cases_to_region_cases_per_1000(polio_data, region_population)
    pv_malaria_region_cases = country_cases_to_region_cases_per_1000(pv_malaria_data, region_population)
    tuberculosis_region_cases = country_cases_to_region_cases_per_1000(tuberculosis_data, region_population)
    escherichia_coli_region_cases = country_cases_to_region_cases_per_1000(escherichia_coli_data, region_population)
    hepatitis_c_region_cases = country_cases_to_region_cases_per_1000(hepatitis_c_data, region_population)
    staphylococcus_aureus_region_cases = country_cases_to_region_cases_per_1000(staphylococcus_aureus_data, region_population)


    all_region_cases = [pf_malaria_region_cases, hiv_region_cases, polio_region_cases, pv_malaria_region_cases, tuberculosis_region_cases, escherichia_coli_region_cases, hepatitis_c_region_cases, staphylococcus_aureus_region_cases]

    df = pd.DataFrame()
    
    extracted_cases_per_1000 = [df['Cases per 1000'] for df in all_region_cases]
    df = pd.DataFrame(extracted_cases_per_1000, index=disease_names)

    for i in range(len(df.columns)):
        df.rename(columns={i: region_codes[i]}, inplace=True)

    return df

def get_research_per_region(disease, df, sources):
    df['country_code'] = df['Country'].map(country_to_code)
    df['region_code'] = df['country_code'].map(code_to_region_WHO)

    df = df[df["Target Source Organism According to Curator or DataSource"].isin(sources)]
    value_counts_region = df["region_code"].value_counts()
    return value_counts_region

def get_research_per_country(disease, df, sources):
    df['country_code'] = df['Country'].map(country_to_code)
    df = df[df["Target Source Organism According to Curator or DataSource"].isin(sources)]
    value_counts_country = df["country_code"].value_counts()
    return value_counts_country

### OLS functions

def normalize_rows_of_df(df):
    mean = df.mean(axis=1)
    std = df.std(axis=1)
    df = df.apply(lambda x: (x - mean) / std, axis=0)
    return df

def research_of_diseases_per_country_df(df_research):
    hiv = get_research_per_country('HIV',df_research, HIV_sources)
    tuberculosis = get_research_per_country('Tuberculosis',df_research, Tuberculosis_sources)
    hepatitis = get_research_per_country('Hepatitis',df_research, Hepatitis_C_sources)
    plasmodium_falciparum = get_research_per_country('Plasmodium falciparum',df_research, Plasmodium_falciparum_sources)
    poliovirus = get_research_per_country('Poliovirus',df_research, Poliovirus_sources)
    plasmodium_vivax = get_research_per_country('Plasmodium vivax',df_research, Plasmodium_vivax_sources)
    escherichia_coli = get_research_per_country('Escherichia coli',df_research, Escherichia_coli_sources)
    staphylococcus_aureus = get_research_per_country('Staphylococcus aureus',df_research,Staphylococcus_aureus_sources)

    diseases_data = [plasmodium_falciparum, hiv, poliovirus, plasmodium_vivax, tuberculosis, escherichia_coli, hepatitis, staphylococcus_aureus]
    df_country = pd.DataFrame()

    for disease in disease_names:
        df_country[disease] = diseases_data[disease_names.index(disease)]
    df_country.fillna(0, inplace=True)
    df_country = df_country.transpose()
    df_country = normalize_rows_of_df(df_country)
    return df_country

def research_of_diseases_per_region_df(df_research):
    hiv = get_research_per_region('HIV',df_research, HIV_sources)
    tuberculosis = get_research_per_region('Tuberculosis',df_research, Tuberculosis_sources)
    hepatitis = get_research_per_region('Hepatitis',df_research, Hepatitis_C_sources)
    plasmodium_falciparum = get_research_per_region('Plasmodium falciparum',df_research, Plasmodium_falciparum_sources)
    poliovirus = get_research_per_region('Poliovirus',df_research, Poliovirus_sources)
    plasmodium_vivax = get_research_per_region('Plasmodium vivax',df_research, Plasmodium_vivax_sources)
    escherichia_coli = get_research_per_region('Escherichia coli',df_research, Escherichia_coli_sources)
    staphylococcus_aureus = get_research_per_region('Staphylococcus aureus',df_research,Staphylococcus_aureus_sources)

    diseases_data = [plasmodium_falciparum, hiv, poliovirus, plasmodium_vivax, tuberculosis, escherichia_coli, hepatitis, staphylococcus_aureus]
    region_codes = ['AFR', 'AMR', 'EMR', 'EUR', 'SEAR', 'WPR']
    df = pd.DataFrame(columns = region_codes)
    for disease in disease_names:
        df.loc[disease] = diseases_data[disease_names.index(disease)]
    df.fillna(0, inplace=True)
    df = normalize_rows_of_df(df)
    return df

def concatenate_disease_and_research_dfs(diseases_df, research_df):
    diseases_df_suffix = diseases_df.add_suffix('_disease')
    research_df_suffix = research_df.add_suffix('_research')
    df = pd.concat([diseases_df_suffix, research_df_suffix], axis=1)
    return df

def region_region_smf_model(df):

    region_codes_disease = 'AFR_disease + AMR_disease + EMR_disease + EUR_disease + SEAR_disease + WPR_disease'
    region_codes_research = ['AFR_research', 'AMR_research', 'EMR_research', 'EUR_research', 'SEAR_research', 'WPR_research']
    coef_df = pd.DataFrame()
    p_value_df = pd.DataFrame()
    for region_research in region_codes_research:
        region_mod = smf.ols(formula=f'{region_research} ~ {region_codes_disease}', data=df, missing='drop')
        region_res = region_mod.fit()
        region_coef = region_res.params
        region_coef = region_coef[1:]
        region_p_values = region_res.pvalues[1:]
        coef_df[region_research] = region_coef
        p_value_df[region_research] = region_p_values
      
    return coef_df, p_value_df

def region_region_heatmap(coef_df, p_value_df):

    fig, ax = plt.subplots(figsize=(8, 6))
    coef_df = coef_df.transpose()
    p_value_df = p_value_df.transpose()
    # Plot the heatmap from coef_df
    sns.heatmap(coef_df, 
                annot=True, 
                fmt=".2f", 
                cmap="coolwarm", 
                cbar=True, 
                linewidths=0.5, 
                linecolor='gray', 
                ax=ax)

    # Highlight cells where p_value > 0.95 with thicker edges
    for i in range(p_value_df.shape[0]):
        for j in range(p_value_df.shape[1]):
            if p_value_df.iloc[i, j] < 0.1:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=3))

    # Adjust axis labels
    ax.set_xlabel("Rows")
    ax.set_ylabel("Columns")
    ax.set_xticklabels(coef_df.columns)
    ax.set_yticklabels(coef_df.index, rotation=0)

def region_country_smf_model(df, research_column_names):

    region_codes_disease = 'AFR_disease + AMR_disease + EMR_disease + EUR_disease + SEAR_disease + WPR_disease'
    coef_df = pd.DataFrame()
    p_value_df = pd.DataFrame()
    for country_research in research_column_names:
        region_mod = smf.ols(formula=f'{country_research} ~ {region_codes_disease}', data=df, missing='drop')
        region_res = region_mod.fit()
        region_coef = region_res.params
        region_coef = region_coef[1:]
        region_p_values = region_res.pvalues[1:]
        coef_df[country_research] = region_coef
        p_value_df[country_research] = region_p_values

    region_codes = ['AFR', 'AMR', 'EMR', 'EUR', 'SEAR', 'WPR']
    total_region_cases = np.zeros(len(region_codes))
    significant_research_in_same_region = np.zeros(len(region_codes))

    return coef_df, p_value_df

def region_country_heatmap(coef_df, p_value_df):
    coef_df = coef_df.transpose()
    p_value_df = p_value_df.transpose()
    fig, ax = plt.subplots(figsize=(8, 6))
    # Plot the heatmap from coef_df
    sns.heatmap(coef_df, 
                annot=True, 
                fmt=".2f", 
                cmap="coolwarm", 
                cbar=True, 
                linewidths=0.5, 
                linecolor='gray', 
                ax=ax)

    # Highlight cells where p_value > 0.95 with thicker edges
    for i in range(p_value_df.shape[0]):
        for j in range(p_value_df.shape[1]):
            if p_value_df.iloc[i, j] < 0.1:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=3))

    # Adjust axis labels
    ax.set_xlabel("Rows")
    ax.set_ylabel("Columns")
    ax.set_xticklabels(coef_df.columns)
    ax.set_yticklabels(coef_df.index, rotation=0)