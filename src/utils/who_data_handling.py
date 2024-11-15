import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

disease_names = ["Plasmodium Falciparum", "Human Immunodeficiency Virus",
                  "Poliovirus", "Plasmodium vivax", 
                  "Mycobacterium tuberculosis", "Escherichia coli", "Hepatitis C"]

disease_indicators_names = ["Number of indigenous P. falciparum malaria cases",
                            "Number of new HIV infections",
                            "Poliomyelitis - number of reported cases",
                            "Number of indigenous P. vivax malaria cases",
                            "Incidence of tuberculosis (per 100 000 population per year)",
                            "Proportion of bloodstream infection due to Escherichia coli resistant to third-generation cephalosporins (%)",
                            "Prevalence of chronic hepatitis C (HCV) ) in the general population (%)"
                            ]

disease_indicator_codes = ["MALARIA_PF_INDIG", "HIV_0000000026", "WHS3_49", "MALARIA_PV_INDIG",
                           "MDG_0000000020", "AMR_INFECT_ECOLI", "HEPATITIS_HCV_PREVALENCE_PER100"]

country_dimensions = pd.read_pickle("who_country_region_codes.pkl")
country_populations = pd.read_pickle("who_country_populations.pkl")

country_codes_no_population_data = country_dimensions[~country_dimensions['Code'].isin(country_populations['SpatialDim'])]['Code']
region_population = country_populations.groupby('ParentLocationCode')['NumericValue'].sum()

diseases_data = pd.read_pickle("who_diseases_data.pkl")

pf_malaria_data = diseases_data['pf_malaria']
hiv_data = diseases_data['hiv']
polio_data = diseases_data['polio']
pv_malaria_data = diseases_data['pv_malaria']
tuberculosis_data = diseases_data['tuberculosis']
escherichia_coli_data = diseases_data['escherichia_coli']
hepatitis_c_data = diseases_data['hepatitis_c']

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

all_region_cases = [pf_malaria_region_cases, hiv_region_cases, polio_region_cases, pv_malaria_region_cases, tuberculosis_region_cases, escherichia_coli_region_cases, hepatitis_c_region_cases]

for i in range(len(all_region_cases)):
    all_region_cases[i] = all_region_cases[i].rename(columns={'Cases per 1000': f'{disease_names[i]}'})

all_region_cases = pd.concat(all_region_cases, axis=1)
all_region_cases.drop(columns=['Cases'], inplace=True)

code_column = all_region_cases['Code']
first_code_column = code_column.iloc[:, 0]
all_region_cases.drop(columns=['Code'], inplace=True)

scaler = StandardScaler()
all_region_cases = pd.DataFrame(scaler.fit_transform(all_region_cases), columns=all_region_cases.columns)

all_region_cases.insert(0, 'Cases', first_code_column)
print(all_region_cases)

heatmap_data = all_region_cases.set_index('Cases').T
plt.figure(figsize=(10, 7))
sns.heatmap(heatmap_data, fmt=".2f", cmap='coolwarm')
plt.show()