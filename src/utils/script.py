# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt

HIV_sources = ['Human immunodeficiency virus 1', 'Human immunodeficiency virus', 'Human immunodeficiency virus type 1 group M subtype B (isolate BRU/LAI)', 'Human immunodeficiency virus type 1 group M subtype B (isolate YU-2)', 'Human immunodeficiency virus type 1 group M subtype B (isolate HXB2)', 'Human immunodeficiency virus 2', 'Human immunodeficiency virus type 1 group M subtype B (isolate PCV12)', 'Human immunodeficiency virus type 1 group M subtype B (isolate MN)']

Plasmodium_falciparum_sources = ['Plasmodium falciparum', 'Plasmodium falciparum (isolate 3D7)', 'Plasmodium falciparum (isolate K1 / Thailand)', 'Plasmodium falciparum (isolate FcB1 / Columbia)']

Poliovirus_sources = ['Poliovirus type 1 (strain Mahoney)']

Plasmodium_vivax_sources = ['Plasmodium vivax']

Tuberculosis_sources = ['Mycobacterium tuberculosis', 'Mycobacterium tuberculosis H37Rv']

Hepatitis_C_sources = ['Hepatitis C virus', 'Hepatitis C virus genotype 1a (isolate H)', 'Hepatitis C virus genotype 1b (isolate BK)', 'Hepatitis C virus genotype 3a (isolate NZL1)', 'Hepatitis C virus genotype 2b (isolate HC-J8)', 'Hepatitis C virus genotype 1b (isolate Taiwan)', 'Hepatitis C virus genotype 2', 'Hepatitis C virus genotype 4a (isolate ED43)', 'Hepatitis C virus genotype 6a (isolate EUHK2)']

Escherichia_coli_sources = ['Escherichia coli', 'Escherichia coli str. K-12 substr. MG1655', 'Escherichia coli (strain UTI89 / UPEC)', 'Escherichia coli (strain K12)', 'Escherichia coli O157:H7', 'Escherichia coli O6']

Staphylococcus_aureus_sources = ['Staphylococcus aureus', 'Staphylococcus aureus (strain MRSA252)', 'Staphylococcus aureus (strain Mu50 / ATCC 700699)', 'Staphylococcus aureus (strain MW2)', 'Staphylococcus aureus (strain NCTC 8325)']

def read_tsv(file_path) :
    'Returns the data in path as a dataframe'
    return pd.read_csv(file_path, sep='\t',encoding='utf-8',on_bad_lines='skip')

def read_pickle(file_path) :
    return pd.read_pickle(file_path)

def describe(df) :
    "Prints the head, the shape the columns of the dataset"
    print("The 5 first rows of the data:")
    print(df.head())
    print('')
    print(f'The shape of the data: {df.shape}')
    print('')
    print('The columns in the data:')
    for col in df.columns :
        print(col)
    print('')
    print(f'The number of columns in the data: {len(df.columns)}')

def ligands_counts(df) :
    'Returns the value counts in the column "BindingDB Ligand Name"'
    return df["BindingDB Ligand Name"].value_counts()

def target_counts(df) : 
    'Returns the value counts in the column "Target Name"'
    return df["Target Name"].value_counts()

def target_source_counts(df) :
    'Returns the value counts in the column "Target Source Organism According to Curator or DataSource"'
    return df["Target Source Organism According to Curator or DataSource"].value_counts()
    
def df_with_target_sources(df) :
    "Returns the data where we know the source of the target"
    df_with_target_sources = df.dropna(subset = ["Target Source Organism According to Curator or DataSource"])
    return df_with_target_sources

def targets_with_and_without_source(df) :
    "Returns two lists, one of targets that have a source, and one of targets that don't have a source"
    targets_with_source = df_with_target_sources(df)["Target Name"].value_counts().index.tolist()
    targets = target_counts(df).index.tolist()
    targets_without_sources = list()
    for target in targets :
        if target not in targets_with_source :
            targets_without_sources.append(target)
    return (targets_with_source,targets_without_sources)

def sources(df) :
    "Returns the list of sources in the dataset"
    return df["Target Source Organism According to Curator or DataSource"].value_counts().index.tolist()
    
def targets_for_sources(df) :
    "Return a dictionary where each source is mapped to a list of targets for which it is a source"
    df_sources =  df_with_target_sources(df)
    source_targets_dict = dict()
    for source in sources(df) :
        source_targets_dict[source] = df_sources[df_sources["Target Source Organism According to Curator or DataSource"]==source]["Target Name"].value_counts().index.tolist()
    return source_targets_dict

def targets_with_multiple_sources(targets_with_source,targets_for_sources) :
    "Returns True if there are targets that have multiple sources, False otherwise"
    for target in targets_with_source :
        num_sources = 0
        for targets in targets_for_sources.values() :
            if target in targets :
                num_sources += 1
        if num_sources>1 :
            print("There are targets with multiple sources")
            return 
    print("There aren't targets with multiple sources")

def plot_proportions_of_experiments_given_TSO(df) :
    "Plots a pie chart of the proportions of experiments for each TSO from df with target sources"
    source_research_num = df_with_target_sources(df)["Target Source Organism According to Curator or DataSource"].value_counts()
    df_sources_reasearch = pd.DataFrame(list(source_research_num.items()), columns=['Source', 'Research count'])
    sources_reasearch_counts_chart = dict()
    for index, row in df_sources_reasearch.head(5).iterrows() :
        sources_reasearch_counts_chart[row["Source"]] = row["Research count"]
    sources_reasearch_counts_chart["Other"] = df_sources_reasearch.iloc[5:]["Research count"].sum() 
    # Create a pie chart
    plt.figure(figsize=(6, 6))  # Optional: set the figure size
    plt.pie(sources_reasearch_counts_chart.values(), labels=sources_reasearch_counts_chart.keys(), autopct='%1.1f%%', startangle=90)
    # Show the plot
    plt.title("Proportions of experiments given TSO")
    plt.show()

def plot_Number_of_experiments_w_r_t_sources(df) :
    plt.bar(['HIV','Tuberculosis','Hepatitis','Plasmodium falciparum','Poliovirus','Staphylococcus aureus','Plasmodium vivax','Escherichia coli'],

        [len(df[df["Target Source Organism According to Curator or DataSource"].isin(HIV_sources)]),
         len(df[df["Target Source Organism According to Curator or DataSource"].isin(Tuberculosis_sources)]),
         len(df[df["Target Source Organism According to Curator or DataSource"].isin(Hepatitis_C_sources)]),
         len(df[df["Target Source Organism According to Curator or DataSource"].isin(Plasmodium_falciparum_sources)]),
         len(df[df["Target Source Organism According to Curator or DataSource"].isin(Poliovirus_sources)]),
         len(df[df["Target Source Organism According to Curator or DataSource"].isin(Staphylococcus_aureus_sources)]),
         len(df[df["Target Source Organism According to Curator or DataSource"].isin(Plasmodium_vivax_sources)]),
        len(df[df["Target Source Organism According to Curator or DataSource"].isin(Escherichia_coli_sources)])]
    )
    plt.xticks(rotation=-45,ha='left')
    plt.xlabel('Source')
    plt.ylabel('Number of experiments')
    plt.title('Number of experiments w.r.t. sources')
    plt.show()
    
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
    'Egypt': 'EGY',
    'New Zealand': 'NZL',
    'Lithuania': 'LTU',
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
    'Iran': 'IRN',
    'Saudi Arabia': 'SAU',
    'Portugal': 'PRT',
    'Malaysia': 'MYS',
    'Thailand': 'THA',
    'Norway': 'NOR',
    'Hungary': 'HUN',
    'Bangladesh': 'BGD',
    'Mexico': 'MEX',
    'Argentina': 'ARG',
    'Jordan': 'JOR',
    'United Arab Emirates': 'ARE',
    'Slovakia': 'SVK',
    'Ireland': 'IRL',
    'Romania': 'ROU',
    'Chile': 'CHL',
    'Serbia': 'SRB',
    'Philippines': 'PHL',
    'Estonia': 'EST',
    'Ukraine': 'UKR',
    'Monaco': 'MCO',
    'Colombia': 'COL',
    'Croatia': 'HRV',
    'Morocco': 'MAR',
    'Tunisia': 'TUN',
    'Indonesia': 'IDN',
    'Latvia': 'LVA',
    'Vietnam': 'VNM',
    'Ecuador': 'ECU',
    'Bulgaria': 'BGR',
    'Zambia': 'ZMB',
    'Venezuela': 'VEN',
    'Lebanon': 'LBN',
    'Rwanda': 'RWA',
    'Oman': 'OMN',
    'Luxembourg': 'LUX',
    'Cyprus': 'CYP',
    'Uruguay': 'URY',
    'Botswana': 'BWA',
    'Algeria': 'DZA',
    'Grenada': 'GRD',
    'Jamaica': 'JAM',
    'Kuwait': 'KWT',
    'Sudan': 'SDN',
    'Peru': 'PER',
    'Qatar': 'QAT',
    'Zimbabwe': 'ZWE',
    'Bosnia and Herzegovina': 'BIH',
    'Sri Lanka': 'LKA',
    'Nepal': 'NPL',
    'Armenia': 'ARM',
    'Kenya': 'KEN',
    'Cuba': 'CUB',
    'Ethiopia': 'ETH',
    'Nigeria': 'NGA'
}

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


def plot_disease(disease, df, sources) :  

    df = df[df["Target Source Organism According to Curator or DataSource"].isin(sources)]
    value_counts_country = df["Country"].value_counts().head(10)
    value_counts_continent = df["Continent"].value_counts().head(10)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    axes[0].bar(value_counts_country.index, value_counts_country.values, color='skyblue')
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=-45, ha='left')
    axes[0].set_xlabel('Country')
    axes[0].set_ylabel(f'Number of experiments on {disease}')
    axes[0].set_title(f'Number of experiments on {disease} w.r.t. Country')

    axes[1].bar(value_counts_continent.index, value_counts_continent.values)
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=-45, ha='left')
    axes[1].set_xlabel('Continent')
    axes[1].set_ylabel(f'Number of experiments on {disease}')
    axes[1].set_title(f'Number of experiments on {disease} w.r.t. Countinent')

    plt.tight_layout()
    plt.show()