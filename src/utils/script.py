# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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


def plot_disease(disease, df, sources) :  

    df = df[df["Target Source Organism According to Curator or DataSource"].isin(sources)]
    value_counts_country = df["Code"].value_counts().head(10)
    value_counts_continent = df["Continent"].value_counts().head(10)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # Plot for countries
    axes[0].bar(value_counts_country.index, value_counts_country.values, color='skyblue')
    axes[0].set_xticks(range(len(value_counts_country.index)))  # Set tick positions explicitly
    axes[0].set_xticklabels(value_counts_country.index, rotation=-45, ha='left')  # Set tick labels
    axes[0].set_xlabel('Country')
    axes[0].set_ylabel(f'Number of experiments on {disease}')
    axes[0].set_title(f'Number of experiments on {disease} w.r.t. Country')

    # Plot for continents
    axes[1].bar(value_counts_continent.index, value_counts_continent.values)
    axes[1].set_xticks(range(len(value_counts_continent.index)))  # Set tick positions explicitly
    axes[1].set_xticklabels(value_counts_continent.index, rotation=-45, ha='left')  # Set tick labels
    axes[1].set_xlabel('Continent')
    axes[1].set_ylabel(f'Number of experiments on {disease}')
    axes[1].set_title(f'Number of experiments on {disease} w.r.t. Continent')

    plt.tight_layout()
    plt.show()

def heatmap_research(df_research) :
    
    regions = ['AFR', 'AMR', 'EMR', 'EUR', 'SEAR', 'WPR']

    HIV_sources = ['Human immunodeficiency virus 1', 'Human immunodeficiency virus', 'Human immunodeficiency virus type 1 group M subtype B (isolate BRU/LAI)', 'Human immunodeficiency virus type 1 group M subtype B (isolate YU-2)', 'Human immunodeficiency virus type 1 group M subtype B (isolate HXB2)', 'Human immunodeficiency virus 2', 'Human immunodeficiency virus type 1 group M subtype B (isolate PCV12)', 'Human immunodeficiency virus type 1 group M subtype B (isolate MN)']
    df_HIV = df_HIV = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(HIV_sources)]

    Plasmodium_falciparum_sources = ['Plasmodium falciparum', 'Plasmodium falciparum (isolate 3D7)', 'Plasmodium falciparum (isolate K1 / Thailand)', 'Plasmodium falciparum (isolate FcB1 / Columbia)']
    df_Plasmodium_falciparum = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Plasmodium_falciparum_sources)]

    Poliovirus_sources = ['Poliovirus type 1 (strain Mahoney)']
    df_Poliovirus = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Poliovirus_sources)]

    Plasmodium_vivax_sources = ['Plasmodium vivax']
    df_Plasmodium_vivax = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Plasmodium_vivax_sources)]

    Tuberculosis_sources = ['Mycobacterium tuberculosis', 'Mycobacterium tuberculosis H37Rv']
    df_Tuberculosis = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Tuberculosis_sources)]

    Hepatitis_C_sources = ['Hepatitis C virus', 'Hepatitis C virus genotype 1a (isolate H)', 'Hepatitis C virus genotype 1b (isolate BK)', 'Hepatitis C virus genotype 3a (isolate NZL1)', 'Hepatitis C virus genotype 2b (isolate HC-J8)', 'Hepatitis C virus genotype 1b (isolate Taiwan)', 'Hepatitis C virus genotype 2', 'Hepatitis C virus genotype 4a (isolate ED43)', 'Hepatitis C virus genotype 6a (isolate EUHK2)']
    df_Hepatitis_C = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Hepatitis_C_sources)]

    Escherichia_coli_sources = ['Escherichia coli', 'Escherichia coli str. K-12 substr. MG1655', 'Escherichia coli (strain UTI89 / UPEC)', 'Escherichia coli (strain K12)', 'Escherichia coli O157:H7', 'Escherichia coli O6']
    df_Escherichia_coli = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Escherichia_coli_sources)]

    Staphylococcus_aureus_sources = ['Staphylococcus aureus', 'Staphylococcus aureus (strain MRSA252)', 'Staphylococcus aureus (strain Mu50 / ATCC 700699)', 'Staphylococcus aureus (strain MW2)', 'Staphylococcus aureus (strain NCTC 8325)']
    df_Staphylococcus_aureus = df_research[df_research["Target Source Organism According to Curator or DataSource"].isin(Staphylococcus_aureus_sources)]

    HIV_research_counts = dict(df_HIV['Continent'].value_counts())
    Plasmodium_falciparum_research_counts = dict(df_Plasmodium_falciparum['Continent'].value_counts())
    Poliovirus_research_counts = dict(df_Poliovirus['Continent'].value_counts())
    Plasmodium_vivax_research_counts = dict(df_Plasmodium_vivax['Continent'].value_counts())
    Tuberculosis_research_counts = dict(df_Tuberculosis['Continent'].value_counts())
    Hepatitis_C_research_counts = dict(df_Hepatitis_C['Continent'].value_counts())
    Escherichia_coli_research_counts = dict(df_Escherichia_coli['Continent'].value_counts())
    Staphylococcus_aureus_research_counts = dict(df_Staphylococcus_aureus['Continent'].value_counts())

    for region in regions :
        if region not in HIV_research_counts :
            HIV_research_counts[region] = 0
        if region not in Plasmodium_falciparum_research_counts :
            Plasmodium_falciparum_research_counts[region] = 0
        if region not in Poliovirus_research_counts :
            Poliovirus_research_counts[region] = 0
        if region not in Plasmodium_vivax_research_counts :
            Plasmodium_vivax_research_counts[region] = 0
        if region not in Tuberculosis_research_counts :
            Tuberculosis_research_counts[region] = 0
        if region not in Hepatitis_C_research_counts :
            Hepatitis_C_research_counts[region] = 0
        if region not in Escherichia_coli_research_counts :
            Escherichia_coli_research_counts[region] = 0
        if region not in Staphylococcus_aureus_research_counts :
            Staphylococcus_aureus_research_counts[region] = 0

    

# Assuming the data dictionary is already created
    data = {
    'HIV': HIV_research_counts,
    'Plasmodium_falciparum': Plasmodium_falciparum_research_counts,
    'Poliovirus': Poliovirus_research_counts,
    'Plasmodium_vivax': Plasmodium_vivax_research_counts,
    'Tuberculosis': Tuberculosis_research_counts,
    'Hepatitis_C': Hepatitis_C_research_counts,
    'Escherichia_coli': Escherichia_coli_research_counts,
    'Staphylococcus_aureus': Staphylococcus_aureus_research_counts
    }

    # Convert to DataFrame 
    df_heatmap = pd.DataFrame(data)

    df_normalized = (df_heatmap - df_heatmap.min()) / (df_heatmap.max() - df_heatmap.min())
    df_standardized = (df_heatmap - df_heatmap.mean()) / df_heatmap.std()

    # Transpose the DataFrame to invert axess
    df_standardized_transposed = df_standardized.T
    df_normalized_transposed = df_normalized.T
    df_standardized_transposed = df_standardized_transposed[regions]
    df_normalized_transposed = df_normalized_transposed[regions]

    # Plot the heatmap with the 'RdBu' color palette
    plt.figure(figsize=(12, 8))
    sns.heatmap(df_normalized_transposed, annot=False, cmap='coolwarm', fmt='g')

    # Set plot title and labels
    plt.title('Number of Research Studies by Disease and Continent')
    plt.xlabel('Continent')
    plt.ylabel('Disease')

    # Ensure the x-axis labels are horizontal
    plt.xticks(rotation=0)

    # Display the heatmap
    plt.show()
