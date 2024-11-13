# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt

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
            return True
    return False

def plot_proportions_of_experiments_given_TSO(df) :
    "Plots a pie chart of the proportions of experiments for each TSO from df with target sources"
    source_research_num = df_with_target_sources(df)["Target Source Organism According to Curator or DataSource"].value_counts()
    df_sources_reasearch = pd.DataFrame(list(source_research_num.items()), columns=['Source', 'Research count'])
    sources_reasearch_counts_chart = dict()
    for index, row in df_sources_reasearch.head(10).iterrows() :
        sources_reasearch_counts_chart[row["Source"]] = row["Research count"]
    sources_reasearch_counts_chart["Other"] = df_sources_reasearch.iloc[10:]["Research count"].sum() 
    # Create a pie chart
    plt.figure(figsize=(6, 6))  # Optional: set the figure size
    plt.pie(sources_reasearch_counts_chart.values(), labels=sources_reasearch_counts_chart.keys(), autopct='%1.1f%%', startangle=90)
    # Show the plot
    plt.title("Proportions of experiments given TSO")
    plt.show()