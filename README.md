# Predicting Drug-Drug Interactions Based on Likeness and Geographical Prevalence

### Website: [https://ekaterina-borisova.github.io/ada_website/](https://ekaterina-borisova.github.io/ada_website/)

## Abstract 
What impact does geography have on the development of drugs? Specifically, does drug's origin matter in drug-drug interaction? The development of pharmaceuticals is a collaborative effort that often reflects the expertise and innovation of the institutions behind them. This project aims to investigate whether drugs developed within the same institutions exhibit superior interaction profiles compared to those developed independently. One may hypothesise that shared research environments, methodologies, and collaborative networks contribute to a higher likelihood of favorable interactions among drugs originating from the same institution. 

Does the prevalence of a disease in a region have an effect on local drug development? The world is highly interconnected and an outbreak of any disease anywhere is today a menace for the whole world. Is there however any trends one can observe in the location of drug research with respect to the geographical prevalence of a disease?

## Research Questions 
1. **Institutional Influence on DDIs**: Do drugs developed in the same universities interact positively?
2. **Regional Disease Influence on Research**: Does the prevalence of a disease in a region drive the amount and focus of research conducted there?  
3. **Compatibility of Popular Drugs**: Do highly "likable" drugs, based on QED, exhibit good interactions?
4. **Predictive DDI Modeling**: Can we predict interaction outcomes for novel drug pairs, and propose new candidates for future research? 


## Proposed Additional Datasets
- **DrugBank**: To cross-reference known DDIs and drug classification. Use the data on drug-interactions to classify based on description. Use the `DrugBank ID of Ligand` to map the two databases. From our initial data preprocessing we saw that the main DrugBank database provides only the data on `antagonistic` interaction of drugs which can be then sorted by severity. We reached out to the support team of DrugBank to get access if any to third-party databases which could provide more insights on the DDI.
- **WHO Database**: To find data on diseases researched in experiments in the BindingDB dataset, we first made a list of all unique words used for target sources in BindingDB. Then we searched WHO's databases for these same words. After doing this we identified 7 diseases for which BindingDB has exepriment data for and WHO has data on their prevalences worldwide. We then used the prevalences data to calculate the number of people affected by each disease both in each country and in each region as defined by WHO.
- **Geolocation data**: Using the google API, we try to link institutions to geographical coordinates and from these coordinates to countries and continents.  

## Methods

- Data cleaning: removing unusable samples in the data, for example: samples for which we don't have any  
- Data augmentation

**Analysing Institutions** 
1. Map institutions to their respective countries and continents
2. Investigate global research patterns and trends across regions
3. Compare the volume of research conducted by countries and continents
4. Compare research efforts on major diseases with their geographical prevalence, analyzing how research activity aligns with disease presence across countries and continents

**Analysing Likeness of Drug**
1. Extract from the SMILES all relevant properties of the drug: molecular weight, Log P, hydrogen bond donors, acceptors, number of rotatable bonds and polar surface area.
2. Apply Lipinski’s Rule of Five, QED scores, and other criteria (e.g., Veber’s Rule) to filter compounds for basic drug-likeness.
3. Focus the project on the most ‘likable’ drugs identified by the intersection of criteria results

**Drug-Drug Interaction**
1. Extract a more extensive database on DDI with examples of `synergetic` interactions.
2. Extract data on DDI from DrugBank and based on description cluster the interaction.
   1) First embed the description possibly with BioBert 
   2) Cluster the data using K-Means into 3 clusters based on severity `Major`, `Moderate` and `Minor`
3. Visualize the drug interaction based on institution location:
	By analyzing DDI data from BindingDB and enriching it with institutional information, we will assess patterns and correlations in drug compatibility across origins.
4. Try predicting the interaction by embeding smiles and training on DrugBank labeled data.
5. Use the developped model to predict the DDI of `likable` drugs in specific geographical locations


## Proposed Timeline
| Milestone         | Task Description                         | Deadline   |
|-------------------|------------------------------------------|------------|
| **Week 1**        | Finalize data preprocessing and cleaning | 2024-11-19 |
| **Week 2**        | Labeling Drug Interaction                | 2024-11-26 |
| **Week 3**        | Drug likability and interaction analysis | 2024-11-26 |
| **Week 4**        | Predictive modeling                      | 2024-12-03 |
| **Week 5**        | Finalize the results and build a website | 2024-12-10 |

## Organization within the Team
1. **Drug-Drug Interaction**: Ekaterina Borisova, Alexandra Krasnova
2. **Drug Likeness**: Ekaterina Borisova, Alexandra Krasnova
3. **Institution Analysis**: Lucas Nicolas, Svierre Djuve, Amine Zghal
