# Drug Interaction Analysis with respect to its origin and likeness Project

## Abstract
Does drug's origin matter in drug-drug interaction? The development of pharmaceuticals is a collaborative effort that often reflects the expertise and innovation of the institutions behind them. This project aims to investigate whether drugs developed within the same institutions exhibit superior interaction profiles compared to those developed independently. One may hypothesise that shared research environments, methodologies, and collaborative networks contribute to a higher likelihood of favorable interactions among drugs originating from the same institution. 

By analyzing DDI data from BindingDB and enriching it with institutional information, we will assess patterns and correlations in drug compatibility across origins. Beyond providing insights into institutional influence, our findings could guide future collaborative efforts in drug design, encouraging synergistic development and improving therapeutic outcomes. This research has the potential to reshape our understanding of how institutional dynamics influence the broader pharmaceutical landscape.

## Research Questions
1. **Institutional Influence on DDIs**: Do drugs developed in the same universities interact positively?
2. **Compatibility of Popular Drugs**: Do highly "likable" drugs, based on QED, exhibit good interactions?
3. **Predictive DDI Modeling**: Can we predict interaction outcomes for novel drug pairs, and propose new candidates for future research?

## Proposed Additional Datasets
- **DrugBank**: To cross-reference known DDIs and drug classification. Use the data on drug-interactions to classify based on description to either `synergetic` or `antagonistic` interaction. Use the `DrugBank ID of Ligand` to map the two databases.

## Methods
**Analysing Institutions** 
1.
**Analysing Likeness of Drug**
1.
**Drug-Drug Interaction**
1. Extract data on DDI from DrugBank and based on description cluster the interaction.
   1) First embed the description possibly with BioBert
   2) Cluster the data using K-Means 
2. Visualize the drug interaction based on institution location.
3. Try predicting the interaction by embeding smiles and training on DrugBank labeled data.
4. Use the developped model to predict the DDI of likable drugs in specific geographical locations


## Proposed Timeline
| Milestone         | Task Description                         | Deadline   |
|--------------------|-----------------------------------------|------------|
| **Week 1**        | Finalize data preprocessing and cleaning | 2024-11-19 |
| **Week 2**        | Labeling Drug Interaction                | 2024-11-26 |
| **Week 3**        | Drug likability and interaction analysis | 2024-11-26 |
| **Week 4**        | Predictive modeling                      | 2024-12-03 |
| **Week 5**        | Finalize the results and build a website | 2024-12-10 |

## Organization within the Team
1. **Drug-Drug Interaction**: Ekaterina Borisova, Alexandra Krasnova
2. **Drug Likeness**: Ekaterina Borisova, Alexandra Krasnova
3. **Institution Analysis**: 
4. 
