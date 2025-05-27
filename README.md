# Code and data for *Homologous nodes in annotated complex networks*

## About:
This is a repository for all code, analysis and data accompanying the Moon & Ahnert paper available [here](https://arxiv.org/). 

## Installation:
- Run:
```
conda env create -f vectorisation.yml
conda activate vectorisation_env
```
- Download the data folder from [here] and add to this directory. 
- Then run Jupyter notebooks to do analysis and generate figures into the ```figures``` folder. 

Note:
- Requires ```Graphviz```  to be installed for ```dot``` layout in ```networkx``` ([see here](https://graphviz.org/download/)). 


## Contents:
### ```data```
####  ```/grn```: Graph regulatory network and associated metadata
- ```/gene_association.tair``` is the gene ontology (GO) annotations associated with each gene downloaded from [The Arabidopsis Information Resource (TAIR)](https://www.arabidopsis.org/download/list?dir=GO_and_PO_Annotations%2FGene_Ontology_Annotations), Berardini et al., *Plant Physiology* (2004). Version generated: 2025-01-01. 
- ```/AtRegNet.csv``` is the gene regulatory network of the *Arabidopsis Thalania*, available from [agris](https://agris-knowledgebase.org/downloads.html), Palaniswamy et al., *Plant Physiology* (2006). Version dated: 2019-03-11. 
- ```/tf_family.csv``` Transcription factor gene list for family annotation metadata from AtTFDB in [agris](https://agris-knowledgebase.org/AtTFDB), Palaniswamy et al., *Plant Physiology* (2006). 
 - ```/ATH_GO_GOSLIM.txt``` is the gene ontology annotation explanatory text file, downloaded from [The Arabidopsis Information Resource (TAIR)](https://www.arabidopsis.org/download/list?dir=GO_and_PO_Annotations%2FGene_Ontology_Annotations), Berardini et al., *Plant Physiology* (2004). Version generated: 2025-03-01. 

#### ```/foodweb``` : ecological food web network and associated metadata
- ```/florida.bare_upto122_fullnames``` is a text file of the prey to predator edge list, from Heymans et al., *Ecological Modelling*, (2002).  
- ```florida.bare_upto122_fullnames.newcat.terms``` is a text file of the organism and organism type annotations, from Heymans et al., *Ecological Modelling*, (2002). 

#### ```/recipe``` : cuisine, ingredient and recipe mulitpartite network
- ```/cuisine_recipe``` is a text file of cuisines types and the recipe IDs that belong to the cuisine, from allrecipes recipe database, via Ahn et al., *Scientific Reports*, (2011).
- ```/recipe_ingredients_bipartite``` is a text file containing the edge list of recipe IDs to cuisine labels, from allrecipes recipe database, via Ahn et al., *Scientific Reports*, (2011).
- ```/foodtype_categorised.csv``` is a csv of manual categorisation of ingredient type. 
- ```/cuisine_geo_labels.csv``` is a csv of manual cuisine categorisation into geographical regions of origin. 

### Enriched clusters in the GRN
- The full ```.csv``` file of the enriched clusters, their TF family labels, and the GO slim descriptors can be found in ```grn_GO_enrichment.csv```, (also in ```data/grn/grn_GO_enrichment.csv```).

### ```src```
- ```ipynb``` notebooks are here. 
    - ```00_visualisation.ipynb``` visualises the multipartite recipe network. **Figure 2C, 3B** 
    - ```01_vectorisation.ipynb``` performs the connectivity aggregation and vectorisation for the three networks. Saves these vectors (optional) into a ```processed``` subfolder in ```data/*/```. 
    - ```02_clustering.ipynb``` clusters and visualises these clusters into a dendrogram and associated vector heatmap. **Figures 1(A,B), 2(A,B), 3(A), 4, S1, S2, S3.**
        - For the gene regulatory network, it performs the enrichment analysis for the clusters of transcription factors. 
    - ```03_null_models_nb.ipynb``` includes null models of the specialization-diversity entropic vector euclidean pairwise distances for cell types, serial homologues and left-right pairs. **Figures 1D, 2D, 3B, S4**

#### ```/module/```
- Scripts and functions to run the analysis:
    - ```fvec.py```: a fast vectorisation routine. 
    - ```graphpeeler.py```: helpers for probabilistic topological layer sorting. 
    - ```enrichement_utils.py```: helpers for performing GO term enrichment analysis. 
    - ```data_prep.py```: preparation of data from raw sources. 
    - ```null_helpers.py```: null model wrappers used in ```src/03_null_models.ipynb```.

### ```figures```
- (Sub)figures are organised into their respective network folders. 
- Fully annotated clustering figures are found in:
    - ```/fw/f2a_dendrogram_full.pdf``` for the food web network. 
    - ```/rn/f3a_dendrogram_full.pdf``` for the recipe network. 
    - ```/grn/s1_TF_GRN_GO_clustering_BH_all_GO_annotated.pdf``` for the gene regulatory network. 

**References**

Ahn, Y.-Y., Ahnert, S. E., Bagrow,  J. P., Barab´asi, A.-L., (2011) **Flavor network and the principles of food pairing**, *Scientific Reports 1, 196*

Berardini, T.Z., Mundodi, S., Reiser, R., Huala, E., Garcia-Hernandez, M., Zhang, P., Mueller, L.M., Yoon, J., Doyle, A., Lander, G., Moseyko, N., Yoo, D., Xu, I., Zoeckler, B., Montoya, M., Miller, N., Weems, D., Rhee, S.Y. (2004) **Functional annotation of the Arabidopsis genome using controlled vocabularies**. *Plant Physiology 135(2):745–755*.

Heymans, J.J., Ulanowicz, R.E., Bondavalli, C., (2003) **Network analysis of the South Florida Everglades graminoid marshes and comparison with nearby cypress ecosystems**, *Ecological Modelling 149(1–2):5-23*.


Palaniswamy, S.K., James, S., Sun, H., Lamb, R.S., Davuluri, R.V., Grotewold, E. (2006) **AGRIS and AtRegNet: A platform to link cis-regulatory elements and transcription factors into regulatory networks.** *Plant Physiology, 140(3):818-829*.