# ProteinNetworks

The library provides convenient tools for rapid analysis of gene ontology data, enrichment, and protein-protein interactions. It is based on [`stringdb`](https://pypi.org/project/stringdb/). Some functions require [R](https://www.r-project.org/) to be installed (see [`EnrichmentAnalysis.prioretizingGO()`](#prioretizingGO)).

## Library Modules

- **Enrichment** – functions for enrichment analysis
- **Protein Networks** – functions for working with and analyzing protein-protein interaction networks
- **Interactors** – functions for searching interactors for target proteins

## Installation

```bash
pip install git+https://github.com/skewer33/ProteinNetworks.git
```

## Contents

- [Enrichment](#Enrichment)
- [Protein Networks](#ProteinNetworksAnalysis)
- [Mapping Tools](#MappingTools)
- [Interaction](#InteractionTools)
- [Helper Functions](#Wrappers)
- [Usage Examples](#usage-example)
- [Links](#links)

---

# <a name='Enrichment'></a> Enrichment

Module: [`ProteinNetworks.enrichment`](ProteinNetworks/enrichment.py)

### Class [`EnrichmentAnalysis`](#classEnrichmentAnalysis)

**Constructor parameters:**
- `data`: DataFrame with protein identifiers (column "Gene" or "UniProtID")
- `enrichment`: DataFrame with results of previous enrichment analysis (optional)
- `protein_id_type`: type of protein identifier (`'Gene'` or `'UniProtID'`)

#### Methods:

- **create_subframe_by_names(df, column, names, add='first')**  
  Returns a subtable by a list of values in the selected column.

- **drop_duplicated_genes(subset=None, silent=False)**  
  Removes duplicate genes.

- **get_category_terms(category, term_type='id')**  
  Returns a set of all terms in the selected category (`id` or `description`).

- **get_enrichment()**  
  Performs enrichment analysis, result is saved in `self.enrichment`.

- **get_genes_by_localization(compartments, set_operation, save=False)**  
  Gets proteins localized in specified compartments, with set operations.

- **get_genes_of_term(term)**  
  Returns a list of genes associated with the selected term.

- **get_mapped(species=9606)**  
  Gene mapping: search for STRINGid by protein id.

- **prioretizingGO(terms, organism='Human', domain='BP')**  
  GO term prioritization using an R script (GOxploreR).

- **proteins_participation_in_the_category(df, category, term_type='id', term_sep='\n')**  
  Statistics on protein participation in category terms.

- **save_table(table, name, saveformat='xlsx', index=True)**  
  Saves DataFrame to file.

- **show_category_terms(category, show=10, sort_by='genes', save=False, savename='terms', saveformat='xlsx')**  
  Shows all terms and the number of associated genes in the category.

- **show_enrichest_terms_in_category(category, count=10, sort_by='fdr', save=False, savename='enrichment', saveformat='xlsx')**  
  Shows the top-%count most enriched terms in the category.

- **show_enrichment_categories()**  
  Shows available enrichment categories for the current dataset.

---

# <a name='ProteinNetworks'></a> Protein Networks

Module: [`ProteinNetworks.networks`](ProteinNetworks/networks.py)

### Class [`NetworkAnalysis`](#classNetworkAnalysis)

**Constructor parameters:**
- `data`: DataFrame with proteins and their interactions
- `directed`: directed graph (bool)
- `weighted`: weighted graph (bool)

#### Methods:

- **create_graph(data, directed=False, weighted=False)**  
  Creates a protein interaction graph (networkx.Graph).

- **get_subgraph_by_genes(genes)**  
  Returns a subgraph by a list of genes.

- **get_connected_components()**  
  Returns connected components of the graph.

- **get_degree_centrality()**  
  Node degree centrality.

- **get_betweenness_centrality()**  
  Betweenness centrality.

- **get_closeness_centrality()**  
  Closeness centrality.

- **draw_network(layout='spring', node_color='skyblue', with_labels=True, figsize=(10, 10))**  
  Network visualization.

---

# <a name='MappingTools'></a> Mapping Tools

Module: [`ProteinNetworks.mapping`](ProteinNetworks/mapping.py)

#### Main functions:

- **get_mapping(df, from_id, to_id, species=9606)**  
  Protein identifier mapping (get STRINGid by Gene or UniProtID). 
---

# <a name='InteractionTools'></a> Interaction Tools

Module: [`ProteinNetworks.interactions`](ProteinNetworks/interactions.py)

#### Main functions:

- **get_interactors_from_biogrid(genes, species=9606)**  
  Get interactions from the BioGRID database.

- **get_interactors_from_stringdb(genes, species=9606)**  
  Get interactions from STRINGdb.

- **get_interactionsTable_from_biogrid(genes, species=9606)**  
  Interaction table from BioGRID.

- **get_interactionsTable_from_stringdb(genes, species=9606)**  
  Interaction table from STRINGdb.

- **merging_interactors_stringdb_and_biogrid(genes, species=9606)**  
  Merge interaction data from both databases.

- **get_interactors(genes, source='stringdb', species=9606)**  
  Universal interface for obtaining interactions.

---

# <a name='Wrappers'></a> Helper Functions

Module: [`ProteinNetworks.wrappers`](ProteinNetworks/wrappers.py)

- **Check_Value(value, valid_values, name)**  
  Value validation.

- **save_table(table, name, saveformat='xlsx', index=True)**  
  Save table.

- **create_subframe_by_names(df, column, names, add='first')**  
  Search for rows by a list of names.

---


## Usage Example

```python
import pandas as pd
from ProteinNetworks import get_enrichment

df = pd.read_csv('your_proteins.csv')
enrich_obj = get_enrichment(df, protein_id_type='Gene', species=9606)
enrich_obj.show_enrichment_categories()
```
See more examples in [Google Colab](https://drive.google.com/file/d/1JlcrtDNwOVLuKmwDy4apfIpt7Mheu4cF/view?usp=sharing)

---

## Links

- [stringdb documentation](https://pypi.org/project/stringdb/)
- [GOxploreR](https://cran.r-universe.dev/GOxploreR/doc/manual.html)
- [Colab Notebook (example)](https://drive.google.com/file/d/1JlcrtDNwOVLuKmwDy4apfIpt7Mheu4cF/view?usp=sharing)