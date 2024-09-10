# ProteinNetworks

The library contains convenient tools for rapid analysis of gene ontology, enrichment and protein-protein interaction data. Based on the [`stringdb`](https://pypi.org/project/stringdb/) library. Some features require you to install [R](https://www.r-project.org/) to work (see [`EnrichmentAnalysis.prioretizingGO()`](#prioretizingGO))

### The module will contain 4 sets of tools:
  * **Enrichment Analysis** 
  * **Protein networks Analysis**
  * **Group comparing tools**
  * **Visualization tools**

## Get Started

`pip install -i https://test.pypi.org/simple/ ProteinNetworks==0.1.3`

## Contents:

* [Enrichment Analysis](#EnrichmentAnalysis)

  * module: [`ProteinNetworks.STRING_enrichment`](#STRING_enrichment)
    * class:  [`EnrichmentAnalysis`](#classEnrichmentAnalysis)
    
      methods:
      * [`EnrichmentAnalysis.create_subframe_by_names()`](#create)
      * [`EnrichmentAnalysis.drop_duplicated_genes()`](#drop_duplicated_genes)
      * [`EnrichmentAnalysis.get_category_terms()`](#get_category_terms)
      * [`EnrichmentAnalysis.get_enrichment()`](#get_enrichment)
      * [`EnrichmentAnalysis.get_genes_by_localization()`](#get_genes_by_localization)
      * [`EnrichmentAnalysis.get_genes_of_term()`](#get_genes_of_term)
      * [`EnrichmentAnalysis.get_mapped()`](#get_mapped)
      * [`EnrichmentAnalysis.prioretizingGO()`](#prioretizingGO)
      * [`EnrichmentAnalysis.proteins_participation_in_the_category()`](#proteins_participation_in_the_category)
      * [`EnrichmentAnalysis.save_table()`](#save_table)
      * [`EnrichmentAnalysis.show_category_terms()`](#show_category_terms)
      * [`EnrichmentAnalysis.show_enrichest_terms_in_category()`](#show_enrichest_terms_in_category)
      * [`EnrichmentAnalysis.show_enrichment_categories()`](#show_enrichment_categories)


_________________________


# <a name='EnrichmentAnalysis'></a> Enrichment Analysis
Contains a set of functions based on the stringdb library for gene ontology analysis and enrichment analysis
Look examples in [Colab Notebook](https://drive.google.com/file/d/1JlcrtDNwOVLuKmwDy4apfIpt7Mheu4cF/view?usp=sharing)


## <a name='STRING_enrichment'></a> ProteinNetworks.STRING_enrichment module


### <a name="classEnrichmentAnalysis"></a> *class* ProteinNetworks.STRING_enrichment.EnrichmentAnalysis *(data, enrichment=None, protein_id_type='UniProtID')*

Bases: `object`

EnrichmentAnalysis class.
* **Parameters:**
  * **data:** Dataframe containing the protein ID for analysis. It must contain either a “Gene” or “UniProtID” column’
  * **enrichment:** Dataframe containing the results of previous enrichment analysis
  * **protein_id_type:** type of protein ID. Valid Types

#### <a name="create"></a>*static* create_subframe_by_names(df, column: str, names: [<class 'list'>, <class 'tuple'>, <class 'set'>], add: str = 'first')

function finds rows in original dataset and returns sub-dataframe including input names in selected column

* **Parameters:**
  * **df** – target DataFrame
  * **column** – the selected column in which names will be searched
  * **names** – list of target names whose records need to be found in the table
  * **add** – [‘first’, ‘last’, ‘all’] parameter of adding found rows.
    ‘first’ - add only the first entry
    ‘last’ - add only the last entry
    ‘all’ - add all entries
* **Returns:**
  sub-dataframe including input names in selected column

#### <a name="drop_duplicated_genes"></a> drop_duplicated_genes(silent=False)

function for droppig dublicated genes
* **Parameters:**
  * **subset:** (list) Only consider certain columns for identifying duplicates, by default use all columns.
return: df of dropped genes

#### <a name="get_category_terms"></a> get_category_terms(category: str, term_type: str = 'id')

function returns set of all terms in chosen category
* **Parameters:**
  * **category:** Name of category
  * **term_type:** ‘id’ or ‘description’.

    > id - returns terms IDs of category (for example, GO terms) 
    > 
    > description - returns Description of IDs of category
* **Returns:**
  set of terms

#### <a name="get_enrichment"></a> get_enrichment()

function performs enrichment analysis. Results store in self.enrichment
* **Returns:** None

#### <a name="get_genes_by_localization"></a> get_genes_by_localization(compartments: list, set_operation: str, save=False)

function for getting proteins localized in target compartments. You also can do common set operations
under compartments genes
> Example: *get_genes_by_localization([Nucleus, Cytosol], ‘union’)*  - return proteins localized in Nucleus or Cytosol

* **Parameters:**
  * **compartments:** list of compartments. **Will be attention**:
    1. Capitalization of letters matters. Get available compartment names by calling *get_components_list()*.

    2. Order of compartments matter if you want to get sets difference.
  * **set_operation:** operation between sets. This means that the operations will be applied sequentially to all
         sets from the compartments. *[**A**, **B**, **C**], 'intersection' **->** **A** and **B** and **C***

    >  For example:
    > 
      > *get_genes_by_localization([‘Nucleus’, ‘Cytosol’], ‘difference’)* -  return just nucleus proteins,
       *get_genes_by_localization([‘Cytosol’, ‘Nucleus’], ‘union’)* - return cytosol and nucleus proteins.
       *get_genes_by_localization([‘all’, ‘Nucleus’], ‘difference’)* - return all proteins except nucleus proteins.

#### <a name="get_genes_of_term"></a> get_genes_of_term(term: str)

function get genes from enrichment table by target term
* **Parameters:**
  * **term**: target GO term from column ‘term’ in enrichment table
* **Returns:** list of genes associated with target term

#### <a name="get_mapped"></a> get_mapped(species=9606)
function makes gene mapping, it finds STRINGids by protein ids. It`s important for future analysis
* **Parameters:**
  * **species:** ID of organism. For example, Human species=9606
* **Returns:** None

#### <a name="prioretizingGO"></a> prioretizingGO(terms: [<class 'list'>, <class 'set'>], organism='Human', domain='BP')

function for prioretizing GO-terms using R script with [GOxploreR](https://cran.r-universe.dev/GOxploreR/doc/manual.html) package ([doi:10.1038/s41598-020-73326-3](https://www.nature.com/articles/s41598-020-73326-3))
See ‘RScript Prioretizing_GO.R’
work with R.4-3.x. Yoy need to add RScript in PATH

If you use this function in google-collab, you will have to install R-packages at the first launch.
This may take a long time (up to 20 minutes)

* **Parameters:**
  * **terms** – list of GO-terms
  * **organism** – name of target organism
  * **domain** – name of domain in GO-graph. Available inputs: ‘BP’ - Biological Process
    ‘CC’ - Cellular Component
    “MF” - Molecular Functions
* **Returns:**
  list of Prioretized GO terms

#### <a name="proteins_participation_in_the_category"></a> proteins_participation_in_the_category(df, category, term_type='id', term_sep='\\n')

function check terms that proteins participated and make statistics table
* **Parameters:**
  * **df:** target DataFrame
  * **category:** Name of category
  * **term_type:** ‘id’ or ‘description’.

    > id - returns terms IDs of category (for example, GO terms) 
    > 
    > description - returns Description of IDs of category
  * **term_sep:** terms connected with each protein will save in one cell. Choose separator beetwen terms
* **Returns:** None

#### <a name="save_table"></a> *static* save_table(table, name, saveformat='xlsx', index: bool = True)

function for saving DataFrame tables
* **Parameters:**
  * **table**: DataFrame
  * **name**: name of file
  * **saveformat**: format of saving file: ‘xlsx’ or ‘csv’
  * **index**: show indexes in saved table?
* **Returns:** None

#### <a name="show_category_terms"></a> show_category_terms(category: str, show: [<class 'int'>, <class 'str'>] = 10, sort_by='genes', save: bool = False, savename='terms', saveformat='xlsx')

function displays  all terms and number of associated genes in category
* **Parameters:**
  * **category:** Name of category. You can check available category by calling ‘show_enrichment_categories’ method
  * **show:** “all” or integer number. Number of strings to display
  * **sort_by:** [“genes”, “term”] - sort by number of genes (by descending) or term names (by ascending)
  * **save:** Need to save? Choose True. By default, save in .xlsx format
  * **savename:** work with save=True, name of file
  * **saveformat:** format of saving file: ‘xlsx’ or ‘csv’
* **Returns:** None

#### <a name="show_enrichest_terms_in_category"></a> show_enrichest_terms_in_category(category: str, count: int = 10, sort_by='fdr', save: bool = False, savename='enrichment', saveformat='xlsx')

function shows top-%count of most enriched terms in %category
* **Parameters:**
  * **category:** Name of category. You can check available category by calling ‘show_enrichment_categories’ method
  * **count:** count of terms you need to show
  * **sort_by:** you can sort target list by one of ‘fdr’, ‘p_value’, ‘number_of_genes’ parameters
  * **save:** Need to save? Choose True. By default, save in .xlsx format
  * **savename:** work with save=True, name of file
  * **saveformat:** format of saving file: ‘xlsx’ or ‘csv’
* **Returns:** None

#### <a name="show_enrichment_categories"></a> show_enrichment_categories()

function shown available enrichment categories for current dataset
* **Returns:** None


