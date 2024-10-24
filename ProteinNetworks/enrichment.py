from   .R_requests import Check_R_packages, short_R_output
from   .wrappers import Check_Value, display_df, titler, print_upline, print_downline, save_table, create_subframe_by_names

from   datetime import datetime
import functools
from   math import log
import os
import pandas as pd
import stringdb
from   subprocess import Popen, PIPE
from   tabulate import tabulate


# path to directory contains all RScripts
RSCRIPTS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'RScripts')


def get_enrichment(proteins, protein_id_type='Gene', species=9606, silent=True):
    """
    Function for one-click presetted enrichment analysis. It creates EnrichmentAnalysis object with your data,
    drops duplicates, maps protein IDs to STRING IDs, makes enrichment analysis and shows
    enrichment categories.

    Parameters
    ----------
    proteins : pd.DataFrame
        DataFrame with protein IDs. It must contain either a "Gene" or "UniProtID" column
    protein_id_type : str
        Type of protein ID. Valid types: 'Gene', 'UniProtID'
    species : int
        Species ID. For example, human species ID = 9606
    silent : bool
        If True, then function will not print anything

    Returns
    -------
    EnrichmentAnalysis
        EnrichmentAnalysis object
    """
    obj = EnrichmentAnalysis(proteins, protein_id_type=protein_id_type)
    obj.drop_duplicated_genes(silent=silent) # drop duplicates from your protein id set
    obj.get_mapped(species=species) # find STRINGid for each protein id
    obj.get_enrichment() # make enrichment
    obj.show_enrichment_categories() # get enrichment category
    return obj


class EnrichmentAnalysis:
    types = {'UniProtID': 'queryItem', 'Gene': 'preferredName'}

    def __init__(self, data, enrichment = None, protein_id_type='Gene'):
        """
        EnrichmentAnalysis class conctructor.

        Parameters
        ----------
        data : pd.DataFrame
            Dataframe containing the protein ID for analysis. It must contain either a "Gene" or "UniProtID" column
        enrichment : pd.DataFrame, optional
            Dataframe containing the results of previous enrichment analysis
        protein_id_type : str, optional
            type of protein ID. Valid Types

        Returns
        -------
        None
        """
        #check correctness of inputs
        self.protein_id_type = protein_id_type
        self._check_proteins_column(data)

        self.orig_data = data
        self.proteins = self.orig_data[self.protein_id_type]
        self.enrichment = enrichment # enrichment from previous analysis

    def _check_proteins_column(self, data):
        """
        The function checks the presence of the 'Gene' and 'UniProtID' columns in data

        Parameters
        ----------
        data : pd.DataFrame
            data contains proteins` IDs

        Returns
        -------
        None
        """

        valid_cols = set(self.types.keys()).intersection(data.columns)
        if len(valid_cols) == 0:
            message = 'The protein data must contain either a "Gene" or "UniProtID" column'
            del self
            raise Exception(message)
        elif len(valid_cols) == 1:
            p_id = next(iter(valid_cols))
            if self.protein_id_type != p_id:
                self.protein_id_type = p_id
                print(f'You choose "protein_id_type" that wasn`t contained in your data. '
                      f'"protein_id_type" is changed to "{self.protein_id_type}"\n')
        elif len(valid_cols) == 2:
            pass

    def _get_valid_category(self)->set:
        """
        function return set of valid category names for current enrichment analysis
        :return:
        """
        return set(self.enrichment.category.unique())

    def _find_nomapped_genes(self):
        """
        check genes in dataset which didn't find by STRING (nomapped genes) and
        found by STRING but wasn`t in dataset (overmapped genes)
        
        Parameters
        ----------
        None

        Returns
        -------
        tuple
            tuple of two sorted lists. First list contains nomapped genes and second list contains overmapped genes.
        """
        nomapped = sorted(
            set(self.proteins.unique()).difference(set(self.genes_mapped[self.types[self.protein_id_type]].unique())))
        overmapped = sorted(
            set(self.genes_mapped[self.types[self.protein_id_type]].unique()).difference(set(self.proteins.unique())))
        return nomapped, overmapped
    
    @titler('DISCARDING DUPLICATES')
    def drop_duplicated_genes(self, silent=False):
        """
        Function for droppig dublicated genes
        
        Parameters
        ----------
        subset : list, optional
            Only consider certain columns for identifying duplicates, by default use all columns.

        Returns
        -------
        pd.DataFrame
            df of dropped genes
        """
        
        subset = self.protein_id_type
        len_orig_set = len(self.orig_data)
        duplicates = self.orig_data[self.orig_data.duplicated(subset=subset)]
        self.orig_data.drop_duplicates(subset=subset, inplace=True)
        self.proteins = self.orig_data[self.protein_id_type]
        if not silent:
            print(f'{len(duplicates)} of {len_orig_set} genes were dropped from original set')
            if len(duplicates) < 1:
                return duplicates
            elif len(duplicates) < 20:
                print('Dropped rows from original set:')
                print(duplicates[self.protein_id_type])
            elif len(duplicates) < 80:
                print('Dropped genes from original set:\n', *list(duplicates.protein_id_type))
        return duplicates

    def get_category_terms(self, category:str, term_type:str='id')->set:
        """
        Parameters
        ----------
        category : str
            Name of category
        term_type : str
            'id' or 'description'.
            id - returns terms IDs of category (for example, GO terms)
            description - returns Description of IDs of category

        Returns
        -------
        set
            set of terms
        """
        d_term = {'id': 'term', 'description': 'description'} # dict associate term_type and colnames of enrichment table
        valid_category = self._get_valid_category()
        Check_Value(category, valid_category, 'category')
        Check_Value(term_type, {'description', 'id'}, 'term_type')
        return set(self.enrichment[d_term[term_type]][self.enrichment.category == category])

    def get_enrichment(self):
        """
        Function performs enrichment analysis. Results store in self.enrichment
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        self.enrichment = stringdb.get_enrichment(self.genes_mapped.queryItem) #get enrichment
        self.enrichment['enrich_score'] = self.enrichment.fdr.apply(lambda x: round(-log(x, 2), 1)) #get enrichment score

    def get_genes_of_term(self, term:str)-> list:
        """
        Parameters
        ----------
        term : str
            target GO term from column 'term' in enrichment table

        Returns
        -------
        list
            list of genes associated with target term
        """
        try:
            return self.enrichment.inputGenes[self.enrichment.term == term].to_list()[0].rstrip().strip().split(',')
        except: print('Term not found')

    def get_genes_by_localization(self, compartments: list, set_operation: str, save=False):
        """
        Function for getting proteins localized in target compartments. You also can do common set operations
        under compartments genes
        Example: get_genes_by_localization([Nucleus, Cytosol], 'union') return proteins localized in Nucleus or Cytosol

        Parameters
        ----------
        compartments : list
            list of compartments. Will be attention:
                1) Capitalization of letters matters. Get available compartment names by calling "get_components_list()".
                2) Order of compartments matter if you want to get sets difference.
        set_operation : str
            operation between sets. This means that the operations will be applied sequentially to all
            sets from the compartments.
            For example:
                get_genes_by_localization(['Nucleus', 'Cytosol'], 'difference') return just nucleus proteins,
                get_genes_by_localization(['Cytosol', 'Nucleus'], 'union') return cytosol and nucleus proteins.
                get_genes_by_localization(['all', 'Nucleus'], 'difference') return all proteins except nucleus proteins.

        Returns
        -------
        set of proteins localized in target compartments
        """
        # check set_operation and compartments error
        Check_Value(set_operation, {'union', 'intersection', 'difference', 'symmetric_difference'}, 'set_operation')
        components_list = self.get_category_terms('Component', term_type='description').union({'all'})
        for c in compartments:
            Check_Value(c, components_list, 'Compartments',
                        message='There is no such compartment. To display a list of available compartments, '
                                'call <<show_category_terms("Components")>>. '
                                'If you want to get all genes, use tag "all" in compartments list')

        # define common set operations
        def union(a: set, b: set) -> set:
            return a.union(b)

        def intersection(a: set, b: set) -> set:
            return a.intersection(b)

        def difference(a: set, b: set) -> set:
            return a.difference(b)

        def symmetric_difference(a: set, b: set) -> set:
            return a.symmetric_difference(b)

        operations = {'union': union, 'intersection': intersection, 'difference': difference,
                      'symmetric_difference': symmetric_difference}

        # choose all strings in enrichment data connected with 'Component'
        component_data = self.enrichment[self.enrichment.category == 'Component']

        # create location genes set and apply set_operation for each compartment gene set
        if compartments[0] == 'all':
            loc_genes = set(self.proteins)
        else:
            loc_genes = set(component_data.inputGenes[component_data.description == compartments[0]]
                            .to_list()[0].rstrip().strip().split(','))
        for i in range(1, len(compartments)):
            if compartments[i] == 'all':
                compartment_genes = set(self.proteins)
            else:
                compartment_genes = set(component_data.inputGenes[component_data.description == compartments[i]]
                                        .to_list()[0].rstrip().strip().split(','))
            loc_genes = operations[set_operation](loc_genes, compartment_genes)
        print(f'{len(loc_genes)} genes were founded\n')

        if save: # save genes in txt format (1 gene on 1 string)
            filename = 'Genes_' + '_'.join(compartments)
            if len(filename) > 45:
                filename = 'Genes_' + compartments[0][:62] + '_and_' + str(len(compartments) - 1) + '_compartments'
            filename += '_' + set_operation + '.txt'
            with open(filename, 'w+') as f:
                for term in loc_genes:
                    f.write(term + '\n')
            print(f'File {filename} successfully saved in {os.path.abspath(os.getcwd())}\n')

        #return list(loc_genes)
        return pd.DataFrame({self.protein_id_type: list(loc_genes)})

    @titler('MAPPING GENES IN STRING')
    def get_mapped(self, species=9606):
        """
        Parameters
        ----------
        species : int, optional
            ID of organism. For example, Human species=9606. Default is 9606.

        Returns
        -------
        None
        """

        self.genes_mapped = stringdb.get_string_ids(self.proteins, species=species)
        self.nomapped_genes, self.overmapped_genes = self._find_nomapped_genes()
        print(
            f'{len(self.genes_mapped.queryItem.unique())} of {len(set(self.proteins.unique()))} unique genes were mapped\n')
        if len(self.nomapped_genes) < 80:
            print('List of nomapped genes:\n', list(self.nomapped_genes))
        if len(self.overmapped_genes) < 80:
            print('List of overmapped genes:\n', list(self.overmapped_genes))

    def get_terms_prioretizing(self, category, sorting_results=True, ascending=False):
        
        """
        Function for prioretizing GO-terms from choosen category
        
        Parameters
        ----------
        obj : EnrichmentAnalysis
            EnrichmentAnalysis object
        category : str
            name of category

        Returns
        -------
        **pd.DataFrame** enrichment table of prioretized Terms, stored in *self.prior_enrichment[category]*
        """
        #check validness of category
        valid_category = self._get_valid_category()
        Check_Value(category, valid_category, 'category')
        
        GO_terms = self.get_category_terms(category)
        prior_GO_terms = self.prioretizingGO(GO_terms)
        try:
            self.prior_enrichment[category] = create_subframe_by_names(self.enrichment, column='term', names=prior_GO_terms)
        except:
            self.prior_enrichment = {}
            self.prior_enrichment[category] = create_subframe_by_names(self.enrichment, column='term', names=prior_GO_terms)
        if sorting_results:
            self.prior_enrichment[category].sort_values(by=['enrich_score'], ascending=ascending, inplace=True)
        return self.prior_enrichment[category]
        

    def prioretizingGO(self, terms: [list, set], organism='Human', domain='BP'):
        """
        Function for prioretizing GO-terms using R script with GOxploreR package (doi:10.1038/s41598-020-73326-3)
        See 'RScript Prioretizing_GO.R'
        work with R.4-3.x. You need to add RScript in PATH

        If you use this function in google-collab, you will have to install R-packages at the first launch.
        This may take a long time (up to 20 minutes)

        Parameters
        ----------
        terms : list or set
            list of GO-terms
        organism : str
            name of target organism
        domain : str
            name of domain in GO-graph. Available inputs: 'BP' - Biological Process
                                                     'CC' - Cellular Component
                                                     "MF" - Molecular Functions

        Returns
        -------
        list of Prioretized GO terms
        """
        valid_organisms = {"Homo Sapiens", "Human", "Rattus Norvegicus", "Rat", "Mus Musculus", "Mouse",
                           "Danio Rerio", "Zebrafish", "Caenorhabditis Elegans", "Worm", "Arabidopsis Thaliana",
                           "Cress", "Saccharomyces Cerevisiae", "Yeast", "Schizosaccharomyces Pombe",
                           "Fission Yeast", "Drosophila Melanogaster", "Fruit Fly", "Escherichia Coli", "E.Coli"}
        Check_Value(organism, valid_organisms, 'organism')
        Check_Value(domain, {'BP', 'MF', 'CC'}, 'domain')

        installing = Check_R_packages(CRAN_packages=["GOxploreR", "data.table", "BiocManager", "utils", "ggplot2"],
                         BiocManager_packages=["GO.db", "annotate", "biomaRt"])

        save_table(pd.DataFrame(terms, columns=['Term']), 'input_priority_terms.csv', saveformat='csv', index=False)

        # Request to CMD to execute RScript
        command = 'Rscript'
        path2script = os.path.join(RSCRIPTS_PATH, 'Prioretizing_GO.R')
        path2file = os.path.abspath('input_priority_terms.csv')

        # Variable number of args in a list
        args = [path2file, 'Human', 'BP']
        # Build subprocess command
        cmd = [command, path2script] + args
        # check_output will run the command and store to result
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, error = p.communicate()

        # PRINT R CONSOLE OUTPUT (ERROR OR NOT)
        if p.returncode == 0:
            if installing:
                print('All R-packages were installed successfully')
            s_output = short_R_output(output.decode("utf8")) # if all is OK, then makes short output
            print(f'R OUTPUT:\n {s_output}')
        else:
            print(f'R ERROR:\n {error.decode("utf8")}')

        prior_terms = pd.read_csv(os.path.abspath('output_priority_terms.csv'))
        return list(prior_terms.Term)

    def proteins_participation_in_the_category(self, df, category, term_type='id', term_sep='\n'):
        """
        Function check terms that proteins participated and make statistics table

        Parameters
        ----------
        df : pd.DataFrame
            target DataFrame
        category : str
            Name of category
        term_type : str
            'id' or 'description'.
            id - returns terms IDs of category (for example, GO terms)
            description - returns Description of IDs of category
        term_sep : str
            terms connected with each protein will save in one cell. Choose separator beetwen terms

        Returns
        -------
        pd.DataFrame
            table with protein`s participation in category
        """
        d_term = {'id': 'term', 'description': 'description'} # dict associate term_type and colnames of enrichment table
        valid_category = self._get_valid_category()
        Check_Value(category, valid_category, 'category')
        Check_Value(term_type, {'description', 'id'}, 'term_type')

        prot_participation = pd.DataFrame(columns=[self.protein_id_type, 'number_of_terms', 'terms'])
        sub_df_category = df[df.category == category]
        for prot in self.proteins:
            prot_participation.loc[len(prot_participation)] = {
                self.protein_id_type: prot,
                'number_of_terms': len(sub_df_category[sub_df_category.inputGenes.str.contains(prot)]),
                'terms': sub_df_category[d_term[term_type]][sub_df_category.inputGenes.str.contains(prot)].apply(
                    lambda x: str(x) + term_sep).sum()}

        prot_participation.sort_values('number_of_terms', ascending=False, inplace=True)
        return prot_participation

    @titler('CATEGORY TERMS')
    def show_category_terms(self, category:str, show:[int, str]=10, sort_by='genes',
                            save:bool = False, savename='terms', saveformat='xlsx')->None:
        """
        Function displays all terms and number of associated genes in category

        Parameters
        ----------
        category : str
            Name of category. You can check available category by calling 'show_enrichment_categories' method
        show : int or str
            "all" or integer number. Number of strings to display
        sort_by : str
            ["genes", "term"] - sort by number of genes (by descending) or term names (by ascending)
        save : bool
            Need to save? Choose True. By default, save in .xlsx format
        savename : str
            work with save=True, name of file
        saveformat : str
            format of saving file: 'xlsx' or 'csv'
        """

        if type(show) != int and show != 'all':
            raise Exception('Error of "show" variable. Choose "all" or integer number')
        valid_category = self._get_valid_category()
        Check_Value(category, valid_category, 'category')
        Check_Value(sort_by, {'genes', 'term'}, 'sort_by')

        table = []
        category_data = self.enrichment[self.enrichment.category == category]
        terms = self.get_category_terms(category, term_type='description')
        
        for term in terms:
            string = category_data[category_data.description == term]
            table.append([term, list(string.number_of_genes)[0]])
        df = pd.DataFrame(table, columns=['Term', '# Genes'])
        if sort_by == 'genes':
            df.sort_values('# Genes', ascending = False, inplace=True)
        elif sort_by == 'term':
            df.sort_values('Term', ascending = True, inplace=True)
        if show == 'all':
            display_df(df)
        else:
            display_df(df.head(show))

        if save:
            if savename == 'terms':
                savename = category + '_' + savename + '_' + datetime.now().strftime('%m-%d-%Y')
            self.save_table(df, savename, saveformat=saveformat, index=False)

    @titler('ENRICHMENT CATEGORIES')
    def show_enrichment_categories(self):
        """
        function shown available enrichment categories for current dataset
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        table = []
        for term in self.enrichment.category.unique():
            table.append([term, len(self.enrichment[self.enrichment.category == term])])
        print(tabulate(table, headers=['Category', 'Number of terms'], tablefmt='orgtbl'))

    def show_enrichest_terms_in_category(self, category: str, count: [int, str]='all', sort_by='fdr',
                                         save: bool = False, savename='enrichment', saveformat='xlsx'):
        """
        Function shows top-%count of most enriched terms in %category

        Parameters
        ----------
        category : str
            Name of category. You can check available category by calling 'show_enrichment_categories' method
        count : int, str
            count of terms you need to show. Choose "all" or integer number
        sort_by : str
            you can sort target list by one of 'fdr', 'p_value', 'number_of_genes' parameters
        save : bool
            Need to save? Choose True. By default, save in .xlsx format
        savename : str
            work with save=True, name of file
        saveformat : str
            format of saving file: 'xlsx' or 'csv'
        Returns
        -------
        None
        """
        
        if type(count) != int and count != 'all':
            raise Exception('Error of "count" variable. Choose "all" or integer number')
        valid_category = self._get_valid_category()
        Check_Value(category, valid_category, 'category')
        Check_Value(sort_by, {'fdr', 'p_value', 'number_of_genes'}, 'sort_by')

        table = self.enrichment[self.enrichment.category == category].sort_values(by=sort_by)
        if count == 'all':
            count = len(table)
            
        if save:
            if savename == 'enrichment':
                savename += '_' + category + '_' + datetime.now().strftime('%m-%d-%Y')
            save_table(table.head(count), savename, saveformat=saveformat, index=False)
        print(f'ENRICHEST TERMS IN CATEGORY "{category}"')
        display_df(table.head(count).drop(['number_of_genes_in_background', 'ncbiTaxonId', 'preferredNames', 'p_value'], axis=1))
        return table
