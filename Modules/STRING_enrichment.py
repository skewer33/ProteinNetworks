from   datetime import datetime
import functools
from   math import log
import os
import pandas as pd
import stringdb
import subprocess
from   tabulate import tabulate


def display_df(df):
    """
    function for displaying DataFrames (df). If IPNB is used, df will display with common IPNB function 'Display', else:
    it will display by print() function
    :param df: DataFrame
    :return:
    """
    try: display(df)
    except: print(df)

def Check_Value(val:[str, float, int], valid_values:set, valname:str, message='Wrong value123'):
    """
    function for check correctness of input value
    :param val: input value
    :param valid_values: set of valid values
    :param valname: group name of valid_values set. Or name of val variable
    :param message: Error message
    :return:
    """
    if val not in valid_values:
        if message == 'Wrong value123':
            message = f'Wrong value of "{valname}" variable! Choose one of {valid_values}'
        raise Exception(message)

def print_downline(line_length:int=40):
    """
    function prints line for titler decorator
    :param line_length:
    :return:
    """
    line = '_'*line_length
    print(f'{line}\n\n')

def print_upline(title:str, line_length:int=40):
    """
    function prints line and adds title for titler decorator
    :param line_length:
    :return:
    """
    line = '_'*line_length
    print(f'\t{title}\n{line}')


def titler(title: str, line_length=40):
    """
    Decorator added Title and edges of Paragraph
    :param title: title
    :param line_length: length of line (number of '_' symbols)
    :return:
    """
    def titler_decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            print_upline(title, line_length)
            original_result = func(*args, **kwargs)
            print_downline(line_length)

            return original_result
        return wrapper
    return titler_decorator


class EnrichmentAnalysis:
    types = {'UniProtID': 'queryItem', 'Gene': 'preferredName'}

    def __init__(self, data, enrichment = None, protein_id_type='UniProtID'):
        """
        EnrichmentAnalysis class conctructor.
        :param data: Dataframe containing the protein ID for analysis. It must contain either a "Gene" or "UniProtID" column'
        :param enrichment: Dataframe containing the results of previous enrichment analysis
        :param protein_id_type: type of protein ID. Valid Types
        """
        #check correctness of inputs
        self.protein_id_type = protein_id_type
        self._check_proteins_column(data)

        self.orig_data = data
        self.proteins = self.orig_data[self.protein_id_type]
        self.enrichment = enrichment # enrichment from previous analysis

    def _check_proteins_column(self, data):
        """
        the function checks the presence of the 'Gene' and 'UniProtID' columns in data
        :param data: data contains proteins` IDs
        :return:
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
        check genes in dataset which didn`t find by STRING (nomapped genes) and
        found by STRING but wasn`t in dataset (overmapped genes)

        :return: sorted list of nomapped and overmapped genes
        """
        nomapped = sorted(
            set(self.proteins.unique()).difference(set(self.genes_mapped[self.types[self.protein_id_type]].unique())))
        overmapped = sorted(
            set(self.genes_mapped[self.types[self.protein_id_type]].unique()).difference(set(self.proteins.unique())))
        return nomapped, overmapped

    @titler('DISCARDING DUPLICATES')
    def drop_duplicated_genes(self, silent=False):
        """
        function for droppig dublicated genes

        subset: (list) Only consider certain columns for identifying duplicates, by default use all columns.
        return: df of dropped genes
        """
        subset = self.protein_id_type
        len_orig_set = len(self.orig_data)
        duplicates = self.orig_data[self.orig_data.duplicated(subset=subset)]
        self.orig_data.drop_duplicates(subset=subset, inplace=True)
        self.proteins = self.orig_data[self.protein_id_type]
        if not silent:
            print(f'{len(duplicates)} of {len_orig_set} genes was dropped from original set')
            if len(duplicates) < 1:
                return duplicates
            elif len(duplicates) < 20:
                print('Dropped rows from original set:')
                print(duplicates[self.protein_id_type])
            elif len(duplicates) < 80:
                print('Dropped genes from original set:\n', *list(duplicates.protein_id_type))
        return duplicates

    def get_category_terms(self, category:str, term_type:set='id')->set:
        """
        function returns set of all terms in chosen category
        :param category: Name of category
        :param term_type: 'id' or 'description'.
                id - returns terms IDs of category (for example, GO terms)
                description - returns Description of IDs of category
        :return: set of terms
        """
        d_term = {'id': 'term', 'description': 'description'} # dict associate term_type and colnames of enrichment table
        valid_category = self._get_valid_category()
        Check_Value(category, valid_category, 'category')
        Check_Value(term_type, {'description', 'id'}, 'term_type')
        return set(self.enrichment[d_term[term_type]][self.enrichment.category == category])

    def _get_script_path(self):
        return os.path.dirname(os.path.abspath(__file__))

    def get_enrichment(self):
        """
        function performs enrichment analysis. Results store in self.enrichment
        :return:
        """
        self.enrichment = stringdb.get_enrichment(self.genes_mapped.queryItem) #get enrichment
        self.enrichment['enrich_score'] = self.enrichment.fdr.apply(lambda x: round(-log(x, 2), 1)) #get enrichment score

    def get_genes_of_term(self, term:str)-> list:
        """
        function get genes from enrichment table by target term
        :param term: target GO term from column 'term' in enrichment table
        :return: list of genes associated with target term
        """
        try:
            return self.enrichment.inputGenes[self.enrichment.term == term].to_list()[0].rstrip().strip().split(',')
        except: print('Term not found')

    def get_genes_by_localization(self, compartments: list, set_operation: str, save=False):
        """
        function for getting proteins localized in target compartments. You also can do common set operations
        under compartments genes
        Example: get_genes_by_localization([Nucleus, Cytosol], 'union') return proteins localized in Nucleus or Cytosol

        compartments: list of compartments. Will be attention:
            1) Capitalization of letters matters. Get available compartment names by calling "get_components_list()".
            2) Order of compartments matter if you want to get sets difference.
            For example:
                get_genes_by_localization(['Nucleus', 'Cytosol'], 'difference') return just nucleus proteins,
                get_genes_by_localization(['Cytosol', 'Nucleus'], 'union') return cytosol and nucleus proteins.
                get_genes_by_localization(['all', 'Nucleus'], 'difference') return all proteins except nucleus proteins.
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
        print(f'{len(loc_genes)} genes was founded\n')

        if save: # save genes in txt format (1 gene on 1 string)
            filename = 'Genes_' + '_'.join(compartments)
            if len(filename) > 45:
                filename = 'Genes_' + compartments[0][:62] + '_and_' + str(len(compartments) - 1) + '_compartments'
            filename += '_' + set_operation + '.txt'
            with open(filename, 'w+') as f:
                for term in loc_genes:
                    f.write(term + '\n')
            print(f'File {filename} successfully saved in {os.path.abspath(os.getcwd())}\n')

        return list(loc_genes)

    @titler('MAPPING GENES IN STRING')
    def get_mapped(self, species=9606):

        self.genes_mapped = stringdb.get_string_ids(self.proteins, species=species)
        self.nomapped_genes, self.overmapped_genes = self._find_nomapped_genes()
        print(
            f'{len(self.genes_mapped.queryItem.unique())} of {len(set(self.proteins.unique()))} unique genes were mapped\n')
        if len(self.nomapped_genes) < 80:
            print('List of nomapped genes:\n', list(self.nomapped_genes))
        if len(self.overmapped_genes) < 80:
            print('List of overmapped genes:\n', list(self.overmapped_genes))

    def prioretizingGO(self, terms: [list, set], organism='Human', domain='BP'):
        """
        function for prioretizing GO-terms using R script with GOxploreR package (doi:10.1038/s41598-020-73326-3)
        See 'RScript Prioretizing_GO.R'
        work with R.4-3.x. Yoy need to add RScript in PATH

        :param terms: list of GO-terms
        :param organism: name of target organism
        :param domain: name of domain in GO-graph. Available inputs: 'BP' - Biological Process
			            											 'CC' - Cellular Component
            														 "MF" - Molecular Functions
        :return: list of Prioretized GO terms
        """
        valid_organisms = {"Homo Sapiens", "Human", "Rattus Norvegicus", "Rat", "Mus Musculus", "Mouse",
                           "Danio Rerio", "Zebrafish", "Caenorhabditis Elegans", "Worm", "Arabidopsis Thaliana",
                           "Cress", "Saccharomyces Cerevisiae", "Yeast", "Schizosaccharomyces Pombe",
                           "Fission Yeast", "Drosophila Melanogaster", "Fruit Fly", "Escherichia Coli", "E.Coli"}
        Check_Value(organism, valid_organisms, 'organism')
        Check_Value(domain, {'BP', 'MF', 'CC'}, 'domain')

        self.save_table(pd.DataFrame(terms, columns=['Term']), 'temp_enrich_terms.csv', saveformat='csv', index=False)

        # Request to CMD to execute RScript
        command = 'Rscript'
        path2script = f'{self._get_script_path()}\\Prioretizing_GO.R'
        path2file = f'{os.getcwd()}\\temp_enrich_terms.csv'

        # Variable number of args in a list
        args = [path2file, 'Human', 'BP']
        # Build subprocess command
        cmd = [command, path2script] + args
        # check_output will run the command and store to result
        x = subprocess.check_output(cmd, universal_newlines=True)
        print(x)

        prior_terms = pd.read_csv(f'{os.getcwd()}\\temp_enrich_terms_output.csv')
        return list(prior_terms.Term)

    def proteins_participation_in_the_category(self, df, category, term_type='id', term_sep='\n'):
        """
        function check terms that proteins participated and make statistics table
        :param df: target DataFrame
        :param category: Name of category
        :param term_type: 'id' or 'description'.
                id - returns terms IDs of category (for example, GO terms)
                description - returns Description of IDs of category
        :param term_sep:
        :return:
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
        function displays  all terms and number of associated genes in category
        :param category: Name of category. You can check available category by calling 'show_enrichment_categories' method
        :param show: "all" or integer number. Number of strings to display
        :param sort_by: ["genes", "term"] - sort by number of genes (by descending) or term names (by ascending)
        :param save: Need to save? Choose True. By default, save in .xlsx format
        :param savename: work with save=True, name of file
        :param saveformat: format of saving file: 'xlsx' or 'csv'
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
        else: display_df(df.head(show))

        if save:
            if savename == 'terms':
                savename = category + '_' + savename + '_' + datetime.now().strftime('%m-%d-%Y')
            self.save_table(df, savename, saveformat=saveformat, index=False)

    @titler('ENRICHMENT CATEGORIES')
    def show_enrichment_categories(self):
        """
        function shown available enrichment categories for current dataset
        :return: None
        """
        table = []
        for term in self.enrichment.category.unique():
            table.append([term, len(self.enrichment[self.enrichment.category == term])])
        print(tabulate(table, headers=['Category', 'Number of terms'], tablefmt='orgtbl'))

    def show_enrichest_terms_in_category(self, category: str, count: int = 10, sort_by='fdr',
                                         save: bool = False, savename='enrichment', saveformat='xlsx'):
        """
        function shows top-%count of most enriched terms in %category
        :param category: Name of category. You can check available category by calling 'show_enrichment_categories' method
        :param count: count of terms you need to show
        :param sort_by: you can sort target list by one of 'fdr', 'p_value', 'number_of_genes' parameters
        :param save: Need to save? Choose True. By default, save in .xlsx format
        :param savename: work with save=True, name of file
        :param saveformat: format of saving file: 'xlsx' or 'csv'
        """
        valid_category = self._get_valid_category()
        Check_Value(category, valid_category, 'category')
        Check_Value(sort_by, {'fdr', 'p_value', 'number_of_genes'}, 'sort_by')

        table = self.enrichment[self.enrichment.category == category].sort_values(by=sort_by)
        if save:
            if savename == 'enrichment':
                savename += '_' + category + '_' + datetime.now().strftime('%m-%d-%Y')
            self.save_table(table.head(count), savename, saveformat=saveformat, index=False)
        print(f'ENRICHEST TERMS IN CATEGORY "{category}"')
        display_df(table.head(count).drop(['number_of_genes_in_background', 'ncbiTaxonId', 'preferredNames', 'p_value'], axis=1))
        return table


    @staticmethod
    def create_subframe_by_names(df, column: str, names: [list, tuple, set], add: str = 'first'):
        """
        function finds rows in original dataset and returns sub-dataframe including input names in selected column

        :param df: target DataFrame
        :param column: the selected column in which names will be searched
        :param names: list of target names whose records need to be found in the table
        :param add: ['first', 'last', 'all'] parameter of adding found rows.
                    'first' - add only the first entry
                    'last' - add only the last entry
                    'all' - add all entries
        :return: sub-dataframe including input names in selected column
        """
        Check_Value(add, {'first', 'last', 'all'}, add)

        def add_all(table, rows):
            return pd.concat([table, rows])

        def add_first(table, rows):
            table.loc[len(table)] = rows.iloc[0]
            return table

        def add_last(table, rows):
            table.loc[len(table)] = rows.iloc[-1]
            return table

        adding_method = {'first': add_first,
                         'last': add_last,
                         'all': add_all}

        new_df = pd.DataFrame(columns=df.columns)
        not_found_names = []
        for name in names:
            rows = df[df[column] == name]
            if len(rows) > 0:
                new_df = adding_method[add](new_df, rows)
            else: not_found_names.append(name)
        print(f'{len(not_found_names)} names were not found in the dataframe:\n')
        print('[', end='')
        print(*not_found_names, sep=', ', end='')
        print(']')

        return new_df

    @staticmethod
    def save_table(table, name, saveformat='xlsx', index:bool = True):
        """
        function for saving DataFrame tables
        :param table: DataFrame
        :param name: name of file
        :param saveformat: format of saving file: 'xlsx' or 'csv'
        :param index: show indexes in saved table?
        :return:
        """
        Check_Value(saveformat, {'csv', 'xlsx'}, 'saveformat')
        try:
            if saveformat == 'xlsx':
                if name[-5:] != '.xlsx' and name[-4:] != '.xls':
                    name += '.xlsx'
                table.to_excel(name, index=index)
            elif saveformat == 'csv':
                if name[-4:] != '.csv':
                    name += '.csv'
                table.to_csv(name, index=index, header=True)
            print(f'File {name} successfully saved in {os.path.abspath(os.getcwd())}\n')
        except PermissionError:
            print('Permission Denied Error: Access is denied. Close file if it`s open and try again')
        except:
            print('Saving file isn`t complete. If you rewrite file, close it and try again')