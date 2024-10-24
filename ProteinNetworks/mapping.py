from   .wrappers import Check_Value, display_df, titler, print_upline, print_downline

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


def get_mapping(ids, return_nomapped=False, protein_id_type='Gene', species=9606):
    """
    Mapping protein IDs to STRING IDs

    Parameters
    ----------
    ids : pd.DataFrame
        pd.DataFrame contains ids in *protein_id_type* column
    return_nomapped : bool
        return list of not mapped genes
    protein_id_type : str
        type of protein ID. Valid Types: 'UniProtID', 'Gene'
    species : int
        species id, for example, human species id = 9606

    Returns
    -------
    mapped_genes : list
        list of mapped genes (or list of not mapped genes if return_nomapped = True)
    """
    
    obj = STRING_mapping(ids, protein_id_type)
    obj.drop_duplicated_genes(silent=True)
    obj.get_mapped(species)
    if return_nomapped:
        return obj.genes_mapped, obj.nomapped_genes
    else:
        return obj.genes_mapped


class STRING_mapping():
    types = {'UniProtID': 'queryItem', 'Gene': 'preferredName'}

    def __init__(self, data, protein_id_type='Gene'):
        """
        Mapping class conctructor.
        
        Parameters
        ----------
        data : pd.DataFrame
            Dataframe containing the protein ID for analysis. It must contain either a "Gene" or "UniProtID" column
        protein_id_type : str
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
        
    def _find_nomapped_genes(self):
        """
        Function for finding genes in dataset which didn't find by STRING (nomapped genes) and
        found by STRING but wasn't in dataset (overmapped genes)

        Parameters
        ----------
        None

        Returns
        -------
        tuple
            tuple of two lists. First list contains nomapped genes and second list contains overmapped genes.
        """

        nomapped = sorted(
            set(self.proteins.unique()).difference(set(self.genes_mapped[self.types[self.protein_id_type]].unique())))
        overmapped = sorted(
            set(self.genes_mapped[self.types[self.protein_id_type]].unique()).difference(set(self.proteins.unique())))
        return nomapped, overmapped

    @titler('DISCARDING DUPLICATES')
    def drop_duplicated_genes(self, silent=False):
        """
        Function for dropping duplicated genes.

        Parameters
        ----------
        silent : bool
            If True, function will not print information about dropped genes

        Returns
        -------
        DataFrame
            DataFrame with dropped genes
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
    
    @titler('MAPPING GENES IN STRING')
    def get_mapped(self, species=9606):
        """
        Function for mapping genes in STRING database

        Parameters
        ----------
        species : int
            ID of organism. For example, Human species=9606

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
                
    @staticmethod
    def create_subframe_by_names(df, column: str, names: [list, tuple, set], add: str = 'first'):
        """
        Function finds rows in original dataset and returns sub-dataframe including input names in selected column

        Parameters
        ----------
        df : pd.DataFrame
            target DataFrame
        column : str
            the selected column in which names will be searched
        names : list of str
            list of target names whose records need to be found in the table
        add : str
            ['first', 'last', 'all'] parameter of adding found rows.
            'first' - add only the first entry
            'last' - add only the last entry
            'all' - add all entries

        Returns
        -------
        pd.DataFrame
            sub-dataframe including input names in selected column
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
