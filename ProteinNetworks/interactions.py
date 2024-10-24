from   .networks import NetworkAnalysis
from   .wrappers import Check_Value, Check_type, Check_kwargs

from biogridpy.biogrid_client import BioGRID
import os
import pandas as pd
import stringdb
import numpy as np
from io import StringIO


BIOGRID_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BioGRID')

conf_path = os.path.join(BIOGRID_PATH, 'biogridpyrc')
INTERACTIONS_TYPE = {'STRINGdb':{
                        'physical': {'included_metrics': ['dscore', 'escore', 'tscore']}, 
                        'genetic': {'included_metrics': ['nscore', 'fscore', 'pscore', 'ascore']},
                        'all':{'included_metrics': ['nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore']}},
                    'BioGRID':{
                        'physical': {'EXPERIMENTAL_SYSTEM_TYPE': 'physical'}, 
                        'genetic': {'EXPERIMENTAL_SYSTEM_TYPE': 'genetic'},
                        'all': {'EXPERIMENTAL_SYSTEM_TYPE': 'any'}} 
                    }

def get_interactionsTable_from_biogrid(geneList:list,
                                       taxId=9606,
                                       interactions_type = None,
                                       EXPERIMENTAL_SYSTEM_TYPE = 'any',
                                       required_evidence = 1):
    """
    Function returns interaction table from BioGRID for given list of genes and organism id (taxId).

    Parameters
    ----------
    geneList : list
        list of genes
    taxId : int
        organism id. Human = 9606
    EXPERIMENTAL_SYSTEM_TYPE : str
            type of experimental system. Valid types are 'genetic', 'physical' or 'any'.
            Default is 'any'
    required_evidence : int
        required number of evidence for interaction
        Default is 1. It means that BioGRID has 1 mention about interaction between two proteins
    interactions_type : str
        None by default. Which type of interactions you want to search. It`s ready to use presets. 'all', 'physical' or 'genetic'.
        Activated only if interactions_type is not None. When activated, EXPERIMENTAL_SYSTEM_TYPE is ignored

    Returns
    -------
    pd.DataFrame
        interaction table
    """
    if interactions_type is not None:
        Check_Value(interactions_type, {'genetic', 'physical', 'all'}, 'interactions_type')
        EXPERIMENTAL_SYSTEM_TYPE = INTERACTIONS_TYPE['BioGRID'][interactions_type]['EXPERIMENTAL_SYSTEM_TYPE']
    BGI = BiogridInteractions(geneList, taxId)
    BGI.get_interaction_partners(EXPERIMENTAL_SYSTEM_TYPE, required_evidence = required_evidence)
    return BGI.interactionTable

def get_interactors_from_biogrid(geneList:list, 
                                 taxId=9606, 
                                 EXPERIMENTAL_SYSTEM_TYPE = 'any',
                                 required_evidence = 1,
                                 interactions_type = None)->list:
    
    """
    Function returns list of interactors from BioGRID for given list of genes and organism id (taxId).

    Parameters
    ----------
    geneList : list
        list of genes
    taxId : int
        organism id. Human = 9606
    EXPERIMENTAL_SYSTEM_TYPE : str
        type of experimental system. Valid types are 'genetic', 'physical' or 'any'.
        Default is 'any'
    required_evidence : int
        required number of evidence for interaction
        Default is 1. It means that BioGRID has 1 mention about interaction between two proteins
    interactions_type : str
        None by default. Which type of interactions you want to search. It`s ready to use presets. 'all', 'physical' or 'genetic'.
        Activated only if interactions_type is not None. When activated, EXPERIMENTAL_SYSTEM_TYPE isignored

    Returns
    -------
    pd.DataFrame
        1-column DataFrame of interactors (header = Gene)
    """
    Check_Value(EXPERIMENTAL_SYSTEM_TYPE, {'any', 'genetic', 'physical'}, 'EXPERIMENTAL_SYSTEM_TYPE')
    if interactions_type is not None:
        Check_Value(interactions_type, {'genetic', 'physical', 'all'}, 'interactions_type')
        EXPERIMENTAL_SYSTEM_TYPE = INTERACTIONS_TYPE['BioGRID'][interactions_type]['EXPERIMENTAL_SYSTEM_TYPE']
        
    BGI = BiogridInteractions(geneList, taxId)
    BGI.get_interaction_partners(EXPERIMENTAL_SYSTEM_TYPE, required_evidence = required_evidence)
    return BGI.interactorsList


def get_interactionsTable_from_stringdb(geneList:list, 
                                        required_score=400,
                                        included_metrics=['nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore'],
                                        limit=None,
                                        taxId=9606,
                                        interactions_type = None)->pd.DataFrame:
  
    """
    Function returns interaction table from STRINGdb for given list of genes and organism id (taxId) and
    required score of interaction.

    Parameters
    ----------
    geneList : list
        list of genes
    required_score : int
        required score of interaction. From 150 to 999
    limit : int or None
        limit of interactors for each gene from geneList
    taxId : int
        organism id. Human = 9606
    interactions_type : str
        None by default. Which type of interactions you want to search. It`s ready to use presets. 'all', 'physical' or 'genetic'.
        Activated only if interactions_type is not None. When activated, included metrics is ignored

    Returns
    -------
    pd.DataFrame
        interaction table
    """
    if interactions_type is not None:
        Check_Value(interactions_type, {'genetic', 'physical', 'all'}, 'interactions_type')
        included_metrics = INTERACTIONS_TYPE['STRINGdb'][interactions_type]['included_metrics']
    STRINGi = StringdbInteractions(geneList, taxId)
    STRINGi.get_interaction_partners(required_score=required_score, limit=limit, included_metrics=included_metrics)
    return STRINGi.interactionTable

def get_interactors_from_stringdb(geneList:list,
                                  required_score=400,
                                  included_metrics=['nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore'],
                                  limit=None,
                                  taxId=9606,
                                  interactions_type = None)->pd.DataFrame:
    """
    Function returns list of interactors from STRINGdb for given list of genes.

    Parameters
    ----------
    geneList : list
        list of genes
    required_score : int
        required score of interaction. From 150 to 999
    limit : int or None
        limit of interactors for each gene from geneList
    taxId : int
        organism id. Human = 9606
    interactions_type : str
        None by default. Which type of interactions you want to search. It`s ready to use presets. 'all', 'physical' or 'genetic'.
        Activated only if interactions_type is not None. When activated, included metrics is ignored

    Returns
    -------
    pd.DataFrame
        1-column DataFrame of interactors (header = Gene)
    """

    if interactions_type is not None:
        Check_Value(interactions_type, {'genetic', 'physical', 'all'}, 'interactions_type')
        included_metrics = INTERACTIONS_TYPE['STRINGdb'][interactions_type]['included_metrics']    
    
    STRINGi = StringdbInteractions(geneList, taxId)
    STRINGi.get_interaction_partners(required_score=required_score, limit=limit, included_metrics=included_metrics)
    return STRINGi.interactorsList

def merging_interactors_stringdb_and_biogrid(geneList1, geneList2, method= 'union', mergeCol='Gene'):
    
    """
    Function merges two lists of interactors (from STRINGdb and BioGRID) using union or intersection operation.

    Parameters
    ----------
    geneList1 : pd.DataFrame, list, tuple, pd.Series
        list of genes from STRINGdb
    geneList2 : pd.DataFrame, list, tuple, pd.Series
        list of genes from BioGRID
    method : str
        method of merging. 'union' or 'intersection'. 'union' by default
    mergeCol : str
        column name to merge. 'Gene' by default

    Returns
    -------
    pd.DataFrame
        merged list of interactors
    """
    Check_Value(method, {'union', 'intersection'}, 'method')
    Check_type(geneList1, (pd.DataFrame, pd.Series, list, tuple), 'geneList1')
    Check_type(geneList2, (pd.DataFrame, pd.Series, list, tuple), 'geneList2')
    
    if isinstance(geneList1, (list, pd.Series, tuple)):
        geneList1 = pd.DataFrame({'Gene': geneList1})  
        
    if isinstance(geneList2, (list, pd.Series, tuple)):
        geneList2 = pd.DataFrame({'Gene': geneList2})
        
    method_key = {'union': 'outer', 'intersection': 'inner'}
    return pd.merge(geneList1, geneList2, how=method_key[method], on=mergeCol)

def get_interactors(target_genes, 
                    method= 'union',
                    interactions_type = None, 
                    **kwargs):
    
    """
    Function returns list of interactors for given list of genes.
    
    Parameters
    ----------
    target_genes : list
        list of genes
    method : str
        method of merging. 'union' or 'intersection'. 'union' by default
    interactions_type : str
        None by default. Which type of interactions you want to search. It`s ready to use presets. 'all', 'physical' or 'genetic'.
        Activated only if interactions_type is not None. When activated, included metrics and EXPERIMENTAL_SYSTEM_TYPE are ignored
    **kwargs :
        parameters for fine-tuning the interactome search
    
    kwargs
    ----------
    **included_metrics** : *list*, optional.
        list of metrics you want to include in total score. Default is ['nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore']
    **taxId** : *int, optional*.
        organism id. Human = 9606. Default is 9606
    **EXPERIMENTAL_SYSTEM_TYPE** : *str, optional*.
        type of experimental system. 'physical' or 'genetic'. 'physical' by default
    **required_evidence** : *int, optional*.
        required number of evidence for interaction
        Default is 1. It means that BioGRID has 1 mention about interaction between two proteins
              
    **required_score** : *int, optional*.
        required score of interaction. From 150 to 999
    **limit** : *int or None, optional*.
        limit of interactors for each gene from geneList
    **taxId** : *int, optional*.
        organism id. Human = 9606   
         
    **mergeCol** : *str*.
        column name to merge. 'Gene' by default
    
    Returns
    -------
    pd.DataFrame
        1-column DataFrame of interactors (header = 'Gene' or choosen value of 'mergeCol')
    """
    # check correctness of kwargs
    valid_kwargs = {'EXPERIMENTAL_SYSTEM_TYPE', 'taxId', 'required_score', 'required_evidence', 'included_metrics', 'limit', 'mergeCol'}
    Check_kwargs(kwargs, valid_kwargs)
 
    EXPERIMENTAL_SYSTEM_TYPE = kwargs.get('EXPERIMENTAL_SYSTEM_TYPE', 'any')
    required_score = kwargs.get('required_score', 400)
    taxId = kwargs.get('taxId', 9606)
    required_evidence = kwargs.get('required_evidence', 1)
    limit = kwargs.get('limit', None)
    included_metrics = kwargs.get('included_metrics', ['nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore'])
    mergeCol = kwargs.get('mergeCol', 'Gene')
    
    #find interactors from biogrid and stringdb
    geneList_biogrid = get_interactors_from_biogrid(target_genes, 
                                                    EXPERIMENTAL_SYSTEM_TYPE=EXPERIMENTAL_SYSTEM_TYPE,
                                                    interactions_type=interactions_type, 
                                                    taxId=taxId,
                                                    required_evidence=required_evidence)
    geneList_string = get_interactors_from_stringdb(target_genes, 
                                                    required_score=required_score, 
                                                    included_metrics=included_metrics, 
                                                    taxId=taxId,
                                                    limit=limit,
                                                    interactions_type=interactions_type)
    # merge lists of interactors
    interactorsList = merging_interactors_stringdb_and_biogrid(geneList_biogrid, geneList_string, method=method, mergeCol=mergeCol)
    print(f'The list of interactors was successfully created.\nNumber of participants: {len(interactorsList)}')
    return interactorsList

class BiogridInteractions():
    
    def __init__(self, geneList:list, taxId=9606):
        
        """
        BiogridInteractions class conctructor.
        
        Parameters
        ----------
        geneList : list
            list of genes
        taxId : int
            organism id. Human = 9606
        """
        if not isinstance(geneList, list):
            raise Exception('Wrong type of "geneList". It must be "list"')
        self.geneList = geneList
        self.taxId = taxId
        self.interactionTable = None
        self.interactorsList = None
    
    @staticmethod    
    def _calculate_evidence(table):
        table['Evidence'] = 0
        for id in table.OFFICIAL_SYMBOL_A.unique():
            subframe = table[table.OFFICIAL_SYMBOL_A == id]
            evidence_col = subframe.apply(lambda x: subframe.OFFICIAL_SYMBOL_B[subframe.OFFICIAL_SYMBOL_B == x.OFFICIAL_SYMBOL_B].count(), axis=1)
            table.update({'Evidence': evidence_col})
        return table['Evidence']    
    
    def get_interaction_partners(self, EXPERIMENTAL_SYSTEM_TYPE = 'any', required_evidence = 1):
        
        """
        Function for getting interaction data from BioGRID database
        
        Parameters
        ----------
        EXPERIMENTAL_SYSTEM_TYPE : str, optional
            type of experimental system. Valid types are 'genetic', 'physical' or 'any'.
            Default is 'any'
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing all possible interactions between genes in geneList
        """
        
        Check_type(required_evidence, (int, float), 'required_evidence')
        Check_Value(EXPERIMENTAL_SYSTEM_TYPE, {'any', 'genetic', 'physical'}, 'EXPERIMENTAL_SYSTEM_TYPE')
        if EXPERIMENTAL_SYSTEM_TYPE == 'any':
            EXPERIMENTAL_SYSTEM_TYPE = ['genetic', 'physical']
        elif EXPERIMENTAL_SYSTEM_TYPE == 'genetic':
            EXPERIMENTAL_SYSTEM_TYPE = ['genetic']
        elif EXPERIMENTAL_SYSTEM_TYPE == 'physical':
            EXPERIMENTAL_SYSTEM_TYPE = ['physical']

        BG = BioGRID(config_filepath=conf_path)
        bg_interact = BG.interactions('json', geneList=self.geneList,
                                        includeEvidence='true',
                                        taxId=self.taxId)
        
        interactionTable = pd.read_json(StringIO(bg_interact.toDataFrame()), orient='index')
        
        interactionTable = self._target_gene_1st_col(interactionTable)
        interactionTable = interactionTable[interactionTable.OFFICIAL_SYMBOL_A.isin(self.geneList) | interactionTable.OFFICIAL_SYMBOL_B.isin(self.geneList)]
        interactionTable['Evidence'] = self._calculate_evidence(interactionTable) # calculation of evidence
        interactionTable.sort_values(by='QUANTITATION', ascending=False, inplace=True) # sorting by QUANTITATION
        interactionTable = interactionTable[(interactionTable.ORGANISM_A == self.taxId) & (interactionTable.ORGANISM_B == self.taxId)] # keep only right taxId
        interactionTable.drop_duplicates(subset=['OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B'], inplace=True) # drop duplicates
        
        # filter by EXPERIMENTAL_SYSTEM_TYPE and evidence
        interactionTable = interactionTable[interactionTable.EXPERIMENTAL_SYSTEM_TYPE.isin(EXPERIMENTAL_SYSTEM_TYPE)]
        interactionTable = interactionTable[interactionTable.Evidence >= required_evidence]
        
        self.interactionTable = interactionTable
        self.interactorsList = pd.DataFrame({'Gene': np.unique(np.concatenate((interactionTable.OFFICIAL_SYMBOL_A.unique(), \
            interactionTable.OFFICIAL_SYMBOL_B.unique())))})
                        
    def _swap_cols(self,row):
        
        """
        Function for swapping columns in a row of interaction table
        
        Parameters
        ----------
        row : pd.Series
            row of interaction table
            
        Returns
        -------
        pd.Series
            row of interaction table with swapped columns
        """
        
        A_cols = ['ENTREZ_GENE_A', 'BIOGRID_ID_A', 'SYSTEMATIC_NAME_A', 'OFFICIAL_SYMBOL_A', 'SYNONYMS_A', 'ORGANISM_A']
        B_cols = ['ENTREZ_GENE_B', 'BIOGRID_ID_B', 'SYSTEMATIC_NAME_B', 'OFFICIAL_SYMBOL_B', 'SYNONYMS_B', 'ORGANISM_B']
        
        k = {}
        for a, b in zip(A_cols, B_cols):
            k[b] = a
            k[a] = b
        new_row = row.rename(k)
        return new_row
    
    def _target_gene_1st_col(self, table):
        
        """
        The function moves target genes to column "*_A" and their interactors to column "*_B" in interaction table
        
        Parameters
        ----------
        table : pd.DataFrame
            interaction table
            
        Returns
        -------
        pd.DataFrame
            table with target genes in column A
        """

        all_cols = table.columns
        table = table.apply(lambda x: self._swap_cols(x) if x.OFFICIAL_SYMBOL_A not in self.geneList else x, axis=1)
        return table.reindex(columns=all_cols)


class StringdbInteractions():
    
    def __init__(self, geneList:list, taxId:int=9606):
                
        """
        StringdbInteractions class conctructor.
        
        Parameters
        ----------
        geneList : list
            list of genes
        taxId : int
            organism id. Human = 9606
        """
        if not isinstance(geneList, list):
            raise Exception('Wrong type of "geneList". It must be "list"')
        
        self.geneList = geneList
        self.taxId = taxId
        self.interactionTable = None
        self.interactorsList = None
        self.parent = NetworkAnalysis(data=pd.DataFrame(), child_class=True)
    
    def _compute_prior_away(self, score, prior:float=0.041):
        """
        Compute a score with a prior probability subtracted.

        If the score is below the prior probability, it will be set to the prior probability.
        The score is then divided by (1 - prior) to scale.

        Parameters
        ----------
        score : float
            The score to be adjusted
        prior : float, optional
            The prior probability to subtract. Defaults to 0.041.

        Returns
        -------
        float
            The adjusted score
        """
        return self.parent._compute_prior_away(score=score, prior=prior)
        
    
    def calculate_total_score(self, data, included_metrics=['nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore'], p:float=0.041):
        """
        Function for calculate total score
        Check list of metrix on: 
        https://string-db.org/help/faq/#the-protein-interactions-from-the-string-website-via-web-api-calls-what-do-the-score-columns-mean-for-example-nscore-fscore-tscore-etc

        Parameters
        ----------
        data : pd.DataFrame
            pd.DataFrame containing metrics
        included_metrics : list
            list of metrics you want to include in total score
        p : float
            prior probability, taken from STRING website. I think, better don`t touch this parameter

        Returns
        -------
        pd.Series
            pd.Series with total score
        """
        return self.parent.calculate_total_score(data=data, included_metrics=included_metrics, p=p) 
        
    def get_interaction_partners(self, 
                                 required_score:int=400,
                                 limit:int=None,
                                 included_metrics=['nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore'],
                                 p:float=0.041):
        
        
        """
        Function returns interaction table from stringdb for given list of genes and organism id (taxId) and
        required score of interaction.
        
        Parameters
        ----------
        required_score : int
            required score of interaction. From 150 to 999
        limit : int or None
            limit of interactors for each gene from geneList
        included_metrics : list
            list of metrics you want to include in total score
            For example, physical metrix: ['dscore', 'escore', 'ascore', 'tscore']
        p : float
            prior probability, taken from STRING website. I think, better don`t touch this parameter
        
        Returns
        -------
        None
        """
        # for some reason, if you search partners only for 1 gene, limit=10 by default, but there are no limits for several genes
        if len(self.geneList) == 1 and limit is None:
            limit = 1000000
        
        interactionTable = stringdb.get_interaction_partners(self.geneList, required_score=required_score, limit=limit)
                
        #cutoff genes with new_score less than required_score. new_score formed by included_metrix
        interactionTable['new_score'] = self.calculate_total_score(data=interactionTable.iloc[:, 6:], included_metrics=included_metrics, p=p)
        interactionTable.rename(mapper={'score': 'old_score', 'new_score': 'score'}, axis=1, inplace=True)
        interactionTable = interactionTable[interactionTable.score >= (required_score / 1000)]
        
        self.interactionTable=interactionTable
        self.interactorsList = pd.DataFrame({'Gene': np.unique(np.concatenate((self.interactionTable.preferredName_A.unique(), \
            self.interactionTable.preferredName_B.unique())))})  # create DataFrame of unique Genes. It`s strange, but faster
    
    
    


