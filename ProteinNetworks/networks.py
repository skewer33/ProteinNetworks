import ProteinNetworks as PN
import pandas as pd
import stringdb
import networkx as nx
import leidenalg as la
import igraph as ig

from umap import UMAP
import matplotlib.pyplot as plt
from math import log
from .wrappers import display_df, Check_Value, titler, Check_kwargs
from .mapping import get_mapping

INTERACTIONS_TYPE = {'STRINGdb':{
                        'physical': {'included_metrics': ['dscore', 'escore', 'tscore']}, 
                        'genetic': {'included_metrics': ['nscore', 'fscore', 'pscore', 'ascore']},
                        'all':{'included_metrics': ['nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore']}},
                    }


def create_graph(geneList, interactions_type=None, **kwargs):
    """
    Function to create graph from gene list. It downloads interactions from STRING database and then create graph.
    For each gene, it creates node in graph. Then it creates edges between nodes if there are interactions between this genes.
    
    Parameters
    ----------
    geneList : pd.DataFrame
        pd.DataFrame contains ids in *protein_id_type* column
    interactions_type : str
        None by default. Which type of interactions you want to search. It`s ready to use presets. 'all', 'physical' or 'genetic'.
        Activated only if interactions_type is not None. When activated, included metrics is ignored    
    
    kwargs
    ----------
    **required_score** : *int.*
        Required score of interaction. From 150 to 999
    **taxId** : *int, optional.*
        Organism id. Human = 9606
    **included_metrics** : *list*, optional.
        List of metrics you want to include in total score. Default is ['nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore']
    **base_degree** : *str or int or float, optional.*
        Base of degree in weights(x) function (see *get_edge_list()* function)
    **neg_exponent** : *str or int or float, optional.*
        Exponent of degree in weights(x) function (see *get_edge_list()* function)
        
    Returns
    ----------
    **NetworkAnalysis() object** contains graph with nodes from 'geneList'
    """
    
    # check correctness of kwargs
    valid_kwargs = {'taxId', 'required_score', 'included_metrics', 'base_degree', 'neg_exponent'}
    Check_kwargs(kwargs, valid_kwargs)
 
    required_score = kwargs.get('required_score', 400)
    taxId = kwargs.get('taxId', 9606)
    included_metrics = kwargs.get('included_metrics', ['nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore'])
    
    mapped_genes = get_mapping(geneList)
    network_obj = NetworkAnalysis(mapped_genes)
    network_obj.get_proteins_network(identifiers=network_obj.orig_data.stringId, 
                                     required_score=required_score, 
                                     included_metrics = included_metrics, 
                                     species=taxId,
                                     interactions_type=interactions_type)
    network_obj.get_edge_list(network_obj.proteins_network, id_type='GeneName')
    network_obj.create_graph_from_edge_list(network_obj.edge_list)
    return network_obj

    
class NetworkAnalysis():
    
    id_type_converter = {'stringId': 'stringId', 'GeneName': 'preferredName'}
    def __init__(self, data, child_class=False):
        self.orig_data = data
        self.data = data.copy()
        if not child_class: # this class is parent for StringdbInteractions class (chech interaction.py). Bad solution, but it works    
            self.data.index = list(self.data.preferredName)
        self.edge_list = None
        self.graph = None
        self.id_type = None
        
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
        
        if score < prior: score = prior
        score_no_prior = (score - prior) / (1 - prior)

        return score_no_prior
    
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

        metrics_list = {'nscore', 'fscore', 'pscore', 'dscore', 'escore', 'ascore', 'tscore'}
        for col in data.columns:
            Check_Value(col, metrics_list.union(['score']), ' ', message = f'Wrong column in data variable! \
                                            It must contain only {metrics_list.union(["score"])} columns')
        
        if type(included_metrics) not in [list]:
            raise Exception('Error of "included_metrics" variable. You must specify list of metrics you want to include in total score')
        for i in included_metrics:    
            Check_Value(i, metrics_list, ' ', message = f'Wrong element in "included_metrics" list! \
                                            It must contain only {metrics_list} metrics')

        cur_data = data[included_metrics]
        
        # calculate total score (https://string-db.org/help/faq/#how-are-the-scores-computed)
        s_tot_nop = cur_data.apply(lambda y: 1 - y.apply(lambda x: 1-self._compute_prior_away(x)).product(axis=0), axis=1)
        s_tot = s_tot_nop + p * (1 - s_tot_nop)

        return s_tot 
    
    def get_proteins_network(self, identifiers, 
                species=9606, 
                required_score=400,
                included_metrics=['dscore', 'escore', 'ascore', 'tscore'], 
                p:float=0.041,
                interactions_type=None):
    
        """
        Function for get proteins network

        Parameters
        ----------
        identifiers : list
            list of identifiers
        species : int
            species
        required_score : int
            required score
        included_metrics : list
            list of metrics you want to include in total score
        p : float
            prior probability, taken from STRING website. I think, better don`t touch this parameter

        Returns
        -------
        pd.DataFrame
            pd.DataFrame with network
        """
        if interactions_type is not None:
            Check_Value(interactions_type, {'genetic', 'physical', 'all'}, 'interactions_type')
            included_metrics = INTERACTIONS_TYPE['STRINGdb'][interactions_type]['included_metrics'] 
        
        proteins_network = stringdb.get_network(identifiers, species=species, required_score=required_score)
        # calculate total score taking into account "included_metrics"
        proteins_network['new_score'] = self.calculate_total_score(proteins_network.iloc[:, 6:], included_metrics, p)
        proteins_network.rename(mapper={'score': 'old_score', 'new_score': 'score'}, axis=1, inplace=True)

        proteins_network = proteins_network[proteins_network.score >= (required_score / 1000)]
        self.proteins_network = proteins_network
        #self.data['node'] = identifiers
        self.data['node'] = self.data['preferredName']

        return proteins_network
    
    
    # def get_edge_list(self, network, base_degree='x', neg_exponent=2, id_type='stringId'):
    #     """
    #     Function for getting adjacency list from STRING network
    #     weights of edges are calculated as a function of 'score'= x: weight(x) = pow(base_degree, -neg_exponent)

    #     Weights of edges can be represented as an inverse power function 'x**(-a)', where 'a' is a real positive number \
    #                     or exponential function 'a**(-x)', where 'a' is a real positive number


    #     Example: get_edge_list(network, base_degree='x', neg_exponent=2) -> weights(x) = x**(-2)
    #                 get_edge_list(network, base_degree= np.e, neg_exponent=x) -> weights(x) = np.e**(x)
        
        
    #     Parameters
    #     ----------
    #     network : pd.DataFrame
    #         pd.DataFrame with STRING network
    #     base_degree : str or int or float
    #         base of degree in weights(x) function
    #     neg_exponent : str or int or float
    #         exponent of degree in weights(x) function
    #     id_type : str
    #         'stringId' or 'GeneName'
        
    #     Returns
    #     -------
    #     pd.DataFrame
    #         pd.DataFrame with adjacency list
    #     """

    #     # check parameters
    #     if not ((isinstance(base_degree, str) and base_degree == 'x' and isinstance(neg_exponent, (int, float)) and neg_exponent > 0) \
    #         or (isinstance(neg_exponent, str) and neg_exponent == 'x' and isinstance(base_degree, (int, float)) and base_degree > 0)):
            
    #         raise ValueError('Wrong parameters: "base_degree" must be "x" or real positive number, \
    #     "neg_exponent" must be real positive number or "x" respectively')

    #     Check_Value(id_type, {'stringId', 'GeneName'}, 'id_type')
        
    #     self.id_type = self.id_type_converter[id_type]
    #     if self.id_type == 'stringId': names = ['stringId_A', 'stringId_B']
    #     else: names = ['preferredName_A', 'preferredName_B']

    #     edge_list = pd.DataFrame()
    #     edge_list[['A', 'B']] = network[names]
    #     #prevent SettingWithCopyWarning message from appearing
    #     pd.options.mode.chained_assignment = None

    #     edge_list['weight'] = network['score'].apply(lambda x: pow(eval(str(base_degree)), eval(str('-') + str(neg_exponent))))
    #     self.edge_list = edge_list
        
    #     return edge_list

    def get_edge_list(self, network, id_type='stringId'):
        """
        Function for getting adjacency list from STRING network
        weights of edges are calculated as a function of 'score'= x: weight(x) = pow(base_degree, -neg_exponent)

        Weights of edges can be represented as an inverse power function 'x**(-a)', where 'a' is a real positive number \
                        or exponential function 'a**(-x)', where 'a' is a real positive number


        Example: get_edge_list(network, base_degree='x', neg_exponent=2) -> weights(x) = x**(-2)
                    get_edge_list(network, base_degree= np.e, neg_exponent=x) -> weights(x) = np.e**(x)
        
        
        Parameters
        ----------
        network : pd.DataFrame
            pd.DataFrame with STRING network
        base_degree : str or int or float
            base of degree in weights(x) function
        neg_exponent : str or int or float
            exponent of degree in weights(x) function
        id_type : str
            'stringId' or 'GeneName'
        
        Returns
        -------
        pd.DataFrame
            pd.DataFrame with adjacency list
        """

        Check_Value(id_type, {'stringId', 'GeneName'}, 'id_type')
        
        self.id_type = self.id_type_converter[id_type]
        if self.id_type == 'stringId': names = ['stringId_A', 'stringId_B']
        else: names = ['preferredName_A', 'preferredName_B']

        edge_list = pd.DataFrame()
        edge_list[['A', 'B']] = network[names]
        #prevent SettingWithCopyWarning message from appearing
        pd.options.mode.chained_assignment = None

        edge_list['weight'] = self.weight_modify(network['score'])
        self.edge_list = edge_list
        
        return edge_list
    
    @staticmethod
    def weight_modify(series, mode='None'):
        if mode == 'None':
            return series
        if mode == 'square':
            return series.apply(lambda x: pow(x, 2))
    
    
    def create_graph_from_edge_list(self, edge_list, weighted:bool=True):
        """
        Create graph from edge list

        Parameters
        ----------
        edge_list : pd.DataFrame
            pd.DataFrame with edges: 1st column - source vertex, 2nd column - target vertex, 3rd column - weight
        weighted : bool
            whether to use weighted edges

        Returns
        -------
        G : networkx.Graph
            networkx graph object
        """
        # Create an empty graph
        G = nx.Graph()

        # Add vertices to the graph
        vertices = set()
        for edge in edge_list.itertuples(index=False):
            vertices.add(edge[0])
            vertices.add(edge[1])
        G.add_nodes_from(vertices)

        # Add edges to the graph
        if weighted:
            for edge in edge_list.itertuples(index=False):
                G.add_edge(edge[0], edge[1], weight=edge[2])
        else: 
            for edge in edge_list.itertuples(index=False):
                G.add_edge(edge[0], edge[1], weight=1)
        
        self.graph = G
        
        return G
    
    def count_interactions(self):
        count_degree = pd.DataFrame.from_dict(dict(self.graph.degree), orient='index', columns=['degree']).reset_index(names=self.id_type)
        self.data = pd.merge(self.data, count_degree, how='outer', on=self.id_type)
        self.data.index = list(self.data.preferredName)
        count_degree = count_degree.rename(columns={'preferredName': 'protein', 'degree': '#interactions'})
        return count_degree.sort_values(by='#interactions', ascending=False)
    
    def clustering(self, method:str='louvain', weight:str='weight'):
        print(2)
        
        """
        Function for clustering proteins in graph using several methods.

        Parameters
        ----------
        method (str): Method for clustering. Available methods: 'louvain', 'infomap', 'edge_betweenness', 'fast_greedy'.
                      Default is 'louvain'.
        weight (str): Name of column in graph with weights. Default is 'weight'.
        
        Returns
        -------
        list of sets: List of clusters. Each cluster is set of protein IDs.
        """
        
        if method == 'louvain':
            clusters = nx.community.louvain_communities(self.graph, weight=weight)
        elif method == 'k_clique':
            clusters = nx.community.k_clique_communities(self.graph, k=4)
        elif method == 'greedy_modularity':
            clusters = nx.community.greedy_modularity_communities(self.graph, weight=weight)
        elif method == 'leiden':
            g = ig.Graph.TupleList([(u, v, float(weight)) for u, v, weight in self.graph.edges(data='weight')], weights=True)
            clusters = la.find_partition(g, la.ModularityVertexPartition, weights='weight')
        else: raise Exception('Wrong method name. Use: "louvain", "k_clique", "greedy_modularity", "leiden')
        
        self.graph.clusters = clusters
        
        self.data['cluster'] = 0
        for i, c in enumerate(clusters):
            #self.data.loc[self.data[self.id_type].isin(clusters[i]), 'cluster'] = i + 1
            self.data.loc[self.data[self.id_type].isin(c), 'cluster'] = i + 1
            
            
        return clusters
    
    @staticmethod
    def _umap_layout(graph):
        dist_matrix = nx.to_numpy_array(graph)
        reducer = UMAP(n_components=2)
        umap_dist = reducer.fit_transform(dist_matrix) # position
        pos = {}
        for i, node in enumerate(graph.nodes()):
            pos[node] = umap_dist[i]
        return pos
    
    @staticmethod
    def _spectral_layout(graph):
        return nx.spring_layout(graph)
    
    @staticmethod
    def _spring_layout(graph):
        if nx.number_connected_components(graph) > 1:  
            return nx.spring_layout(graph, k=0.1)
        else:
            return nx.spring_layout(graph)
    
    def define_node_position(self, method='spring'):
        Check_Value(method, {'spring', 'umap', 'spectral'}, 'method')
        d_method = {'spring': self._spring_layout,
                'spectral': self._spectral_layout,
                'umap': self._umap_layout}
        try:
            self.graph.pos[method] = d_method[method](self.graph)
        except:
            self.graph.pos = {}
            self.graph.pos[method] = d_method[method](self.graph)
        
    def draw_graph(self, pos=None, clusterList=None, save=False, **kwargs):
        """
        Function for fast graph visualization

        Parameters
        ----------
        graph: networkx.Graph
            Graph object
        clusterList: pd.Series object
            List of clusters from NetworkAnalysis.data.cluster object. Contains indexes of nodes in the graph
        **kwargs: dict
            valid kwargs:
    
            *figsize*: tuple (width, height). Size of the figure
    
            *node_size*: int. Size of the nodes
    
            *node_color*: string. Color of node. Default: '#B897D5'
    
            *font_size*: int. Size of the font. Default: 5
    
            *font_color*: string. Color of the font. Default: '#222222'
    
            *edge_color*: string. Color of the edge. Default: 'lightgray'
            
            *palette*: string. Choose one of this color palettes for coloring clusters: https://matplotlib.org/stable/users/explain/colors/colormaps.html
            
            *view_labels*: bool. True = show labels (node names)
            
            *dpi*: int. Work with 'save'=True. Dpi of saving image
            
        Returns
        -----------
        None
        
        """
        valid_kwargs = {'figsize','node_size', 'node_color', 'node_edge_color', 'edge_color', 'font_size', 'font_color', 'view_labels', 'palette', 'pos', 'dpi'}
        Check_kwargs(kwargs, valid_kwargs)
    
        figsize = kwargs.get('figsize', (10, 10))
        node_size = kwargs.get('node_size', 50)
        unicolor_node_color = kwargs.get('node_color', '#B897D5')
        font_size = kwargs.get('font_size', 5)
        font_color = kwargs.get('font_color', '#222222')
        edge_color = kwargs.get('edge_color', 'lightgray')
        view_labels = kwargs.get('view_labels', True)
        palette = kwargs.get('palette', 'rainbow')
        dpi = kwargs.get('dpi', 400)
    
        
        if clusterList is not None:
            clusters = list(map(lambda x: clusterList[clusterList.index == x].iloc[0], self.graph.nodes()))
            cmap = plt.get_cmap(palette)  # Choose a color palette
            num_clusters = len(clusterList.unique())
            node_colors = [cmap(i / num_clusters) for i in clusters]
        else:
            clusters = [1]*len(self.graph.nodes())
            node_colors = [unicolor_node_color]*len(self.graph.nodes())

        fig = plt.figure(figsize=figsize)
        
        if pos is None:
            pos = self._spring_layout(self.graph)
    
        nx.draw_networkx_nodes(self.graph, 
                                pos, 
                                node_size=node_size, 
                                node_color=node_colors)

        nx.draw_networkx_edges(self.graph, pos, alpha=0.5, edge_color=edge_color)
        
        if view_labels:
            for node in self.graph.nodes:
                plt.text(pos[node][0], pos[node][1], 
                        node, 
                        ha='center', va='center',
                        fontsize=font_size,
                        color=font_color,
                        fontweight='bold')

        plt.axis('off')
        if save:
            plt.savefig('graph.png', dpi=dpi)
        plt.show()
        
    def visualize_graph(self, pos_method='spring', view_clusters=False, clustering_method='louvain', **kwargs):
        
        """
        Function for graph visualization
        
        Parameters
        ----------
        pos_method: str
            Method for node positioning. Choose one of: 'spring', 'umap', 'spectral'
        view_clusters: bool
            If True, first, clustering will be performed with 'clustering_method' method and then cluster colors will be used for visualization
        clustering_method: str
            Choose one of the following clustering methods: 'louvain', 'k_clique', 'greedy_modularity'
        **kwargs: dict
            valid kwargs for 'draw_graph' method
        
        Returns
        -------
        None
        """
        try:
            pos = self.graph.pos[pos_method]
        except:
            self.define_node_position(method=pos_method)
            pos = self.graph.pos[pos_method]
        if view_clusters:
            self.clustering(method=clustering_method)
            self.draw_graph(pos=pos, clusterList=self.data.cluster, **kwargs)
        else:
            self.draw_graph(pos=pos, **kwargs)