#from .STRING_enrichment import *
#from .R_requests import *
from . import R_requests, enrichment, interactions, mapping, networks, wrappers
from .wrappers import save_table, display_df
from .enrichment import get_enrichment, EnrichmentAnalysis
from .mapping import get_mapping, STRING_mapping
from .networks import NetworkAnalysis, create_graph
from .interactions import get_interactors_from_biogrid, get_interactors_from_stringdb, \
    get_interactionsTable_from_biogrid, get_interactionsTable_from_stringdb, \
        merging_interactors_stringdb_and_biogrid, get_interactors
        
__version__ = "0.1.6"
