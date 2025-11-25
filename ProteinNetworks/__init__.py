#from .STRING_enrichment import *
#from .R_requests import *
from . import R_requests, enrichment, interactions, mapping, wrappers
from .wrappers import save_table, display_df
from .enrichment import get_enrichment, EnrichmentAnalysis
from .mapping import get_mapping, STRING_mapping
from .interactions import get_interactors_from_biogrid, get_interactors_from_stringdb, \
    get_interactionsTable_from_biogrid, get_interactionsTable_from_stringdb, \
        merging_interactors_stringdb_and_biogrid, get_interactors

__version__ = "0.1.7"

# Fallback for optional heavy dependencies (matplotlib, igraph, leidenalg, umap)
try:
    from . import networks
    from .networks import NetworkAnalysis, create_graph
except ImportError as e:
    _networks_import_error = str(e)
    _networks_available = False
    
    class NetworkAnalysis:
        def __init__(self, *args, **kwargs):
            raise ImportError(
                f"NetworkAnalysis requires optional dependencies. Please install them with:\n"
                f"  pip install 'ProteinNetworks[full]'\n"
                f"or with conda:\n"
                f"  conda install matplotlib igraph leidenalg umap-learn -c conda-forge\n\n"
                f"Original error: {_networks_import_error}"
            )
    
    def create_graph(*args, **kwargs):
        raise ImportError(
            f"create_graph requires optional dependencies. Please install them with:\n"
            f"  pip install 'ProteinNetworks[full]'\n"
            f"or with conda:\n"
            f"  conda install matplotlib igraph leidenalg umap-learn -c conda-forge\n\n"
            f"Original error: {_networks_import_error}"
        )
else:
    _networks_available = True
