import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
from scipy.stats import kendalltau
from tqdm import tqdm

from typing import Any


def kt_cluster(mean_window_expression_ranks: np.ndarray, t: int = 4) -> Any:
    """For clustering points (ideally mean window expression vectors) via 1 - KT, where KT is the Kendall-Tau correlation of those vectors

    Parameters
    ----------
    mean_window_expression_ranks :
        
    t : Number of clusters to look for
        (Default value = 4)
    mean_window_expression_ranks : np.ndarray :
        
    t : int :
        (Default value = 4)
    mean_window_expression_ranks: np.ndarray :
        
    t: int :
         (Default value = 4)

    Returns
    -------

    
    """

    def kt(a: np.ndarray, b: np.ndarray) -> float:
        """

        Parameters
        ----------
        a : First vector
            
        b : Second vector
            
        a : np.ndarray :
            
        b : np.ndarray :
            
        a: np.ndarray :
            
        b: np.ndarray :
            

        Returns
        -------

        
        """
        assert len(a) == len(b)
        pbar.update(1)
        return 1 - kendalltau(a, b).correlation

    with tqdm(desc='Computing 1- Kendall-Tau Correlations') as pbar:
        Z = linkage(mean_window_expression_ranks.T,
                    metric=kt,
                    method='weighted',
                    optimal_ordering=True)
    return fcluster(Z, t=t, criterion='maxclust'), Z


def leiden_with_silhouette_score(X, leiden_nneighbors, skip_silhouette=False, leiden_iterations=10):
    """

    Parameters
    ----------
    X :
        
    leiden_nneighbors :
        
    skip_silhouette :
         (Default value = False)
    leiden_iterations :
         (Default value = 10)

    Returns
    -------

    """
    from sklearn.neighbors import kneighbors_graph
    from panopticon.utilities import get_igraph_from_adjacency
    from panopticon.utilities import import_check
    from sklearn.metrics import silhouette_score
    from collections import namedtuple

    exit_code = import_check("leidenalg",
                             'conda install -c conda-forge leidenalg')
    if exit_code != 0:
        return
    import leidenalg
    A = kneighbors_graph(X,
                         leiden_nneighbors,
                         mode='connectivity',
                         include_self=True,
                         metric='cosine')
    ig = get_igraph_from_adjacency(A)
    part = leidenalg.find_partition(
        ig,
        leidenalg.RBConfigurationVertexPartition,
        n_iterations=leiden_iterations,
        seed=17)
    clustering = part.membership
    if skip_silhouette:
        score = None
    else:
        score = silhouette_score(
            X,
            clustering,
            metric='cosine',
        )
    
    leiden_silhouette_output = namedtuple("LeidenSilhouetteOutput",
                                          "score nneighbors clustering")

    return leiden_silhouette_output(score, leiden_nneighbors, clustering)


def silhouette_optimized_leiden(X,
                                min_neighbors=2,
                                initial_intermediate=128,
                                max_neighbors=1024,
                                verbose=True):
    """

    Parameters
    ----------
    X :
        
    min_neighbors :
         (Default value = 2)
    initial_intermediate :
         (Default value = 128)
    max_neighbors :
         (Default value = 1024)
    verbose :
         (Default value = True)

    Returns
    -------

    """

    from collections import namedtuple
    from panopticon.clustering import leiden_with_silhouette_score


    nneighbors_to_silhouette_score = {}
    for leiden_nneighbors in [
            min_neighbors, initial_intermediate, max_neighbors
    ]:  # 2, 128, 1024 ??
        nneighbors_to_silhouette_score[
            leiden_nneighbors] = leiden_with_silhouette_score(
                X, leiden_nneighbors).score
        if verbose:
            print("Silhoutte score with", leiden_nneighbors, "neighbors: ",
                  nneighbors_to_silhouette_score[leiden_nneighbors])

    converged = False
    lbound_nneighbors = min_neighbors
    rbound_nneighbors = max_neighbors

    lbound_score = nneighbors_to_silhouette_score[lbound_nneighbors]
    rbound_score = nneighbors_to_silhouette_score[rbound_nneighbors]
    verbose = True
    while not converged:

        top_score = np.max([
            nneighbors_to_silhouette_score[x]
            for x in nneighbors_to_silhouette_score.keys()
            if x not in [lbound_nneighbors, rbound_nneighbors]
        ])
        silhouette_score_to_nneighbors = {
            y: x
            for x, y in nneighbors_to_silhouette_score.items()
        }
        top_score_nneighbors = silhouette_score_to_nneighbors[top_score]
        midpoint_nneighbors = np.round(
            np.mean([lbound_nneighbors, rbound_nneighbors])).astype(int)
        if midpoint_nneighbors in nneighbors_to_silhouette_score.keys():
            random_shift = np.random.choice([-1, 1], p=[0.5, 0.5])
            if midpoint_nneighbors + random_shift in nneighbors_to_silhouette_score.keys():
                midpoint_nneighbors -= random_shift
            else:
                midpoint_nneighbors += random_shift
                if midpoint_nneighbors - random_shift in nneighbors_to_silhouette_score.keys():
                    converged = True
        if np.abs(rbound_nneighbors - lbound_nneighbors) <= 2:
            converged = True
            break
        nneighbors_to_silhouette_score[
            midpoint_nneighbors] = leiden_with_silhouette_score(
                X, midpoint_nneighbors).score
        if verbose:
            print("Silhoutte score with", midpoint_nneighbors, "neighbors: ",
                  nneighbors_to_silhouette_score[midpoint_nneighbors])

        midpoint_score = nneighbors_to_silhouette_score[midpoint_nneighbors]
        if midpoint_score < top_score:
            if midpoint_nneighbors < top_score_nneighbors:
                lbound_nneighbors = midpoint_nneighbors
            else:
                rbound_nneighbors = midpoint_nneighbors
        else:
            print(midpoint_nneighbors, top_score_nneighbors)
            if midpoint_nneighbors < top_score_nneighbors:
                rbound_nneighbors = top_score_nneighbors
            else:
                lbound_nneighbors = top_score_nneighbors
            top_score = midpoint_score
        if verbose:
            print('Current nneighbors bounds: ', lbound_nneighbors,
                  rbound_nneighbors)

    top_silhouette_output = leiden_with_silhouette_score(X, top_score_nneighbors)
    silhouette_optimized_leiden_output = namedtuple(
        "SilhouetteOptimizedLeidenOutput",
        "score nneighbors clustering nneighbors_to_silhouette_score")
    return silhouette_optimized_leiden_output(top_silhouette_output.score,
                                              top_silhouette_output.nneighbors,
                                              top_silhouette_output.clustering,
                                              nneighbors_to_silhouette_score)
