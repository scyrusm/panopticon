import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
from scipy.stats import kendalltau
from tqdm import tqdm

from typing import Any


def kt_cluster(mean_window_expression_ranks: np.ndarray, t: int = 4) -> Any:
    """
    For clustering points (ideally mean window expression vectors) via 1 - KT, where KT is the Kendall-Tau correlation of those vectors

    Parameters
    ----------
    mean_window_expression_ranks :

    t : Number of clusters to look for
         (Default value = 4)

    Returns
    -------

    """

    def kt(a: np.ndarray, b: np.ndarray) -> float:
        """

        Parameters
        ----------
        a : First vector

        b :  Second vector


        Returns
        -------

        """
        assert len(a) == len(b)
        pbar.update(1)
        return 1 - kendalltau(a, b).correlation

    with tqdm(desc='Computing 1- Kendall-Tau Correlations') as pbar:
        Z = linkage(
            mean_window_expression_ranks.T,
            metric=kt,
            method='weighted',
            optimal_ordering=True)
    return fcluster(Z, t=t, criterion='maxclust'), Z
