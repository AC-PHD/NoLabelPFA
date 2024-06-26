from umap import UMAP
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# umap paramters
# path_original_data: string path to the original input file
# path_principal_features: string path to the txt file containing the principal features
# n_components: number of components in the embedding output
# n_neighbors: balance between local and global structures. The higher the value the bigger the neighborhoods considered for each point.
# init: initialization of embedding. Default 'random'
# random_state: seed used by random number generator
# plot: flag for visual output

def umap(path_original_data, path_principal_features="principal_features0.txt", n_components=2, n_neighbors=15, init='random', random_state=0, plot=True):

    with open(path_principal_features) as f:
        pfs = f.readlines()
    pfs = [int(x[:len(x)-2]) for x in pfs]

    X = pd.read_csv(path_original_data, sep=',', header=None)

    X_pf = X.T[pfs]

    X_embedded = UMAP(n_components=n_components, init=init,
                      random_state=random_state, n_neighbors=n_neighbors).fit_transform(X_pf).T

    np.savetxt("umap_output.csv", X_embedded, delimiter=",")

    if plot:
        plt.scatter(X_embedded[0], X_embedded[1])
        plt.xlim(min(X_embedded[0]) - 1, max(X_embedded[0]) + 1)
        plt.xticks(fontsize=16, fontweight='bold')
        plt.yticks(fontsize=16, fontweight='bold')
        filename_umap = f"umap_nc{n_components}_nn{n_neighbors}.svg"
        plt.savefig(filename_umap, dpi=600, bbox_inches='tight')
        plt.show()
