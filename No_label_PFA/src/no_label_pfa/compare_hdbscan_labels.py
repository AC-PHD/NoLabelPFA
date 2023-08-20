import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# compare_hdbscan_labels paramters
# path_comparison_labels: string path to the file containing labels to compare
# path_hdbscan_labels: string path to the file containing the labels found by hdbscan, optional
# clusters: list of hdbscan clusters to be considered in the calculation. If empty, all clusters are considered
# cmap: colormap used for plotting (default is 'viridis', other examples are 'PiYG', 'twilight' or 'tab20'. For more details check the matplotlib colomap documentation)

def compare_hdbscan_labels(path_comparison_labels, path_hdbscan_labels="hdbscan_labels.csv", clusters=[], cmap=None):
    comparison_labels = pd.read_csv(
        path_comparison_labels, sep=',', header=None).to_numpy().flatten()
    clustering = pd.read_csv(path_hdbscan_labels, sep=',',
                             header=None).to_numpy().flatten()

    if len(clusters) == 0:
        clusters = np.unique(clustering)

    cluster_bins = []

    for l in clusters:
        tmp = [int(l)]
        for cl in np.unique(comparison_labels):
            count = 0
            for i in range(len(comparison_labels)):
                if clustering[i] == l and comparison_labels[i] == cl and clustering[i] in clusters:
                    count += 1
            tmp.append(count)

        cluster_bins.append(tmp)

    columns = ['HDBSCAN Label']

    for cl in np.unique(comparison_labels):
        columns.append(str(cl))

    df = pd.DataFrame(cluster_bins, columns=columns)

    df.plot(x='HDBSCAN Label', kind='bar', stacked=True, figsize=(10,6), colormap=cmap)
    legend = plt.legend(title="Comparison Label", fancybox=True, fontsize=16, bbox_to_anchor=(1.05, 1.0), loc='upper left')
    legend.get_title().set_fontweight('bold')
    for text in legend.get_texts():
            text.set_fontweight('bold')
    plt.tight_layout()
    plt.xticks(fontsize=16, fontweight='bold')
    plt.yticks(fontsize=16, fontweight='bold')
    plt.savefig('HDB_comparison_plot.svg', dpi=600, bbox_inches='tight')
    plt.show()

