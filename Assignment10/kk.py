# Author: Huiting Xu and Desirée Renschler-Sperl
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import argparse


# create perse
def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('-i', '--input', help='csv_file')
    return p.parse_args()


def parse_csv(filename):
    graph = nx.Graph()
    with open(filename, 'r') as file:
        for row in file.readlines():
            protein1, protein2 = row.split(',')
            graph.add_edge(int(protein1), int(protein2))
    return graph


def main():
    # df = pd.read_csv('test.csv',header=None)
    graph = parse_csv('human_interactome.csv')
    # graph = nx.Graph()
    # graph.add_edges_from(df.values)

    num_nodes = graph.number_of_nodes()
    print("Number of nodes:", num_nodes)

    num_edges = graph.number_of_edges()
    print("Number of edges:", num_edges)

    degree = sorted([d for n, d in graph.degree()], reverse=True)
    degree_count = np.bincount(degree)

    plt.loglog(degree_count, 'o')
    plt.xlabel('k (degree)')
    plt.ylabel('P(k)')
    plt.title('Node Degree Distribution')
    plt.show()



    num_components = nx.number_connected_components(graph)
    print("Number of connected components:", num_components)




if __name__ == '__main__':
    args = create_parser()
    main()
