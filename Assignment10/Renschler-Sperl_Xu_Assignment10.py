# Author: Huiting Xu and Desirée Renschler-Sperl
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.optimize import curve_fit


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


def main(filename):
    graph = parse_csv(filename)

    num_nodes = graph.number_of_nodes()
    print("Number of nodes:", num_nodes)

    num_edges = graph.number_of_edges()
    print("Number of edges:", num_edges)

    degree = sorted([d for n, d in graph.degree()], reverse=True)
    degrees, degree_count = np.unique(degree, return_counts=True)

    log_degree=np.log10(degrees)
    log_prob=np.log10(degree_count / num_nodes)

    plt.loglog(degrees, degree_count / num_nodes, 'o')
    plt.xlabel('k (degree)')
    plt.ylabel('P(k)')
    plt.title('Node Degree Distribution')
    plt.savefig('nodedist.pdf')
    plt.show()

    def func(x, a, b):
        return a * x + b

    gamma=curve_fit(func, log_degree, log_prob)[0][0]
    print('Degree exponent γ: ',gamma)

    num_components = nx.number_connected_components(graph)
    print("Number of connected components:", num_components)

    largest_component = max(nx.connected_components(graph), key=len)
    sub_graph = graph.subgraph(largest_component)

    random_node1 = np.random.choice(list(sub_graph.nodes), size=10000, replace=True)
    random_node2 = np.random.choice(list(sub_graph.nodes), size=10000, replace=True)
    path_lengths = []
    for i, j in zip(random_node1, random_node2):
        if i!=j:
            path_lengths.append(nx.shortest_path_length(sub_graph, i, j))

    plt.hist(path_lengths, bins='auto')
    plt.xlabel('Shortest Path Length')
    plt.ylabel('Count')
    plt.title('Shortest Path Length Distribution')
    plt.savefig('shortestpathdist.pdf')


if __name__ == '__main__':
    args = create_parser()
    main(args.input)
