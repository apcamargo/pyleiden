import argparse
import pathlib
import sys
import random
import igraph as ig
from collections import defaultdict
import pyleiden


def parse_cli():
    parser = argparse.ArgumentParser(
        description="pyLeiden: a CLI tool for clustering with the Leiden algorithm.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version=pyleiden.__version__
    )
    parser.add_argument(
        "INPUT",
        type=pathlib.Path,
        help="""Tabular file containing the edges of the network. The first two columns
                should be elelemnts that are connected via an edge in the graph. A third
                column may be provided with numerical values that represent the weight of
                the edge.""",
    )
    parser.add_argument(
        "OUTPUT",
        type=pathlib.Path,
        help="""Cluster membership. Each line contains all the members of a given cluster
                separated by a tab.""",
    )
    parser.add_argument(
        "-o",
        "--objective-function",
        choices=["CPM", "modularity"],
        default="CPM",
        type=str,
        help="Whether to use the Constant Potts Model (CPM) or modularity.",
    )
    parser.add_argument(
        "-r",
        "--resolution-parameter",
        default=1.0,
        type=float,
        help="""The resolution parameter to use. Higher resolutions lead to more smaller
                clusters, while lower resolutions lead to fewer larger clusters.""",
    )
    parser.add_argument(
        "-b",
        "--beta",
        default=0.01,
        type=float,
        help="""Parameter affecting the randomness in the Leiden algorithm. This affects
                only the refinement step of the algorithm.""",
    )
    parser.add_argument(
        "-n",
        "--number-iterations",
        default=2,
        type=int,
        help="""The number of iterations to iterate the Leiden algorithm. Each iteration
                may improve the partition further. Using a negative number of iterations
                will run until a stable iteration is encountered (i.e. the quality was not
                increased during that iteration).""",
    )
    parser.add_argument(
        "-s",
        "--seed",
        default=1953,
        type=int,
        help="Seed for the random number generator.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Run quietly (will not print output to stout).",
    )
    return parser


def generate_graph(input_file):
    graph = ig.Graph.Load(input_file, format="ncol")
    graph.to_undirected(mode="each")
    return graph


def cluster_graph(
    graph, objective_function, resolution_parameter, beta, number_iterations
):
    return graph.community_leiden(
        objective_function=objective_function,
        resolution=resolution_parameter,
        weights="weight" if graph.is_weighted() else None,
        beta=beta,
        n_iterations=number_iterations,
    )


def write_clusters(output_file, graph, clusters):
    cluster_members = defaultdict(list)
    for i, m in enumerate(clusters.membership):
        cluster_members[m].append(i)
    with open(output_file, "w") as fout:
        for m in sorted(
            cluster_members, key=lambda i: len(cluster_members[i]), reverse=True
        ):
            subgraph = graph.subgraph(cluster_members[m])
            members = sorted(subgraph.vs, key=lambda i: i.degree(), reverse=True)
            fout.write("\t".join([i["name"] for i in members]) + "\n")


def main():
    parser = parse_cli()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args()
    random.seed(args.seed)
    if not args.quiet:
        print("[1/3] Reading input file and building the graph.")
    graph = generate_graph(args.INPUT)
    if not args.quiet:
        print("[2/3] Clustering nodes.")
    clusters = cluster_graph(
        graph,
        args.objective_function,
        args.resolution_parameter,
        args.beta,
        args.number_iterations,
    )
    if not args.quiet:
        print("[3/3] Writing output file.")
    write_clusters(args.OUTPUT, graph, clusters)
    if not args.quiet:
        print(f"Total number of nodes: {len(graph.vs):,}")
        print(f"Total number of clusters: {len(clusters):,}")
