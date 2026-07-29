#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2026 Jiawen Zhao.
# All rights reserved.

"""
Cluster miRNAs based on pairwise relationships provided in a two-column file.
Each line should contain at least two tab-separated fields:
  <miRNA1>  <miRNA2>  ...
Empty lines or lines with fewer than two fields are ignored.
The script identifies connected components and assigns a cluster ID (CL1, CL2, …)
to each miRNA. Output is a three-column file:
  miRNA    cluster_id    intra_cluster_connectivity_ratio
where the ratio is the number of edges the miRNA has to other members of the
same cluster divided by the total number of miRNAs in that cluster.
"""

import argparse
import sys
from collections import defaultdict

class UnionFind:
    """Union-Find data structure for connected components."""
    def __init__(self, elements):
        self.parent = {e: e for e in elements}
        self.rank = {e: 0 for e in elements}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        rx = self.find(x)
        ry = self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            self.parent[rx] = ry
        elif self.rank[rx] > self.rank[ry]:
            self.parent[ry] = rx
        else:
            self.parent[ry] = rx
            self.rank[rx] += 1

def main():
    parser = argparse.ArgumentParser(
        description='Cluster miRNAs based on pairwise relationships, with connectivity ratio.')
    parser.add_argument('-i', '--input', required=True,
                        help='Input file with two columns of miRNA names (tab-separated).')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file with miRNA, cluster ID, and ratio.')
    args = parser.parse_args()

    # Read edges and collect all nodes; also build adjacency list
    nodes = set()
    edges = []
    adjacency = defaultdict(set)   # neighbor sets (undirected)

    with open(args.input, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            node1, node2 = parts[0].strip(), parts[1].strip()
            if not node1 or not node2:
                continue
            nodes.add(node1)
            nodes.add(node2)
            edges.append((node1, node2))
            adjacency[node1].add(node2)
            adjacency[node2].add(node1)

    if not nodes:
        sys.stderr.write("Warning: no valid miRNA pairs found in input.\n")
        with open(args.output, 'w') as out:
            out.write("")
        sys.exit(0)

    # Build union-find
    uf = UnionFind(nodes)
    for u, v in edges:
        uf.union(u, v)

    # Group nodes by their root
    clusters = defaultdict(list)
    for node in nodes:
        root = uf.find(node)
        clusters[root].append(node)

    # Sort clusters deterministically by the minimal miRNA name in each cluster
    sorted_clusters = []
    for root, members in clusters.items():
        sorted_clusters.append((min(members), members))
    sorted_clusters.sort(key=lambda x: x[0])

    # Write output with three columns
    with open(args.output, 'w') as out:
        for idx, (_, members) in enumerate(sorted_clusters, start=1):
            clus_id = f"CL{idx}"
            cluster_size = len(members)
            members_set = set(members)   # for quick lookup

            # For each member, count how many neighbours are in the same cluster
            for mirna in sorted(members):   # sorted alphabetically for stable output
                # Count unique neighbours inside the cluster
                count = sum(1 for neighbor in adjacency[mirna] if neighbor in members_set)
                ratio = count / cluster_size
                out.write(f"{mirna}\t{clus_id}\t{ratio:.4f}\n")

if __name__ == '__main__':
    main()
