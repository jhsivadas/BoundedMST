from typing import List, Dict, Tuple
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class BoundedMSTProblem:
    # number of vertices and edges
    n: int
    m: int

    # bounds[i] is the maximum degree of node i in the MST
    bounds: Dict[int, int]

    # adj[i][j] = w -> there is an edge between (i, j) of weight w
    adj: Dict[int, Dict[int, int]]

    # edges in the original order they appeard in the input
    ee: List[Tuple[int, int, int]]

    @classmethod
    def from_file(cls, file):
        with open(file, 'r') as f:
            n, m = map(int, next(f).strip('\n').split(' '))
            bounds = dict()

            for i in range(n):
                bounds[i] = int(next(f).strip('\n'))

            ee: List[Tuple[int, int, int]] = []
            for _ in range(m):
                u, v, w = map(int, next(f).strip('\n').split(' '))
                u -= 1
                v -= 1
                ee.append((u, v, w))
            ret = cls(bounds, ee)

        return ret

    def add_edge(self, u, v, w):
        if u not in self.adj:
            self.adj[u] = dict()
        self.adj[u][v] = w

    def __init__(self, bounds, edges):
        self.n = len(bounds)
        self.m = len(edges)
        self.ee = edges
        self.bounds = bounds
        self.adj = dict()
        for u, v, w in edges:
            self.add_edge(u, v, w)
            self.add_edge(v, u, w)

    def score_file(self, output_file):
        try:
            with open(output_file, 'r+') as f:
                selected_edges = [int(i) for i in f]
        except:
            return f"Unable to parse output from {output_file}", 0

        # output is of length n-1 and has valid edges
        if len(selected_edges) != self.n - 1:
            return f"Too many edges in output (expected {self.n - 1}, got {len(selected_edges)})", 0

        if not all(1 <= i <= self.m for i in selected_edges):
            return f"Not all edges in range [1, {self.m}]", 0

        # check that output forms a tree
        adj = defaultdict(list)
        cost = 0
        for idx in selected_edges:
            u, v, w = self.ee[idx - 1]
            adj[u].append(v)
            adj[v].append(u)
            cost += w

        vis = set()

        def visit(idx):
            stack = [idx]
            while len(stack):
                u = stack[-1]
                stack.pop()
                if u in vis:
                    continue
                vis.add(u)
                for v in adj[u]:
                    if v not in vis:
                        stack.append(v)

        visit(0)

        if len(vis) != self.n:
            return "Output does not form a tree", 0

        # check bounds
        for i in range(self.n):
            if len(adj[i]) > self.bounds[i]:
                return f"Node {i} exceeds bound of {self.bounds[i]} (actually has degree {len(adj[i])})", 0

        return None, cost
