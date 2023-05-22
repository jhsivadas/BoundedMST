import numpy as np
import cvxpy as cp
from BoundedMSTProblem import BoundedMSTProblem
import heapq


def parse_input(i):
    input_file = ''
  
    if i < 10:
        input_file = f"input0{i}.txt"
    else:
        input_file = f"input{i}.txt"

    print(input_file)    
    inst = BoundedMSTProblem.from_file('inputs/' + input_file)
    N = inst.n
    M = inst.m
    bounds = inst.bounds
    degrees = []
    for val in bounds:
        degrees.append(bounds[val])
    edges = inst.ee

    return N, M, degrees, edges

def output_to_file(solution, i):
    output_file = ''
    if i < 10:
        output_file = f"output0{i}.txt"
    else:
        output_file = f"output{i}.txt"

    file_path = "outputs/" + output_file
    with open(file_path, "w") as file:
        for item in solution:
            file.write(str(item) + "\n")


def kruskal(graph):
    edges = []
    for u in graph:
        for v, w, edge in graph[u]:
            if u < v:
                edges.append((w, u, v, edge))

    heapq.heapify(edges)
    parent = {u: u for u in graph}
    rank = {u: 0 for u in graph}
    mst_weight = 0
    mst_edges = []
    num_edges = 0

    while edges and num_edges < len(graph) - 1:
        w, u, v, edge = heapq.heappop(edges)
        if find(parent, u) != find(parent, v):
            apply_union(parent, rank, u, v)
            mst_weight += w
            mst_edges.append(edge)
            num_edges += 1

    return mst_edges

def find(parent, u):
    if parent[u] != u:
        parent[u] = find(parent, parent[u])
    return parent[u]


def apply_union(parent, rank, u, v):
    u_root, v_root = find(parent, u), find(parent, v)
    if rank[u_root] < rank[v_root]:
        parent[u_root] = v_root
    elif rank[u_root] > rank[v_root]:
        parent[v_root] = u_root
    else:
        parent[v_root] = u_root
        rank[u_root] += 1


def main():

    # for u in graph:
    #         for v, w, edge in graph[u]:
        
    for i in range(1,  21):
        
        # Get Input
        N, M, degrees, edges = parse_input(i)

        graph = {u: [] for u in range(N)}

        for index, edge in enumerate(edges):
            v1, v2, weight = edge
            graph[v1].append((v2, weight, index))

        solution = kruskal(graph)
        solution = [i + 1 for i in solution]
    
        output_to_file(solution, i)



if __name__=="__main__":
    main()

