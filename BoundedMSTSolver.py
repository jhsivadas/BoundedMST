import numpy as np
import cvxpy as cp
from collections import defaultdict
from BoundedMSTProblem import BoundedMSTProblem
import heapq


class BoundedMST:

    def __init__(self, N, M, degrees, edges):
        # Initial variables
        self.N = N
        self.M = M
        self.degrees = degrees
        self.edges = edges 
        self.values = [edge[2] for edge in edges]

        # NxM Matrix
        self.matrix = np.zeros((N, M))
    
        for index, val in enumerate(self.edges):
            v1, v2, weight = val
            self.matrix[v1, index] += 1
            self.matrix[v2, index] += 1

        # LP variables
        self.constraints = []
        self.x = cp.Variable(self.M, boolean=True)  
        self.objective = cp.Minimize(cp.sum(self.values * self.x))

            
    def solve_bounded_mst(self, N, M, degrees, edges):
        self.constraints += [0 <= self.x, self.x <= 1]
        self.constraints += [cp.sum(self.x) == N - 1]

        for index, row in enumerate(self.matrix):
            self.constraints.append(row @ self.x <= degrees[index])
            self.constraints.append(row @ self.x >= 1)


        problem = cp.Problem(self.objective, self.constraints)
        problem.solve()

        solution = []

        if self.x.value is not None:
            solution = [i for i, value in enumerate(self.x.value) if np.isclose(value, 1)]

        return solution


    def connect_components(self, solution):
        components = self.get_connected_components(solution)
        num_components = len(components)

        count = 0
        while num_components > 1 and count < 50:

            new_edges = []
            for edge in solution:
                new_edges.append((edge, self.edges[edge]))

            adjacency_list = defaultdict(list)
            for edge, vals in new_edges:
                node1, node2, weight = vals
                adjacency_list[node1].append((node2, weight, edge))
                adjacency_list[node2].append((node1, weight, edge))


            graph2 = []
            
            for edge in solution:
                v1, v2, weight  = self.edges[edge]
                graph2.append([v1 ,v2])



            new_components = []
            new_solution = []
            for component in components:
                graph = {key: adjacency_list[key] for key in component if key in adjacency_list}
        
                new_mst, new_comps = self.kruskal2(graph)
                new_solution.append(new_mst)
                new_components.append(new_comps)

            # Count Vertexes
            new_edges = [self.edges[s] for s in solution]
            vertex_count = {}
            for edge in new_edges:
                if edge[0] not in vertex_count:
                    vertex_count[edge[0]] = 1
                else:
                    vertex_count[edge[0]] += 1
                if edge[1] not in vertex_count:
                    vertex_count[edge[1]] = 1
                else:
                    vertex_count[edge[1]] += 1

            # MST
            component_graph = {i: [] for i in range(num_components)}
            for i in range(num_components - 1):
                self.get_connected_edges(i, new_components, component_graph, vertex_count)

            

            mst = self.kruskal(component_graph, vertex_count)

            for edge in mst:
                self.constraints.append(self.x[edge] == 1)

            problem = cp.Problem(self.objective, self.constraints)
            problem.solve()
            if self.x.value is not None:
                solution = [i for i, value in enumerate(self.x.value) if np.isclose(value, 1)]

            count+=1

            components = self.get_connected_components(solution)
            num_components = len(components)
                
        
        return solution

    
    def get_connected_edges(self, index, components, component_graph, vertex_count):
        component1 = components[index]
        for j in range(index + 1, len(components)):
            component2 = components[j]
            for edge, (v1, v2, weight) in enumerate(self.edges):
                if v1 in component1 and v2 in component2:
                    component_graph[index].append((j, edge, weight))
                    component_graph[j].append((index, edge, weight))

    def kruskal2(self, graph):
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
        mst_nodes = []
        num_edges = 0


        while edges and num_edges < len(graph) - 1:
            w, u, v, edge = heapq.heappop(edges)
            if self.find(parent, u) != self.find(parent, v):
                self.apply_union(parent, rank, u, v)
                mst_weight += w
                mst_edges.append(edge)

                if u not in mst_nodes:
                    mst_nodes.append(u)
                if v not in mst_nodes:
                    mst_nodes.append(v)

                num_edges += 1

        return mst_edges, mst_nodes
    
    def kruskal(self, graph, vertex_count):
        edges = []
        for u in graph:
            for v, edge_name, edge_weight in graph[u]:
                if u < v:
                    edges.append((edge_weight, u, v, edge_name))

        heapq.heapify(edges)
        parent = {u: u for u in graph}
        rank = {u: 0 for u in graph}
        mst_weight = 0
        mst_edges = []
        num_edges = 0

   

        while edges and num_edges < len(graph) - 1:
            w, u, v, edge_name = heapq.heappop(edges)
            if (
                self.find(parent, u) != self.find(parent, v)
            ):
                self.apply_union(parent, rank, u, v)
                mst_weight += w
                mst_edges.append(edge_name)
                num_edges += 1
                vertex_count[u] += 1
                vertex_count[v] += 1

        return mst_edges


    def find(self, parent, u):
        if parent[u] != u:
            parent[u] = self.find(parent, parent[u])
        return parent[u]


    def apply_union(self, parent, rank, u, v):
        u_root, v_root = self.find(parent, u), self.find(parent, v)
        if rank[u_root] < rank[v_root]:
            parent[u_root] = v_root
        elif rank[u_root] > rank[v_root]:
            parent[v_root] = u_root
        else:
            parent[v_root] = u_root
            rank[u_root] += 1

        
    def dfs(self, node, adjacency_list, visited, component):
        visited[node] = True
        component.append(node)

        for neighbor in adjacency_list[node]:
            if not visited[neighbor]:
                self.dfs(neighbor, adjacency_list, visited, component)

    def get_connected_components(self, solution):
        adjacency_list = defaultdict(list)
     
        new_edges = []
        for edge in solution:
            new_edges.append(self.edges[edge])


        for edge in new_edges:
            node1, node2, _ = edge
            adjacency_list[node1].append(node2)
            adjacency_list[node2].append(node1)

      

        num_nodes = len(adjacency_list)
        visited = [False] * (num_nodes + 1)
        components = []

        for node in adjacency_list:
            if not visited[node]:
                component = []
                self.dfs(node, adjacency_list, visited, component)
                components.append(component)
        return components


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


def main():
        
    for i in range(1,  21):
        
        # Get Input
        N, M, degrees, edges = parse_input(i)
    
        # Solve
        solver = BoundedMST(N, M, degrees, edges)
        solution = solver.solve_bounded_mst(N, M, degrees, edges)
        solution2 = solver.connect_components(solution)
        final_solution = [i + 1 for i in solution2]

        # Output to file 
        output_to_file(final_solution, i)



if __name__=="__main__":

    main()
