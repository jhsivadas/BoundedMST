import numpy as np
import cvxpy as cp
from itertools import combinations
import os
import cv2
from collections import defaultdict
from BoundedMSTProblem import BoundedMSTProblem
import heapq

cyclenumber = 0
class BoundedMST:

    def __init__(self, N, M, degrees, edges):
        # Initial variables
        self.N = N
        self.M = M
        self.degrees = degrees
        self.edges = edges 
        self.values = [edge[2] for edge in edges]

        # NxM Matrix
        # self.matrix = np.zeros((N, M))
    
        # for index, val in enumerate(self.edges):
        #     v1, v2, weight = val
        #     self.matrix[v1, index] += 1
        #     self.matrix[v2, index] += 1

        # # LP variables
        # self.constraints = []
        # self.x = cp.Variable(self.M, boolean=True)  
        # self.objective = cp.Minimize(cp.sum(self.values * self.x))

    
    def get_all_neighbhors(self):

        neighbors = {u:[] for u in range(self.N)}
        for index, edge in enumerate(self.edges):
            v1, v2, weight = edge
            neighbors[v1].append((v2, index)) # add connected_node, edge_num
            neighbors[v2].append((v1, index))

        return neighbors
   
    def greedy_ratio(self):

        neighbors = self.get_all_neighbhors()

        ratios = {u: [0, 0, 0, 0] for u in range(self.N)} # declare list of edges_num, bounds, ratio, connected_edges

        # Keep track of number of edges for each vertex
        for v1, v2, weight in self.edges:
            ratios[v1][0] += 1
            ratios[v2][0] += 1
            ratios[v2][3] = weight
            ratios[v1][3] = weight

 
        # Get bounds and calculate bounds / edges ratio for each vertex
        for vertex in range(self.N):
            bounds = self.degrees[vertex]
            ratios[vertex][1] = bounds
            ratios[vertex][2] = (ratios[vertex][1] / ratios[vertex][0]) / ratios[vertex][3]

        # Get vertex with the maximum ratio
        maximum_ratio = max(ratios.items(), key=lambda x: x[1][2])
        curr_node = maximum_ratio[0]
        


        # Start traversal
        visited = set() # keeps track of vertexes visited
        visited.add(curr_node)
        can_enter = [] # keep track of which to visit next
        solution = [] # Final edges we take
        connections = [0 for i in range(self.N)]
        while len(solution) < self.N - 1:
            neighbor = neighbors[curr_node]  # list of neighbors and edge_num ot that neighbor
            curr_node_bound = self.degrees[curr_node] # integer bound (Gabe thinks error is here)

        
           # Add all neighbors to heap
            neighbor_heap = [] 
            for n, edge_num in neighbor:
                if n not in visited:
                    ratio = ratios[n][2] # gets ratio for that vertex
                    edges = ratios[n][0]
                    degs = ratios[n][1]
                    heapq.heappush(neighbor_heap, (-ratio, degs, edges, n, edge_num))  # keep track of max 
                
            # print(neighbor_heap)   

            # Add connected edges / select which edges to take 
            
            for i in range((min(curr_node_bound, len(neighbor_heap)))):

                # Take heap with largest ratio in heap 
                ratio, edges, v_num, edge_num = heapq.heappop(neighbor_heap)

                # min_values = []

                # while len(neighbor_heap) > 0: 
                #     next = neighbor_heap[0]
                #     next_ratio 
                #     if ratio == next_ratio:
                #         print("hi")
                #         break
                #     else:
                #         break 
                        

                
                

                ratios[v_num][0] -= 1
                ratios[v_num][1] -= 1
                self.degrees[v_num] -= 1 # Gabe's fix
                if connections[v_num] < bounds:
                    connections[v_num] += 1
                    visited.add(v_num)

                    # Check if edges
                    if ratios[v_num][0] == 0:
                        ratios[v_num][2] = 0
                    else:
                        ratios[v_num][2] = (ratios[v_num][1] / ratios[v_num][0]) /ratios[v_num][3]
                    
                    heapq.heappush(can_enter, (-ratios[v_num][2], ratios[degs][0], ratios[v_num][0] ,v_num, edge_num)) # push heap tiebreaker
        

            if can_enter == []:
                print("nope")
                break
            # if can_enter == []:
            #     print("hi")
            #     missing = []
            #     for i in range(self.N):
            #         if i not in visited:
            #             missing.append(i)

            #     ## Change this so that it selects the node missing that has hte highest rati0
            #     missing_max = -1
            #     missing_max_node = -1
            #     for node in missing:
            #         if missing_max < ratios[node][2]:
            #             missing_max_node = node
            #             missing_max = ratios[node][2]
                    
            #     curr_node = missing_max_node
            #     continue
            
            
            ratio, degs, edge_count, v_num, edge_num = heapq.heappop(can_enter) # grab the next vertex with the smallest ratio, then number of edges
            # print(edge_num)
            curr_node = v_num
            
            solution.append(edge_num)
    
            # print(visited)
            # print(count)
        # print(len(solution))

       




        return [i + 1 for i in solution]
            


    def solve_bounded_mst(self, N, M, degrees, edges):
        new_degrees = degrees
        bounds = max(degrees)
        print(bounds)
        for i in range(2, bounds + 1):
            for j in range(len(degrees)):
                if degrees[j] > i:
                    new_degrees[j] = i

        
            self.constraints += [0 <= self.x, self.x <= 1]
            self.constraints += [cp.sum(self.x) == N - 1]

            for index, row in enumerate(self.matrix):
                self.constraints.append(row @ self.x <= new_degrees[index])
                self.constraints.append(row @ self.x >= 1)


            problem = cp.Problem(self.objective, self.constraints)
            problem.solve()

            solution = []

            if self.x.value is not None:
                solution = [i for i, value in enumerate(self.x.value) if np.isclose(value, 1)]
                self.connect_components(solution)
                components = self.get_connected_components(solution)
                if len(components) == 1:
                    break
      
            

            # if solution == []:
            #     print("fail")

        return solution




    def connect_components(self, solution):
    
        print("solution")
        print(solution)
        components = self.get_connected_components(solution)
        num_components = len(components)

        print("___componenets___")
        print(components)

        

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
            # print('newmst')
            # print(new_mst)
        print("______new_componenets_____")
        print(new_components)
        print(len(components[0]))
        print(len(new_components[0]))
        print("___new_solution___")
        print(new_solution)

            
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
        print("____componenet_graph_____")
        print(component_graph)

        

        mst = self.kruskal(component_graph, vertex_count)
        if mst != []:
            return solution
        print("______mst_____")
        print(mst)

            # update solution 
            
        solution += mst

        print('hi')
        print(len(solution))
        print(self.N)
        print()
            # return solution
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

        # print(parent)


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
                and vertex_count[u] < self.degrees[u]
                and vertex_count[v] < self.degrees[v]
            ):
                self.apply_union(parent, rank, u, v)
                mst_weight += w
                mst_edges.append(edge_name)
                num_edges += 1
                vertex_count[u] += 1
                vertex_count[v] += 1

        return mst_edges

        # while edges and num_edges < len(graph) - 1:
        #     w, u, v, edge_name = heapq.heappop(edges)
        #     if self.find(parent, u) != self.find(parent, v):
        #         self.apply_union(parent, rank, u, v)
        #         mst_weight += w
        #         mst_edges.append(edge_name)
        #         num_edges += 1

  


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
        
    for i in range(1, 21):
        
        # Get Input
        N, M, degrees, edges = parse_input(i)
    
        # Solve
        solver = BoundedMST(N, M, degrees, edges)
        solution = solver.greedy_ratio()
        # solution = solver.solve_bounded_mst(N, M, degrees, edges)
        # solution2 = solver.connect_components(solution)
        # final_solution = [i + 1 for i in solution2]

        # Output to file 
        output_to_file(solution, i)














def findNewCycles(path, graph, cycles):
    start_node = path[0]
    next_node= None
    sub = []

    #visit each edge and each node of each edge
    for edge in graph:
        node1, node2 = edge
        if start_node in edge:
                if node1 == start_node:
                    next_node = node2
                else:
                    next_node = node1
                if not visited(next_node, path):
                        # neighbor node not on path yet
                        sub = [next_node]
                        sub.extend(path)
                        # explore extended path
                        findNewCycles(sub, graph, cycles)
                elif len(path) > 2  and next_node == path[-1]:
                        # cycle found
                        p = rotate_to_smallest(path)
                        inv = invert(p)
                        if isNew(p, cycles) and isNew(inv, cycles):
                            cycles.append(p)

def invert(path):
    return rotate_to_smallest(path[::-1])

#  rotate cycle path such that it begins with the smallest node
def rotate_to_smallest(path):
    n = path.index(min(path))
    return path[n:]+path[:n]

def isNew(path, cycles):
    return not path in cycles

def visited(node, path):
    return node in path    

def get_edges(nodes):
    edges = []
    for cycle in nodes:
        # print(cycle)
        edge = []
        for i in range(len(cycle) - 1):
            edge.append((cycle[i], cycle[i + 1]))
            # print(edge)
        edge.append((cycle[len(cycle) - 1], cycle[0]))
        edges.append(edge)
    return edges
        

def remove_max_edge(graph, cycles):
    removed_edges = []
    # print(cycles)

    for cycle in cycles:
        max_weight = float('-inf')
        max_edge = None

        for edge in cycle:
            key = edge[0]
            connected = edge[1]

            connections = graph[key]
            for connection in connections:
                weight = connection[2]
                if weight > max_weight and connection[1] not in removed_edges:
                    max_weight = weight
                    max_edge = connection[1]

        if max_edge is not None:
            removed_edges.append(max_edge)
            # print(connections)
            # print(edge)
        #     node = edge[0]
        #     for item in graph[node]:
        #         if item[0] == edge[1] and item[2] > max_weight:
        #             max_weight = item[2]
        #             max_edge = (node, item[1])

        # if max_edge is not None:
        #     graph[max_edge[0]] = [(dest, name, weight) for (dest, name, weight) in graph[max_edge[0]] if dest != max_edge[1]]
        #     removed_edges.append(max_edge[1])

    return removed_edges



if __name__=="__main__":
    # graph = [[0, 4], [0, 5], [1, 5], [2, 6], [2, 8], [2, 9], [3, 7], [7, 9], [8, 9]]
    # cycles = []   
    # for edge in graph:
    #     for node in edge:
    #         findNewCycles([node], graph, cycles)
    # print(cycles)
    # for cy in cycles:
    #     path = [str(node) for node in cy]
    #     s = ",".join(path)
        # print(s)
    main()




        # for edge in mst:
        #     self.constraints.append(self.x[edge] == 1)

        # problem = cp.Problem(self.objective, self.constraints)
        # problem.solve()
        # if self.x.value is not None:
        #     solution = [i for i, value in enumerate(self.x.value) if np.isclose(value, 1)]
        



        # Do it here 
