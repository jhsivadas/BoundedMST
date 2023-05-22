import numpy as np
import cvxpy as cp
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
                    heapq.heappush(neighbor_heap, (-ratio, edges, n, edge_num))  # keep track of max 
                
            # print(neighbor_heap)   

            # Add connected edges / select which edges to take 
            
            for i in range((min(curr_node_bound, len(neighbor_heap)))):

                # Take heap with largest ratio in heap 
                ratio, edges, v_num, edge_num = heapq.heappop(neighbor_heap)

                ratios[v_num][0] -= 1
                ratios[v_num][1] -= 1
                self.degrees[v_num] -= 1 
                if connections[v_num] < bounds:
                    connections[v_num] += 1
                    visited.add(v_num)

                    # Check if edges
                    if ratios[v_num][0] == 0:
                        ratios[v_num][2] = 0
                    else:
                        ratios[v_num][2] = (ratios[v_num][1] / ratios[v_num][0]) /ratios[v_num][3]
                    
                    heapq.heappush(can_enter, (-ratios[v_num][2], v_num, edge_num))
        

            if can_enter == []:
                print("nope")
                break
            
            
            ratio, v_num, edge_num = heapq.heappop(can_enter)

            curr_node = v_num
            
            solution.append(edge_num)

        return [i + 1 for i in solution]
        

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

        # Output to file 
        output_to_file(solution, i)
