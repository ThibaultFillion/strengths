import sys
sys.path.append("../src/")
from strengths import *

# tests ##########################################################

def test_graph() : 
    nodes = []
    nodes.append(RDGraphSpaceNode(volume = "1 m3",  environment = 1))
    nodes.append(RDGraphSpaceNode(volume = "3 km3", environment = 10))
    nodes.append(RDGraphSpaceNode(volume = "2 dm3", environment = 8))
    
    edges = []
    edges.append(RDGraphSpaceEdge(i = 0, j = 1, distance = "4 cm",  surface = "10 nm2"))
    edges.append(RDGraphSpaceEdge(i = 1, j = 2, distance = "2 mm",  surface = "9 pm2"))
    edges.append(RDGraphSpaceEdge(i = 0, j = 2, distance = "10 cm", surface = "7 m2"))
    
    graph1 = RDGraphSpace(nodes=nodes, edges=edges)
    
    gsd = {
        "type" : "graph",
        "nodes" : [
                {"volume" : "1 m3",  "environment" : 1},
                {"volume" : "3 km3", "environment" : 10},
                {"volume" : "2 dm3", "environment" : 8}
            ],
        "edges" : [
                {"nodes" : [0, 1], "distance" : "4 cm",  "surface" : "10 nm2"},
                {"nodes" : [1, 2], "distance" : "2 mm",  "surface" : "9 pm2"},
                {"nodes" : [0, 2], "distance" : "10 cm", "surface" : "7 m2"}
            ]
        }
        
    graph2 = rdspace_from_dict(gsd)
    
    for i in range(3) :
        assert graph1.nodes[i].volume == graph2.nodes[i].volume
        assert graph1.nodes[i].environment == graph2.nodes[i].environment
        
        assert graph1.edges[i].i == graph2.edges[i].i
        assert graph1.edges[i].j == graph2.edges[i].j
        assert graph1.edges[i].surface == graph2.edges[i].surface
        assert graph1.edges[i].distance == graph2.edges[i].distance
    
# run all #############################################################

# def run_all_tests() : 
test_graph()