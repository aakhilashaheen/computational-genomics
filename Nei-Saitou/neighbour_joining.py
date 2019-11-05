from node import Node
import sys
import random
import utils

# Holds the methods to join neighbours
class Neighbour_Joining:

    #Applies nei_saitou to calculte the distances recursively and get the phylogeny tree
    def nei_saitou(self, original_ids, original_distMatrix, nodeId):
        #Using a global root holder
        global root
        root = None
        # Maintaining a map for the relationships between child and parent
        relationship_map = {}

        # keeps track of id -> index mappings for O(1) lookup
        seqId_index_map = {}
        for i, id in enumerate(original_ids):
            seqId_index_map[id] = str(i + 1)
        
        def calculate_next_neighbours(ids, distMatrix, nodeId, seqId_index_map):
            global root
            N  = len(distMatrix)

            # Base case for recursion
            if N <= 2:
                # Assigning the ids to these nodes according to the requirements
                if ids[0] in seqId_index_map:
                    tip_1 = seqId_index_map[ids[0]]
                else:
                    tip_1 = ids[0]
                if ids[1] in seqId_index_map:
                    tip_2 = seqId_index_map[ids[1]]
                else:
                    tip_2 = ids[1]
                # Adding their relationship to the tree
                relationship_map[tip_1] = (tip_2, distMatrix[0][1])
                # Since these are the last two nodes, any of these can be chosen to be the root
                root = tip_2
                return

            #First step is to get the closest two tips from the distance matrix
            mini, minj, minVal = utils.construct_Q_matrix(distMatrix)
            # get the edge lengths
            edge_i, edge_j = utils.calculate_edge_lengths(distMatrix, mini, minj)

            tip_1 = tip_2 = None
            # Setting the new ids to the tips/ internal node
            if ids[mini] in seqId_index_map:
                tip_1 = seqId_index_map[ids[mini]]
            else:
                tip_1 = ids[mini]
            if ids[minj] in seqId_index_map:
                tip_2 = seqId_index_map[ids[minj]]
            else:
                tip_2 = ids[minj]
            
            # Add them to the tree relationships
            relationship_map[tip_1] = (str(nodeId), edge_i)
            relationship_map[tip_2] = (str(nodeId), edge_j)

            # Remove the used ids and add the new nodeId
            ids = [id for id in ids if id not in (ids[mini], ids[minj])] + [str(nodeId)]

            nodeId -=1
            new_distMatrix = utils.calculate_new_distMatrix(distMatrix, mini, minj)

            # recursively call with the new distance matrix
            calculate_next_neighbours(ids, new_distMatrix, nodeId, seqId_index_map)

        # initial call to the recursive helper function with the original params
        calculate_next_neighbours(original_ids, original_distMatrix, nodeId, seqId_index_map)

        # construct tree from the relationship_map relations
        rootNode = Node(root)
        queue = [rootNode]
        # bfs for tree construction
        while len(queue) > 0:
            nodeList = []
            for node in queue:
                for child, (parent, distance) in relationship_map.items():
                    if parent == node.id:
                        childNode = Node(child)
                        node.add_child(childNode, distance)
                        nodeList.append(childNode)
            for node in nodeList:
                del relationship_map[node.id]
            queue = nodeList
        return rootNode

