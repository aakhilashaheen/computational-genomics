from node import Node
import sys
import random
import utils
# implemetnation of the nei saitou algorithm.
# Takes in the original ids and the distance matrix
# need original id ordering
class Neighbour_Joining:
    def nei_saitou(original_ids, original_distance_matrix, nodeId):
        # global root value and nodeId
        global root
        # nodeId starts at 120 and decrements
        root = None
        # keeps track of child parent relationships for constructing the tree later
        child_parent = {}
        # keeps track of id -> index mappings for O(1) lookup
        id_index_map = {}
        for i, id in enumerate(original_ids):
            id_index_map[id] = str(i + 1)

        # recursive helper function that computes Q matrix, q values,
        # and the new distance matrix
        def helper(ids, distance_matrix, nodeId):
            global root
            N  = len(distance_matrix)

            # ending recursive case when the matrix is size 2x2
            if N <= 2:
                # add the relationship into our dictionary and also set the root
                node_1 = id_index_map[ids[0]] if ids[0] in id_index_map else ids[0]
                node_2 = id_index_map[ids[1]] if ids[1] in id_index_map else ids[1]
                child_parent[node_1] = (node_2, distance_matrix[0][1])
                root = node_2
                return

            # default q matrix initially all 0 (NxN)
            Q_matrix = [[0] * N for _ in range(N)]

            # variables to keep track of the minimum indices in the q matrix
            min_i = min_j = 0
            minVal = float('inf')

            for i in range(N):
                for j in range(N):
                    if i == j:
                        continue
                    # based on the equation in the book
                    Q_matrix[i][j] = (N - 2) * distance_matrix[i][j] - sum(distance_matrix[i]) - sum(distance_matrix[j])
                    # keep track of the minimum value and indices in the Q matrix
                    if Q_matrix[i][j] < minVal:
                        min_i = i
                        min_j = j
                        minVal = Q_matrix[i][j]

            # u is the new node
            # distance from min_i node to u
            # Based on equation from the book
            q_i_u = 1 / 2.0 * distance_matrix[min_i][min_j] + 1 / (2.0 * (N - 2)) * (sum(distance_matrix[min_i]) - sum(distance_matrix[min_j]))
            # distance from min_j node to u
            q_j_u = distance_matrix[min_i][min_j] - q_i_u

            # store child parent relations
            # we want to get the node id. ranging from 1-120
            child_1 = id_index_map[ids[min_i]] if ids[min_i] in id_index_map else ids[min_i]
            child_2 = id_index_map[ids[min_j]] if ids[min_j] in id_index_map else ids[min_j]
            child_parent[child_1] = (str(nodeId), q_i_u)
            child_parent[child_2] = (str(nodeId), q_j_u)
            # remove the used ids from the id list and also add the new nodeid
            ids = [e for e in ids if e not in (ids[min_i], ids[min_j])] + [str(nodeId)]
            # atomically decrement the global nodeId
            nodeId -= 1

            # create a tmp distance matrix with new node at the end (index N)
            # used to calculate distances to the new node and then used to create the new smaller
            # distance matrix
            tmp_matrix = [[0] * (N + 1) for _ in range(N + 1)]
            # copy everything over
            for i in range(N):
                for j in range(N):
                    tmp_matrix[i][j] = distance_matrix[i][j]
            # update the distances to the new node
            for k in range(N):
                # based on equations in the book
                tmp_matrix[N][k] = 1 / 2.0 * (distance_matrix[min_i][k] + distance_matrix[min_j][k] - distance_matrix[min_i][min_j])
                tmp_matrix[k][N] = tmp_matrix[N][k]

            # shrink the tmp matrix to create the new distance matrix
            new_distance_matrix = [[0] * (N - 1) for _ in range(N - 1)]

            # where to move values from tmp values to the new_distance_matrix
            i_pointer = j_pointer = 0
            for i in range(N + 1):
                # if it is the minimum indices, skip to the next iteration
                if i == min_i or i == min_j:
                    continue
                j_pointer = 0
                for j in range(N + 1):
                    if j == min_i or j == min_j:
                        continue
                    # set the corresponding values in the new distance matrix
                    new_distance_matrix[i_pointer][j_pointer] = tmp_matrix[i][j]
                    j_pointer += 1
                i_pointer += 1

            # recursively call with the new distance matrix
            helper(ids, new_distance_matrix, nodeId)

        # initial call to the recursive helper function with the original params
        helper(original_ids, original_distance_matrix, nodeId)

        # construct tree from the child_parent relations
        rootNode = Node(root)
        # queue for breadth first search
        queue = [rootNode]

        # bfs to construct the tree
        # will end when all nodes have been traversed
        while len(queue) > 0:
            tmp = []
            for node in queue:
                print(node.id)
                # for every parent node taht matches the current node, create
                # a new child node and add it to the queue and add it as a child
                for child, (parent, distance) in child_parent.items():
                    if parent == node.id:
                        childNode = Node(child)
                        node.add_child(childNode, distance)
                        tmp.append(childNode)
            # delete the node in the dictionary so we dont reuse it
            for node in tmp:
                del child_parent[node.id]
            # the queue is set to the next level
            queue = tmp

        # return the root of the tree
        return rootNode
