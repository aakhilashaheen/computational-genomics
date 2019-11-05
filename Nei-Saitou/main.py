#!/usr/bin/env python
from node import Node
import sys
import random
import neighbour_joining
import utils


# uses preorder traversal to traverse the tree and generate the edges file
def write_edge_file(root):
    traversal = []
    # recursive preorder traversal
    def preorder_traversal(node):
        # Base case for recursion
        if  node is None or  node.children is None:
            return
        for child, distance in node.children.items():
            # traverse root first
            traversal.append((node.id, child.id, distance))
            # then traverse children recursively
            preorder_traversal(child)
    preorder_traversal(root)
    # Finall write to the edges file
    with open('edges.txt', 'w') as f:
        for (parent, child, distance) in traversal:
            f.write(parent + '\t' + child + '\t' + str(distance) + '\n')


# Using post order traversal of the tree for newick file
def write_newick_file(seqIds, root):
    def postorder_traversal(node):
        if  node.children is None:
            if 1 <= int(node.id) <= 61:
                return seqIds[int(node.id) - 1]
        traversal = []
        # recursively traverse the children first
        for child, distance in node.children.items():
            traversal.append(postorder_traversal(child) + ':' + str(distance))
        result = '(' + ','.join(traversal) + ')'
        return result
    (parentNode, pN), (child1, c1), (child2, c2) = root.children.items()
    # Following the format mentioned. Performing postorder on the children
    prefix = '(' + postorder_traversal(child1) + ':' + str(c1) + ','+ postorder_traversal(child2) + ':' + str(c2) + ')'
    entireOrder =  '(' + prefix + ':' + str(pN) + ');'
    with open('tree.txt', 'w') as f:
        f.write(entireOrder)


# Function that reads the given fasta file and gets the sequences.
def read_fasta_file(filename):
    sequences = {}
    ids = []
    with open(filename) as f:
        lines = f.read().splitlines()
    curr = None
    for index, line in enumerate(lines):
        if index % 2 == 0:
            curr = line[1:]
            ids.append(curr)
        else:
            sequences[curr] = line
    return ids, sequences


# get the leaves under this node using dfs
# returns a set of node.ids
def get_leaves(node):
    leaves = set()
    def dfs(node):
        global count
        if not node or not node.children:
            leaves.add(node.id)
        else:
            for child in node.children.keys():
                dfs(child)
    dfs(node)
    return leaves


# return the dictionary of id -> set of leaves under that node
# using bfs with a queue
def get_partitions(root):
    queue = [root]
    partition = {}
    while len(queue):
        tmp = []
        for node in queue:
            # this is a leaf
            if not node.children:
                continue
            leaves = get_leaves(node)
            partition[node.id] = leaves
            for child in node.children.keys():
                tmp.append(child)
        queue = tmp
    return partition


# get the dfs order of the tree
def get_id_order(root):
    order = []
    def dfs(node):
        if not node or not node.children:
            return
        else:
            order.append(node.id)
            for child in node.children.keys():
                dfs(child)
    dfs(root)
    return order


# do the 100 bootstrap sampling
def bootstrap(original_root, ids, id_sequences):
    # original order of teh ids in the original tree (dfs)
    original_order = get_id_order(original_root)
    original_partitions = get_partitions(original_root)
    partition_count = [0] * 59
    for _ in range(100):
        # new bootstrap sequence
        bootstrap_sequences = {}
        # the columns to swap
        all_indices = [i for i in range(len(id_sequences['152801']))]
        indices = [random.choice(all_indices) for _ in range(len(all_indices))]
        for id in ids:
            # choose randomly the to reconstruct the sequence
            new_sequence = ''
            for index in indices:
                new_sequence += id_sequences[id][index]
            bootstrap_sequences[id] = new_sequence
        # do the nei_saitou and get a new tree
        distMatrix = utils.get_distMatrix(ids, bootstrap_sequences)
        seqCounter = 120
        Neighbour_Joining_instance = neighbour_joining.Neighbour_Joining()
        root = Neighbour_Joining_instance.nei_saitou(ids, distMatrix, seqCounter)
        # get the dictionary of partitions and ids under those partitions
        partitions = get_partitions(root)

        # compare partitions to see which trees make the same partition as the original
        for index, id in enumerate(original_order):
            if original_partitions[id] == partitions[id]:
                partition_count[index] += 1
    percentages = [count / 100.0 for count in partition_count]
    return percentages

# write the percentages in the bootstrap file
def write_bootstrap(percentages):
    with open('bootstrap.txt', 'w') as f:
        for percent in percentages:
            f.write(str(percent) + '\n')


def main(filename):
    # First step is to read the file and get the sequences
    ids, sequences = read_fasta_file(filename)

    # Then generate the distance matrix
    distMatrix = utils.get_distMatrix(ids, sequences)

    # write the distances.txt file (distance matrix)
    utils.write_distMatrix(ids, distMatrix)
  
    seqCounter = 120

    Neighbour_Joining_instance = neighbour_joining.Neighbour_Joining()
    root = Neighbour_Joining_instance.nei_saitou(ids, distMatrix, seqCounter)

    # Writing the edge file
    write_edge_file(root)

    # Writing the newick file
    write_newick_file(ids, root)

    #bootstrap calculations
    percentages = bootstrap(root, ids, sequences)
    write_bootstrap(percentages)


if __name__ == '__main__':
    main(sys.argv[1])
