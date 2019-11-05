from node import Node
import sys
import random
from neighbour_joining import Neighbour_Joining
import utils


# uses preorder traversal to traverse the tree and generate the edges file
def write_edge_file(root):
    traversal = []
    # recursive preorder traversal
    def preorder_traversal(node):
        if not node or not node.children:
            return
        for child, distance in node.children.items():
            print(node.id)
            print(child.id)
            # traverse root first
            traversal.append((node.id, child.id, distance))
            # then traverse children recursively
            preorder_traversal(child)
    preorder_traversal(root)
    # writes the traversal list to file
    with open('edges.txt', 'w') as f:
        for (parent, child, distance) in traversal:
            f.write(parent + '\t' + child + '\t' + str(distance) + '\n')


# uses post order traversal to write the newick file
def write_newick_file(original_ids, root):
    # recursive postorder traversal
    def postorder_traversal(node):
        if not node.children:
            # if it a leaf node then return the original id
            if 1 <= int(node.id) <= 61:
                return original_ids[int(node.id) - 1]
        vals = []
        # recursively traverse the children first
        for child, distance in node.children.items():
            vals.append(postorder_traversal(child) + ':' + str(distance))
        result = '(' + ','.join(vals) + ')'
        return result
    # get the root nodes children.
    # split it up into a 2 pair and 1 extra node for the root
    (first, fd), (second, sd), (third, td) = root.children.items()
    # recursively call the postorder traversal function to generate the newick format
    first_part = '(' + postorder_traversal(second) + ':' + str(sd) + ','+ postorder_traversal(third) + ':' + str(td) + ')'
    newick =  '(' + first_part + ':' + str(fd) + ');'
    with open('tree.txt', 'w') as f:
        f.write(newick)


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
        distance_matrix = get_difference_matrix(ids, bootstrap_sequences)
        root = nei_saitou(ids, distance_matrix)
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
    distance_matrix = utils.get_distance_matrix(ids, sequences)

    # write the distances.txt file (distance matrix)
    utils.write_distance_matrix(ids, distance_matrix)
  
    root = Neighbour_Joining.nei_saitou(ids, distance_matrix, 120)

    # write the edge file using preorder traversal
    write_edge_file(root)

    #write the newick file using postorder traversal
    write_newick_file(ids, root)

    # bootstrap bonus points
    #percentages = bootstrap(root, ids, id_sequences)
    #write_bootstrap(percentages)


if __name__ == '__main__':
    main(sys.argv[1])
