To run the code:
python main.py data/hw3.fna

This generates the files by default:
1. genetic_distances.txt
2. edges.txt
3. tree.txt
4. bootstrap.txt

Running the bootstrap code take a little longer. To avoid that,everything that
follow the comment #bootstrap calculations# in main.py can be commented out.

For visualization:

Rscript hw3-plot-newick.r tree.txt hw3-tip-labels.txt
or
Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt

For bonus visualization:

Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt boot.txt
