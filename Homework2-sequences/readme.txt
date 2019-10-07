Required dependencies for the program are:
1. python 3
2. numpy
3. matplotlib

To test for the questions, use the following commands:

1) To run needleman wunsch, comment out the run_10000_perm line and uncomment the needleman_wunsch line in the file "homework2.py".
	Use python homework2.py -q Human_PAX.fa -r Fly_PAX.fa or python homework2.py -q Human_HOX.fa -r Fly_HOX.fa
2) To run the anchored version, 
	Use python homework2.py -q Human_PAX.fa -r Fly_PAX.fa -m Match_PAX.txt or python homework2.py -q Human_HOX.fa -r Fly_HOX.fa -m Match_HOX.txt
3) Same as the for 2)
4) To run needleman wunsch with 10000 permutation, comment the line for needleman_wunsch and uncomment the line run_10000_perm 
5)  Use python homework2.py -q TITIN_Human.fa -r TITIN_Mouse.fa -m TITIN_Match.txt
