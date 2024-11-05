import sys, os
# add path to root folder to path (change this based on where you're running this script from)
sys.path.append("../../")
# from (name of folder) import (name of binary)
import lib.fold.fold as fold

test_matrix = fold.BasePairMatrix("test.tsv")
test_matrix.toCSV("test_as_matrix.csv")
pairs = []
test_matrix.getBestPairing(pairs)
with open("test_output.txt", 'w') as outfile:
    outfile.write("i coord\tj coord\ti nuc\tj nuc\tz-norm\n")
    for pair in pairs:
        line = str(pair.i_coord) + "\t" + str(pair.j_coord)  + "\t" + str(pair.i_nucleotide) + "\t" + str(pair.j_nucleotide) + "\t" + str(pair.getZNorm()) + "\n"
        outfile.write(line)