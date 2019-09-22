#!/usr/bin/env python
import sys
import csv

def main():
    with open(r'output.txt', 'r') as f:
        total_count = 0
        percentage_match = 0;
        bacterial_dict = {}
        spamreader = csv.reader(f, delimiter='\t')
        total_similarity = 0
        for row in spamreader:
            total_count += 1
            if (float(row[2]) >= 97):
                total_similarity = total_similarity + float(row[2])
                percentage_match+=1
            name = row[12]
            if(name in bacterial_dict):
                bacterial_dict[name]+=1
            else:
                bacterial_dict[name] = 1

        fraction_match = float(percentage_match)/total_count
        key,value = sorted(bacterial_dict.items(), key=lambda x:x[1], reverse=True)[0]
        average_similarity = total_similarity / percentage_match 
                        
        print "Percentage match >= 97:" , fraction_match
        print "Most common species in the query is ", key.split("s__")[1], "with", value, "occurrences."
        print "Average Similarity:" , average_similarity
if __name__ == '__main__':
    main()
