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
            if row[2] >= 97:
                total_similarity = total_similarity + float(row[2])
                percentage_match+=1
            if(row[12] in bacterial_dict):
                bacterial_dict[row[12]]+=1
            else:
                bacterial_dict[row[12]] = 1

        fraction_match = percentage_match/total_count
        key,value = sorted(bacterial_dict.items(), key=lambda x:x[1], reverse=True)[0]
        average_similarity = total_similarity / percentage_match

                
                        
        print(total_count)
        print(fraction_match)
        print(key)
        print(value)
        print(average_similarity)
if __name__ == '__main__':
    main()