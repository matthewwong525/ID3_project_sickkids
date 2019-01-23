import vcf
import itertools
from sklearn import tree

import codecs
from ped_parser import parser

vcf_reader = vcf.Reader(open('chr22.vcf', 'r'))
variant_dict = {}
indiv_list = []
test_ans = []

## fetching 20 records
# TODO: fetch 20 records that are randomly apart not in sequence
for count in range(0,1):
    init_pos = 16050074
    #16056319
    #variants = vcf_reader.fetch('22', init_pos + (600*count), init_pos + 1 + (100000*count))
    variants = vcf_reader.fetch('22', 16050074, 16050075+150000)

    # loops through the variants
    for variant in variants:
        genotype_list = []

        # creates dictionary where:
            #   key is individual id
            #   value is list of binary values where (1 = variant exists, 0 = variant doesn't exist)
        for call in variant.samples:
            if call.sample in variant_dict:
                if call['GT'] == '0|0':
                    variant_dict[call.sample].append(0)
                else:
                    variant_dict[call.sample].append(1)
            else:
                variant_dict[call.sample] = []

## parsing ped file for ancestry classification
## also produces lists/matrices for the tree algorithm

# list of population used as the  
popu_list = [] 
variant_list = []
indiv_list = []

# loops through every line in the ped file
with open('user_map.ped') as file:
    next(file)
    for line in file:
        split_line = line.split('\t')
        indiv_id = split_line[1]
        # checks if the id in ped map exists in vcf file
        if indiv_id in variant_dict:
            population = split_line[6]
            indiv_list.append(indiv_id)
            popu_list.append(population)
            variant_list.append(variant_dict[indiv_id])

## creates tree based on the lists created
clf = tree.DecisionTreeClassifier()
clf = clf.fit(variant_list, popu_list)

print list(set(popu_list))
print clf.predict_proba([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]])
print clf.predict([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]])
