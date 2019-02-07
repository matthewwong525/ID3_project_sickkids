from __future__ import division

import math
import vcf
import itertools
from sklearn import tree
from anytree import Node, RenderTree
from anytree.exporter import DotExporter

import codecs
from ped_parser import parser
from collections import Counter

api_variant_list = []
api_variant_name_list = []
api_indiv_list = []
api_popu_list = []
ancestry_dict = {}

def main():
    global api_indiv_list
    global api_popu_list
    global api_variant_name_list
    global api_variant_list
    global ancestry_dict

    vcf_reader = vcf.Reader(open('chr22.vcf', 'r'))
    variant_dict = {}
    test_ans = []

    # 
    variants = vcf_reader.fetch('22', 16050074, 16050075+150000)
    variant_dict, api_variant_name_list = create_variant_dict(variants)

    # creates a dictionary of ancestries based on metadata
    # key: individual id
    # value: ancestry/population
    with open('user_map.ped') as file:
        next(file)
        for line in file:
            split_line = line.split('\t')
            indiv_id = split_line[1]
            population = split_line[6]
            ancestry_dict[indiv_id] = population
            if indiv_id in variant_dict:
                api_indiv_list.append(indiv_id)
                api_popu_list.append(population)
                api_variant_list.append(variant_dict[indiv_id])

    # API call here to get original subset (split on nothing)
    #subset = split_subset([api_variant_name_list[0]])[0]
    subset =  split_target()
    #print api_variant_name_list
    #print(split_subset(split_var='16050922,T,[G]')[2][15])
    #find_variant_split(subset,[])

    udo = ID3(Node('root'), subset)
    #, split_paths=(['16050075,A,[G]'], [0])
    #print(split_subset(split_var='16050319,C,[T]', split_paths=(['16050075,A,[G]'], [0]))[:2])

    #variant_count_list = [count_by_variant(variant['variant'], ancestry_dict)[1] for variant in variant_list]
    #udo = ID3(variant_list, Node('root\n%s' % str(get_summed_counts(variant_count_list))), ancestry_dict)
    DotExporter(udo).to_picture("udo.png")
    #print(RenderTree(udo))

def split_target():
    ancestry_list = set(ancestry_dict.values())
    #ancestry_list.remove('CHD')
    counts = dict.fromkeys(ancestry_list, 0)
    for i in api_popu_list:
        counts[i] = counts.get(i, 0) + 1
    return counts

def create_variant_dict(variants):
    variant_dict = {}
    variant_name_list = []
    # loops through the variants
    idx = 0
    for variant in variants:
        idx += 1
        variant_name_list.append(','.join([str(variant.POS), str(variant.REF), str(variant.ALT)]))
        # creates dictionary where:
            #   key is individual id
            #   value is list of binary values where (1 = variant exists, 0 = variant doesn't exist)
        for call in variant.samples:
            variant_dict[call.sample] = variant_dict.get(call.sample, [])
            if call['GT'] == '0|0':
                variant_dict[call.sample].append(0)
            else:
                variant_dict[call.sample].append(1)
    return variant_dict, variant_name_list


# input: the variant and ancestry dict
# output: the variant count w/ and w/o variants
def count_by_variant(variant, ancestry_dict):

    # finds ancestries
    ancestry_list = set(ancestry_dict.values())
    #ancestry_list.remove('CHD')

    # creates a list of counts based on ancestry
    count_w_var = dict.fromkeys(ancestry_list, 0)

    for call in variant.samples:
        if call.sample in ancestry_dict and '1' in call['GT']:
            ancestry = ancestry_dict[call.sample]
            count_w_var[ancestry]  += 1
    return count_w_var

# calculates the entropy given one set of variant counts
def entropy_by_count(count):
    total_count = sum(count.values())
    probability_list = []
    entropy = 0
    if total_count == 0:
        return 0
    for value in count.values():
        probability = value / total_count
        if value != 0 : entropy -= (math.log(probability, 2) * probability)
    return entropy

# input: variant_count_list
# outpput: index on list to split on
def find_variant_split(subset, split_path):
    # sum appearance of variant for each ancestry
    variant_list = split_subset(split_path)[2]
    total_count = sum(subset.values())
    ret_index = 0
    final_info_gain = 0
    var_idx_list = [api_variant_name_list.index(var_name) for var_name in split_path[0]]

    if total_count == 0:
        return None

    # loops through all the counts for each variant
    for idx, w_var_counts in enumerate(variant_list):
        wo_var_counts = { k : subset.get(k, 0) - w_var_counts.get(k, 0) for k in subset.keys() + w_var_counts.keys() }
        info_gain = entropy_by_count(subset)
        for popu in w_var_counts.keys():
            info_gain -= wo_var_counts[popu] / total_count * entropy_by_count(w_var_counts)
        #print info_gain
        if final_info_gain < info_gain and idx not in var_idx_list:
            final_info_gain = info_gain
            ret_index = idx
    #print final_info_gain
    return ret_index

    '''
    index = 0
    entropy = entropy_by_count(variant_count_list[0])
    for cur_idx, variant_count in enumerate(variant_count_list):
        curr_entropy = entropy_by_count(variant_count)
        if entropy == 0 : return index
        if curr_entropy < entropy:
            entropy = curr_entropy
            index = cur_idx
    return index
    '''

def get_summed_counts(variant_count_list):
    summed_dict = {}
    max_ancestry = ''
    max_ancestry_count = 0
    summed_dict['Remaining Count Sum'] = 0

    for variant_count in variant_count_list:
        for ancestry in variant_count.keys():
            summed_dict['Remaining Count Sum'] += variant_count[ancestry]

            if max_ancestry_count < variant_count[ancestry]:
                max_ancestry = ancestry
                max_ancestry_count = variant_count[ancestry]


    summed_dict['Remaining Count Sum'] -= max_ancestry_count
    summed_dict[max_ancestry] = max_ancestry_count
    return summed_dict

def get_subset_count(subset):
    max_key = str(max(subset, key=subset.get))
    total_count = str(sum(subset.values()))
    ret_str = "ancestry: %s,%s | total count: %s" % (max_key, subset[max_key], total_count)
    return ret_str

# splits the set given a variant 
# returns 2 subsets of the data
def split_subset(split_paths=([], []), split_var=None):
    # retrieves variant from "API"
    ancestry_list = set(ancestry_dict.values())
    #ancestry_list.remove('CHD')
    w_variant_list = [dict.fromkeys(ancestry_list, 0) for variant_names in api_variant_name_list]
    wo_variant_list = [dict.fromkeys(ancestry_list, 0) for variant_names in api_variant_name_list]
    w_variant_dict = dict.fromkeys(ancestry_list, 0)
    wo_variant_dict = dict.fromkeys(ancestry_list, 0)

    ignore_rows_idxs = []

    # basically find rows to ignore because it has been split upon already
    for exc_var, direction in zip(split_paths[0], split_paths[1]):
        # finding rows to ignore
        var_idx = api_variant_name_list.index(exc_var)
        # looping through all the people
        for idx, variants in enumerate(api_variant_list):
            if idx not in ignore_rows_idxs and variants[var_idx] is not direction:
                ignore_rows_idxs.append(idx)
    #print ignore_rows_idxs

    # create new subset after finding all the rows to ignore
    for idx, variants in enumerate(api_variant_list):
        if idx not in ignore_rows_idxs:
            popu = api_popu_list[idx]
            # check if split_var is null
            if split_var:
                f_var_idx = api_variant_name_list.index(split_var)
                if variants[f_var_idx] is 1:
                    w_variant_dict[popu] += 1
                else:
                    wo_variant_dict[popu] += 1
            # find counts of variants
            for idx2, variant in enumerate(variants):
                if variant is 1: 
                    w_variant_list[idx2][popu] += 1


    return w_variant_dict, wo_variant_dict, w_variant_list, ignore_rows_idxs




def create_split_path(split_path, new_variant_name):
    w_split_path = (list(split_path[0]), list(split_path[1]))
    wo_split_path = (list(split_path[0]), list(split_path[1]))
    w_split_path[0].append(new_variant_name)
    wo_split_path[0].append(new_variant_name)
    w_split_path[1].append(1)
    wo_split_path[1].append(0)

    return w_split_path, wo_split_path
# note: variant list must be same length as count list
def ID3(parent_node, subset, split_path=([], [])):
    # find the attrivute to split on and adds that variant to exclude variant list
    #print get_subset_count(subset)
    split_index = find_variant_split(subset, split_path)

    #len(api_variant_name_list)
    if len(split_path[0]) < len(api_variant_name_list) and not is_leaf_node(subset) and split_index is not None:
        #print split_index
        var_name = api_variant_name_list[split_index]
        

        print split_path, var_name
        w_subset, wo_subset = split_subset(split_path, var_name)[:2]
        #print split_subset(split_path, var_name)[2]

        #print w_subset, wo_subset
        #print

        w_split_path, wo_split_path = create_split_path(split_path, var_name)

        if sum(w_subset.values()) > 0:
            ID3(Node("with:%s\n%s" % (var_name, get_subset_count(w_subset)), parent=parent_node), dict(w_subset), w_split_path)
        if sum(wo_subset.values()) > 0:
            ID3(Node("without:%s\n%s" % (var_name, get_subset_count(wo_subset)), parent=parent_node), dict(wo_subset), wo_split_path)
        return parent_node


    '''# checks if variant list is empty
    if variant_list != []:
        # loops through remaining variants and creates a tree
        pop_index = find_variant_split(variant_count_list)
        var_name = variant_list[pop_index]['name']
        # add variant to exclude to a list
        exclude_variant_list.append(variant_list.pop(pop_index)['variant'])

        # makes the api call here gets lists with and without variant
        with_variant_count_list = [count_by_variant(variant['variant'], ancestry_dict, exclude_variant)[0] for variant in variant_list]
        without_variant_count_list = [count_by_variant(variant['variant'], ancestry_dict, exclude_variant)[1] for variant in variant_list]
        with_sum = str(get_summed_counts(with_variant_count_list))
        without_sum = str(get_summed_counts(without_variant_count_list))

        print with_sum
        # TODO: make it print out the answer
        ID3(variant_list, Node("with:%s\n%s" % (var_name, with_sum), parent=parent_node), ancestry_dict)
        ID3(variant_list, Node("without:%s\n%s" % (var_name, without_sum), parent=parent_node), ancestry_dict)'''

        


# checks if it is a leaf node
def is_leaf_node(subset):
    # check if all variants are of one ancestry (Essentially if all remaining variants(attributes) contains one region(value))
    attr_list = set()
    # deletes all 0 values in dictionary
    upd_var_count = { k : v for k, v in subset.items() if v != 0 }
    attr_list.update(upd_var_count.keys())
    if len(attr_list) == 1:
        return True
    return False


main()
