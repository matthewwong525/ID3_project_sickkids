from __future__ import division
import math
from anytree import Node, RenderTree, LevelOrderIter
from anytree.exporter import DotExporter
from ga4gh_API import API

class ID3:

    def __init__(self):
        """
        Initializes the ID3 class

        Attributes:
        api (API): API object that is used to interact with the virtual API
        root_node (Node): Creates the root node of the tree to be added upon
        """
        self.api = API()
        subset = self.api.get_target_set()
        #print subset
        #print "\n\n\n\n"
        #print self.api.split_subset(split_var="16050654,A,[<CN0>, <CN2>, <CN3>, <CN4>]")[2]
        #print ["%s : %s" %j (sum(varz.values()), self.api.variant_name_list[idx]) for idx, varz in enumerate(self.api.split_subset(split_var="16050654,A,[<CN0>, <CN2>, <CN3>, <CN4>]")[2])]
        self.root_node = self.ID3(Node('root'), subset)

    @staticmethod
    def get_subset_count(subset):
        """
        Args:
            subset (dict): A dictionary containing keys of ancestries and values of the counts for the particular ancestry

        Returns:
            ret_str: A string containing counts of the subset and the name of the most common subset
        """
        max_key = str(max(subset, key=subset.get))
        total_count = str(sum(subset.values()))
        ret_str = "most common ancestry: %s,%s | total count: %s" % (max_key, subset[max_key], total_count)
        return ret_str

    @staticmethod
    def entropy_by_count(subset):
        """
        Gets the entropy given an input of a subset

        Args:
            subset (dict): A dictionary containing keys of ancestries and values of the counts for the particular ancestry

        Returns:
            entropy (float): the entropy value of the subset used for the ID3 algorithm
        """
        total_count = sum(subset.values())
        probability_list = []
        entropy = 0
        if total_count == 0:
            return 0
        for value in subset.values():
            probability = value / total_count
            if value != 0 : entropy -= (math.log(probability, 2) * probability)
        return entropy

    @staticmethod
    def create_split_path(split_path, new_variant_name):
        """
        Creates two new split paths given the variant name and the split path

        Args:
            subset (dict): A dictionary containing keys of ancestries and values of the counts for the particular ancestry

        Returns:
            w_split_path: The split path with the variant
            wo_split_path: The split path without the variant

        """
        w_split_path = (list(split_path[0]), list(split_path[1]))
        wo_split_path = (list(split_path[0]), list(split_path[1]))
        w_split_path[0].append(new_variant_name)
        wo_split_path[0].append(new_variant_name)
        w_split_path[1].append(1)
        wo_split_path[1].append(0)

        return w_split_path, wo_split_path


    def is_leaf_node(self, subset, split_path, split_index):
        """
        Checks if the node is a leaf node given a subset

        Args:
            subset (dict): A dictionary containing keys of ancestries and values of the counts for the particular ancestry
            split_paths (list1, list2): 
                This is the paths of the splits before the current split. The first list
                is the list of variant names and the second list is the direction
                of the split. The direction of the second list is depicted by 1's
                and 0's. Where 1 is splitting in the direction with the variant
                and 0 is splitting in the direction without the variant. 
            split_index (int): the index for the next split to be split on
        Returns:
            Boolean Value: the boolean value represents whether or not the node is a leaf node

        """
        # check if all variants are of one ancestry (Essentially if all remaining variants(attributes) contains one region(value))
        attr_list = set()
        upd_var_count = { k : v for k, v in subset.items() if v != 0 }
        attr_list.update(upd_var_count.keys())
        #if len(attr_list) == 1 or len(split_path[0]) >= len(self.api.variant_name_list) or split_index is None or ID3.entropy_by_count(subset) == 0:
        if len(attr_list) == 1 or len(split_path[0]) >= 5 or split_index is None or ID3.entropy_by_count(subset) == 0:
            return True
        return False

    def print_tree(self, file_name):
        DotExporter(self.root_node).to_picture("%s.png" % file_name)

    def find_variant_split(self, subset, split_path):
        """
        Finds the variant to split on and returns the index where it should be split on.
        This calculation is based on which attribute gives the greatest information gain.

        Args:
            subset (dict): A dictionary containing keys of ancestries and values of the counts for the particular ancestry
            split_paths (list1, list2): 
                This is the paths of the splits before the current split. The first list
                is the list of variant names and the second list is the direction
                of the split. The direction of the second list is depicted by 1's
                and 0's. Where 1 is splitting in the direction with the variant
                and 0 is splitting in the direction without the variant. 

        Returns:
            ret_index: index that yields the greatest information gain
        """
        variant_list = self.api.split_subset(split_path)[2]
        total_count = sum(subset.values())
        ret_index = 0
        final_info_gain = 0
        var_idx_list = [self.api.variant_name_list.index(var_name) for var_name in split_path[0]]

        if total_count == 0:
            return None

        # loops through all the counts for each variant
        for idx, w_var_counts in enumerate(variant_list):
            wo_var_counts = { k : subset.get(k, 0) - w_var_counts.get(k, 0) for k in subset.keys() + w_var_counts.keys() }

            # calculates info gain
            info_gain = ID3.entropy_by_count(subset) - ( sum(wo_var_counts.values()) / total_count * ID3.entropy_by_count(wo_var_counts) + sum(w_var_counts.values()) / total_count * ID3.entropy_by_count(w_var_counts) )
            # finds max info gain and checks if the index is not in excluded variants list
            if final_info_gain < info_gain and idx not in var_idx_list:
                final_info_gain = info_gain
                ret_index = idx

        return ret_index

    # note: variant list must be same length as count list
    def ID3(self, node, subset, split_path=([], [])):
        """
        A recursive function that creates a tree given the root node and a subset

        Note: variant list must be same length as count list

        Args:
            node (Node): A node object from the anytree library
            subset (dict): A dictionary containing keys of ancestries and values of the counts for the particular ancestry
            split_paths (list1, list2): 
                This is the paths of the splits before the current split. The first list
                is the list of variant names and the second list is the direction
                of the split. The direction of the second list is depicted by 1's
                and 0's. Where 1 is splitting in the direction with the variant
                and 0 is splitting in the direction without the variant. 

        TODO:
            * Delete the extra nodes in the tree where the tree doesn't split

        """
        # find the attrivute to split on and adds that variant to exclude variant list
        split_index = self.find_variant_split(subset, split_path)
        if not self.is_leaf_node(subset, split_path, split_index):
            var_name = self.api.variant_name_list[split_index]
            
            w_subset, wo_subset = self.api.split_subset(split_path, var_name)[:2]

            w_split_path, wo_split_path = ID3.create_split_path(split_path, var_name)


            if sum(w_subset.values()) > 0:
                self.ID3(Node("with:%s\n%s" % (var_name, ID3.get_subset_count(w_subset)), parent=node), dict(w_subset), w_split_path)
            if sum(wo_subset.values()) > 0:
                self.ID3(Node("without:%s\n%s" % (var_name, ID3.get_subset_count(wo_subset)), parent=node), dict(wo_subset), wo_split_path)
            return node

if __name__ == "__main__":
    id3_alg = ID3()
    id3_alg.print_tree('udo1')
