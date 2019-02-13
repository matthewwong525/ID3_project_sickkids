import vcf

class API:
    def __init__(self):
        """
        Initializes the API class

        NOTE:
            api_variant_list, api_indiv_list, and api_popu_list should be of the same
            length and every element in this list corresponds to a particular person

        Attributes:
            api_variant_list (list): Represents the variants in each person. 
                                     The value of the variants is a list of 0's and 
                                     1's indicating if the variant exists in the person.
            api_indiv_list (list): Represents the individual code of a person
            api_popu_list (list): Represents the ancestry of the person
            api_variant_name_list (list): All the variant names queried
            api_variant_name_list (list): Names of the variants in the format of
                                          "VARIANT_POS,VARIANT_REF,VARIANT_ALT"
            ancestry_dict (dict): a dictionary which maps indivdual ID to population

        TODO:
            * Add option to query for different variants rather than grabbing a predefined set of variants
        """
        self.variant_list = []
        self.indiv_list = []
        self.popu_list = []
        self.variant_name_list = []
        self.ancestry_dict = {}

        # fetch variants from vcf and create a dictionary
        vcf_reader = vcf.Reader(open('chr22.vcf', 'r'))
        variants = vcf_reader.fetch('22', 16050074, 16050075+150000)
        variant_dict = self.create_variant_dict(variants)

        # updates variables
        self.read_user_mappings(variant_dict)

    def create_variant_dict(self, variants):
        """
        Creates a ditionary of variants from a vcf file which is used to update variables
        within this class

        Args:
            variants (Reader): A Reader Object from the PyVCF library which contains variants

        Returns:
            variant_dict (dict): A dictionary where the
                                    key: is the individual ID
                                    value: is a list of binary values where
                                           (1 = variant exists, 0 = variant doesn't exist)
        """
        variant_dict = {}
        variant_name_list = []
        idx = 0
        # loops through variants
        for variant in variants:
            idx += 1
            self.variant_name_list.append(','.join([str(variant.POS), str(variant.REF), str(variant.ALT)]))
            # loops through people in variants
            for call in variant.samples:
                variant_dict[call.sample] = variant_dict.get(call.sample, [])
                # checks if variant exists in person
                if call['GT'] == '0|0':
                    variant_dict[call.sample].append(0)
                else:
                    variant_dict[call.sample].append(1)
        return variant_dict

    def read_user_mappings(self, variant_dict):
        """
        Reads the usermappings from a file and updates the variables in the class
        """
        with open('user_map.ped') as file:
            next(file)
            for line in file:
                split_line = line.split('\t')
                indiv_id = split_line[1]
                population = split_line[6]

                self.ancestry_dict[indiv_id] = population

                # checks if variant individual is in the variant dict from the vcf file
                if indiv_id in variant_dict:
                    self.indiv_list.append(indiv_id)
                    self.popu_list.append(population)
                    self.variant_list.append(variant_dict[indiv_id])


    # splits the set given a variant 
    # returns 2 subsets of the data
    def split_subset(self, split_paths=([], []), split_var=None):
        """
        Splits the subset given the path of the splits before and a
        variable to split on.

        Attributes:
            split_paths (list1, list2): 
                This is the paths of the splits before the current split. The first list
                is the list of variant names and the second list is the direction
                of the split. The direction of the second list is depicted by 1's
                and 0's. Where 1 is splitting in the direction with the variant
                and 0 is splitting in the direction without the variant.
            split_var (string): The variant name it is now splitting on

        Returns:
            w_variant_dict (dict): The split subset that includes the variant
            wo_variant_dict (dict): The split subset that does not include the variant
            w_variant_list: 
                A list representing of a dictionary of ancestry counts per variant
                [
                    {'GBR': 5, 'CML': 0, 'ABC': 3 ...}
                    ...
                    ..
                    .
                ]

        """
        # retrieves variant from "API"
        ancestry_list = set(self.ancestry_dict.values())
        w_variant_list = [dict.fromkeys(ancestry_list, 0) for variant_names in self.variant_name_list]
        wo_variant_list = [dict.fromkeys(ancestry_list, 0) for variant_names in self.variant_name_list]
        w_variant_dict = dict.fromkeys(ancestry_list, 0)
        wo_variant_dict = dict.fromkeys(ancestry_list, 0)

        ignore_rows_idxs = []

        # basically find rows to ignore because it has been split upon already
        for exc_var, direction in zip(split_paths[0], split_paths[1]):
            # finding rows to ignore
            var_idx = self.variant_name_list.index(exc_var)
            # looping through all the people
            for idx, variants in enumerate(self.variant_list):
                if idx not in ignore_rows_idxs and variants[var_idx] is not direction:
                    ignore_rows_idxs.append(idx)

        # create new subset after finding all the rows to ignore
        for idx, variants in enumerate(self.variant_list):
            if idx not in ignore_rows_idxs:
                popu = self.popu_list[idx]
                # check if split_var is null
                if split_var:
                    f_var_idx = self.variant_name_list.index(split_var)
                    if variants[f_var_idx] is 1:
                        w_variant_dict[popu] += 1
                    else:
                        wo_variant_dict[popu] += 1
                # find counts of variants
                for idx2, variant in enumerate(variants):
                    if variant is 1: 
                        w_variant_list[idx2][popu] += 1


        return w_variant_dict, wo_variant_dict, w_variant_list

    def get_target_set(self):
        """
        Gets the target subset

        Returns:
            counts (dict): A dictionary containing keys of ancestries and values of the counts for the particular ancestry
        """
        ancestry_list = set(self.ancestry_dict.values())
        counts = dict.fromkeys(ancestry_list, 0)
        for i in self.popu_list:
            counts[i] = counts.get(i, 0) + 1
        return counts

    def count_variants(self):
        my_dict = {}
        for variants in (self.variant_list):
            count = 0
            for idx, variant in enumerate(variants):
                if variant is 1:
                    my_dict[self.variant_name_list[idx]] = my_dict.get(self.variant_name_list[idx], 0) + 1
        return my_dict

if __name__ == "__main__":
    api = API()
    print api.variant_name_list



