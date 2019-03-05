import json
import requests

class GA4GH_API:
    host_url = 'http://localhost:8000/'
    def __init__(self, file_path, conf_matrix=False):
        """
        Initializes the GA4GH_API class

        NOTE:
            api_variant_list, api_indiv_list, and api_popu_list should be of the same
            length and every element in this list corresponds to a particular person

        Args:
            file_path (str): Path to json file that contains the variant ranges
            conf_matrix (bool): Initializes API to perform conf_matrix operations

        Attributes:
            variant_name_list (list): Names of the variants in the format of
                                          "VARIANT_POS,VARIANT_REF,VARIANT_ALT"
                                          "CHROMOSOME_#:START_POS:END_POS" (TODO - UPDATE TO THIS)
            ancestry_list (list): A unique list of all the ancestries of people
            is_conf_matrix (bool): tells API to initialize the API for confusion matrix operations

        TODO:
            * Add option to query for different variants rather than grabbing a predefined set of variants
        """
        self.dataset_id = 'WyIxa2dlbm9tZSJd'
        self.variant_set_ids = GA4GH_API.get_variant_set_ids(self.dataset_id)[::10]
        self.variant_name_list = self.fetch_variants(file_path)
        self.is_conf_matrix = conf_matrix
        self.ancestry_list = []

        # updates variables
        #self.read_user_mappings(variant_dict)

    @staticmethod
    def get_variant_set_ids(dataset_id):
        req_body = { "dataset_id": dataset_id }
        r = requests.post('%s%s' %  (GA4GH_API.host_url, 'variantsets/search'), json=req_body).json()
        variant_set_ids = [ variantset['id'] for variantset in r['results']['variantSets'] ]
        return variant_set_ids

    def fetch_variants(self, file_path):
        variant_list = []
        with open(file_path) as f:
            data = json.load(f)
        for var_range in data:
            # TODO: queries for the variants within the range here:
            chrom = str(var_range['chr']).replace('chr', '')
            variant_list.extend(self.query_variants(chrom, str(var_range['start']), str(var_range['end'])))
        return variant_list

    def query_variants(self, chrom, start, end):
        """
        Queries the POS value of the individual variants within a range of variants.
        Returns the variable variant_name_list to the variants found within
        this range of variants.

        Args:
            chrom (str): chromosome number
            start (str): starting position of variant
            end (str): ending position of variant

        Returns: 
            variant_name_list (list): list of variant names formatted in the for of `CHR:POS`
        """
        variant_list = []
        req_body = {
            'variantSetIds' : self.variant_set_ids[:1],
            'start': start,
            'end': end,
            'referenceName': chrom
        }
        r = requests.post('%s%s' %  (GA4GH_API.host_url, 'variants/search'), json=req_body).json()
        for variant in r['results']['variants']:
            variant_list.append(':'.join([chrom, variant['start'], variant['end']]))
        return variant_list


    def craft_api_request(self, split_paths=([], [])):
        """
        Crafts an GA4GH_API body to be sent to the ga4gh server which is intended to filter
        the counts of the variants based on the inclusion or exclusion of particular
        variants

        Attributes:
            split_paths (list1, list2): 
                    This is the paths of the splits before the current split. The first list
                    is the list of variant names and the second list is the direction
                    of the split. The direction of the second list is depicted by 1's
                    and 0's. Where 1 is splitting in the direction with the variant
                    and 0 is splitting in the direction without the variant.
            split_var (string): The variant name it is now splitting on

        Returns:
            req_body (json): Returns a JSON containing the request body
        TODO:
            * Add NAND logic for filtering without a particular variant
        """
        components = []
        logic = { 'and': 
            [ 
                { 'or': [] },
                { 'and': [] }
            ] 
        }
        id_list = []
        # Add to OR list
        for variant_id in self.variant_name_list:
            id_list.append( { 'id': variant_id } )
            CHR, START, END = variant_id.split(':')
            components.append(
                {
                    "id": variant_id,
                    "variants":{
                        "start": START,
                        "end": END,
                        "referenceName": CHR,
                        "variant_set_ids": self.variant_set_ids
                    }
                })
        logic['and'][0]['or'].extend(id_list)

        for variant, direction in zip(split_paths[0], split_paths[1]):
            # remove from OR
            logic['and'][0]['or'] = [ variant_id for variant_id in logic['and'][0]['or'] if variant != variant_id['id']]
            if not direction:
                continue
            # add to AND
            logic['and'][1]['and'].append({"id": variant})

        # Finds index with empty list and removes it
        i = None
        for idx, val in enumerate(logic['and']):
            key = val.keys()[0]
            if val[key] == []:
                i = idx
                break
        if i: 
            del logic['and'][i]

        req_body = {}
        req_body['logic'] = logic
        req_body['components'] = components
        req_body['dataset_id'] = self.dataset_id
        req_body['results'] = [ {
                    "table": "patients",
                    "field": [
                        "ethnicity"
                    ]
                } ]
        req_body['page_size'] = 10000
        return req_body          

    def get_target_set(self):
        """
        Gets the target subset, which is the ancestry counts every variant

        Returns:
            counts (dict): A dictionary containing keys of ancestries and values of the counts for the particular ancestry
        """
        req =  self.craft_api_request()
        ancestry_counts = requests.post('%s%s' %  (GA4GH_API.host_url, 'count'), json=req).json()['results']['patients'][0]['ethnicity']
        if self.ancestry_list == []:
            self.ancestry_list = ancestry_counts.keys()

        return ancestry_counts

    def split_subset(self, node, split_var):
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

        """
        wo_variant_split_path = list(node.split_path)
        wo_variant_split_path.append((split_var, 0))
        w_variant_split_path = list(node.split_path)
        w_variant_split_path.append((split_var, 1))

        w_var_req_body = self.craft_api_request(w_variant_split_path)
        wo_var_req_body = self.craft_api_request(wo_variant_split_path)

        # make query here for
        r_w_var = requests.post('%s%s' %  (GA4GH_API.host_url, 'count'), json=w_var_req_body).json()['results']['patients'][0]['ethnicity']
        r_wo_var = requests.post('%s%s' %  (GA4GH_API.host_url, 'count'), json=w_var_req_body).json()['results']['patients'][0]['ethnicity']

        print("performed split")
        print(r_w_var)
        print(w_var_req_body)

        print(r_wo_var)
        print(wo_var_req_body)
        print("")

        return r_w_var, r_wo_var

        # For w_variant list

    def find_next_variant_counts(self, split_paths):
        """
        Finds the counts of the a potential next variant to perform the
        split on

        Attributes:
            split_paths (list1, list2): 
                This is the paths of the splits before the current split. The first list
                is the list of variant names and the second list is the direction
                of the split. The direction of the second list is depicted by 1's
                and 0's. Where 1 is splitting in the direction with the variant
                and 0 is splitting in the direction without the variant.
        Returns:
            w_variant_list: 
                A list representing of a dictionary of ancestry counts per variant
                [
                    {'GBR': 5, 'CML': 0, 'ABC': 3 ...}
                    ...
                    ..
                    .
                ]
        """
        w_variant_list = []
        exclude_variants = [ variant for variant in split_paths[0] ]
        include_variants = list(self.variant_name_list)

        for var in exclude_variants:
            include_variants.remove(var)

        for var in include_variants:
            split_paths[0].append(var)
            split_paths[1].append(1)

            req_body = self.craft_api_request(split_paths)
            resp = requests.post('%s%s' %  (GA4GH_API.host_url, 'count'), json=req_body).json()
            variant_counts = resp['results']['patients'][0]['ethnicity']

            w_variant_list.append(variant_counts)

            del split_paths[0][-1]
            del split_paths[1][-1]

        return w_variant_list



