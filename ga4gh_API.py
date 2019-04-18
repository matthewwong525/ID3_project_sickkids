import json
import requests

class GA4GH_API:
    def __init__(self, file_path):
        """
        Initializes the GA4GH_API class

        NOTE:
            api_variant_list, api_indiv_list, and api_popu_list should be of the same
            length and every element in this list corresponds to a particular person

        Args:
            file_path (str): Path to json file that contains the variant ranges

        Attributes:
            config (json): loaded config file
            host_url (str): url pointing to ga4gh_server
            dataset_id (str): dataset id of the ga4gh_server that is to b accessed
            variant_name_list (list): Names of the variants in the format of
                                          "VARIANT_POS,VARIANT_REF,VARIANT_ALT"
                                          "CHROMOSOME_#:START_POS:END_POS" (TODO - UPDATE TO THIS)
            ancestry_list (list): A unique list of all the ancestries of people
            is_conf_matrix (bool): tells API to initialize the API for confusion matrix operations

        TODO:
            * Throw error when server gives incorrect response
        """
        with open(file_path) as f:
            self.config = json.load(f)
        self.host_url = self.config['ga4gh_server_url']
        self.dataset_id = self.config['ga4gh_server_dataset_id']
        self.variant_name_list = self.fetch_variants(file_path)
        self.ancestry_list = []

        # updates variables
        #self.read_user_mappings(variant_dict)

    @staticmethod
    def create_split_path(split_path, new_variant_name):
        """
        Creates two new split paths given the variant name and the split path

        Args:
            split_path (list1, list2): 
                This is the paths of the splits before the current split. The first list
                is the list of variant names and the second list is the direction
                of the split. The direction of the second list is depicted by 1's
                and 0's. Where 1 is splitting in the direction with the variant
                and 0 is splitting in the direction without the variant. 
            new_variant_name (str): unique id of variant

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

    def fetch_variants(self, file_path):
        variant_list = []
        for var_range in self.config['variant_ranges']:
            variant_list.extend(self.query_variants(str(var_range['chr']), str(var_range['start']), str(var_range['end'])))
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
            'datasetId' : self.dataset_id,
            'start': start,
            'end': end,
            'referenceName': chrom
        }
        r = requests.post('%s%s' %  (self.host_url, 'variants/search'), json=req_body).json()
        for variant in r['results']['variants']:
            variant_list.append(':'.join([chrom, variant['start'], variant['end']]))
        return variant_list


    def craft_api_request(self, split_path=([], [])):
        """
        Crafts an GA4GH_API body to be sent to the ga4gh server which is intended to filter
        the counts of the variants based on the inclusion or exclusion of particular
        variants

        Attributes:
            split_path (list1, list2): 
                    This is the paths of the splits before the current split. The first list
                    is the list of variant names and the second list is the direction
                    of the split. The direction of the second list is depicted by 1's
                    and 0's. Where 1 is splitting in the direction with the variant
                    and 0 is splitting in the direction without the variant.
            split_var (string): The variant name it is now splitting on

        Returns:
            req_body (json): Returns a JSON containing the request body
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
                    }
                })
        logic['and'][0]['or'].extend(id_list)

        for variant, direction in zip(split_path[0], split_path[1]):
            if direction:
                logic['and'][1]['and'].append({"id": variant})
            else:
                logic['and'][1]['and'].append({"id": variant, "negate": True})

        # Finds index with empty list and removes it
        if logic['and'][1]['and'] == []:
            del logic['and'][1]

        # puts the request together into a form digestable by the API
        req_body = {}
        req_body['logic'] = logic
        req_body['components'] = components
        req_body['dataset_id'] = self.dataset_id
        req_body['results'] = [ {
                    "table": "patients",
                    "fields": [
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
        ancestry_counts = requests.post('%s%s' %  (self.host_url, 'count'), json=req).json()['results']['patients'][0]['ethnicity']
        if self.ancestry_list == []:
            self.ancestry_list = ancestry_counts.keys()

        return ancestry_counts

    def split_subset(self, node, split_var):
        """
        Splits the subset given the path of the splits before and a
        variable to split on.

        Attributes:
            split_path (list1, list2): 
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
        w_variant_split_path, wo_variant_split_path = GA4GH_API.create_split_path(node.split_path, split_var)

        w_var_req_body = self.craft_api_request(w_variant_split_path)
        wo_var_req_body = self.craft_api_request(wo_variant_split_path)

        # make query here for
        r_w_var = requests.post('%s%s' %  (self.host_url, 'count'), json=w_var_req_body).json()['results']['patients'][0]['ethnicity']
        r_wo_var = requests.post('%s%s' %  (self.host_url, 'count'), json=wo_var_req_body).json()['results']['patients'][0]['ethnicity']

        return r_w_var, r_wo_var


    def find_next_variant_counts(self, split_path):
        """
        Finds the counts of the a potential next variant to perform the
        split on

        Attributes:
            split_path (list1, list2): 
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

        for var in self.variant_name_list:
            if var in split_path[0]:
                w_variant_list.append({})
                continue

            split_path[0].append(var)
            split_path[1].append(1)

            req_body = self.craft_api_request(split_path)
            resp = requests.post('%s%s' %  (self.host_url, 'count'), json=req_body).json()
            variant_counts = resp['results']['patients'][0]['ethnicity'] if 'ethnicity' in resp['results']['patients'][0] else {}
            w_variant_list.append(variant_counts)

            del split_path[0][-1]
            del split_path[1][-1]
        return w_variant_list



