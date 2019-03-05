from __future__ import division
from ID3_Class import ID3
from local_API import LOCAL_API
from ga4gh_API import GA4GH_API
from ID3_Node import ID3_Node

class ConfusionMatrix(ID3):

    def __init__(self, file_path):
        '''
        Creates confusion matrix with the first index (Y) as the correct population and the
        second index (X) as the predicted population. THe order of the ancestries is dictated
        by the ancestry_list within the api class

        actual is Y axis
        predicted is X axis

        Args:
            file_path (str): Path to json file that contains the variant rangess

        Attributes:
            api (API): API object that is used to interact with the virtual API
            root_node (Node): Creates the root node of the tree to be added upon

            length (int): length of all the ancestries
            ancestry_list (list): list of all ancestries without repetitions
            conf_matrix (list): the confusion matrix based on the ID3 classifier
            diagonal_sum (int): sum of the diagonals within the matrix
            total (int): total sum of all values in matrix
        '''

        # modify API for conf_matrix
        self.api = LOCAL_API(file_path, conf_matrix=True)

        # initialize ID3 algorithm
        subset = self.api.get_target_set()
        self.root_node = ID3_Node('root', subset, True)
        self.ID3(self.root_node)

        # create conf_matrix and calculate useful attributes
        self.length = len(self.api.ancestry_list)
        self.conf_matrix = [ [0 for x in range(0, self.length)] for y in range(0, self.length) ]

        for variants, popu in zip(self.api.test_variant_list, self.api.test_popu_list):
            include_variants = [self.api.variant_name_list[idx] for idx, is_variant in enumerate(variants) if is_variant == 1]

            # Actual Result
            y = self.api.ancestry_list.index( popu )
            # Predicted Result
            x = self.api.ancestry_list.index( self.predict(include_variants).most_common_ancestry )

            self.conf_matrix[y][x] += 1

        self.diagonal_sum = sum( [self.conf_matrix[i][i] for i in range(0, self.length)] )
        self.total = sum( [sum(self.conf_matrix[i]) for i in range(0, self.length)] )

    def accuracy(self):
        '''
        How often the classifier is correct

        Returns:
            (float): a number between 0 and 1
        '''
        return self.diagonal_sum / self.total

    def misclass_rate(self):
        '''
        How often is the classifier wrong

        Returns:
            (float): a number between 0 and 1
        '''
        return 1 - self.accuracy()

    def true_ancestry_rate(self, ancestry):
        '''
        For a particular ancestry, how often does it predict the correct ancestry

        Returns:
            (float): a number between 0 and 1

        '''
        if ancestry not in self.api.ancestry_list:
            return "Not a valid ancestry"

        popu_i = self.api.ancestry_list.index(ancestry)

        true_anc = self.conf_matrix[popu_i][popu_i]
        sum_actual_anc = sum([ self.conf_matrix[popu_i][i] for i in range(0, self.length) ])

    def false_ancestry_rate(self, ancestry):
        """
        For a particular ancestry, how often does it predict the incorrect ancestry

        Args:
            ancestry (string): a three character code depicting populations of people

        Returns:
            (float): a number between 0 and 1

        """
        return 1 - self.true_ancestry_rate

    def precision(self, ancestry):
        '''
        For a particular ancestry, how often is the prediction correct

        Args:
            ancestry (string): a three character code depicting populations of people

        Returns:
            (float): a number between 0 and 1
        '''
        if ancestry not in self.api.ancestry_list:
            return "Not a valid ancestry"

        popu_i = self.api.ancestry_list.index(ancestry)
        true_anc = self.conf_matrix[popu_i][popu_i]
        sum_pred_anc = sum([ self.conf_matrix[i][popu_i] for i in range(0, self.length) ])

        return true_anc/sum_pred_anc

    def prevalance(self, ancestry):
        '''
        How often does the ancestry appear in the sample relative to the sum of all ancestries

        Args:
            ancestry (string): a three character code depicting populations of people

        Returns:
            (float): a number between 0 and 1
        '''
        if ancestry not in self.api.ancestry_list:
            return "Not a valid ancestry"

        popu_i = self.api.ancestry_list.index(ancestry)
        sum_actual_anc = sum([ self.conf_matrix[popu_i][i] for i in range(0, self.length) ])

        return sum_actual_anc / self.total

    def print_matrix(self):
        for row in self.conf_matrix:
            print row

if __name__ == "__main__":
    c = ConfusionMatrix('test_variant_ranges.json')
    c.print_matrix()
    print c.accuracy()
    print c.prevalance('ESN')


