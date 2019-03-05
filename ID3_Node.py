from anytree import NodeMixin

class ID3_Node(NodeMixin):
    def __init__(self, variant_name, subset, with_variant, split_path=([], []), parent=None, children=None):
        super(ID3_Node, self).__init__()
        self.variant_name = variant_name
        self.with_variant = with_variant
        self.subset = subset
        self.total_count = str(sum(subset.values()))
        self.most_common_ancestry = str(max(subset, key=subset.get))
        self.parent = parent
        self.split_path = split_path
        if children:
            self.children = children

    @staticmethod
    def name_func(node):
        return "%s: %s \n %s" % ( 'with' if node.with_variant else 'w/o', node.variant_name, ID3_Node.nodeattrfunc(node))

    @staticmethod
    def nodeattrfunc(node):
        return "most common ancestry: %s | total count: %s" % (node.most_common_ancestry, node.total_count)