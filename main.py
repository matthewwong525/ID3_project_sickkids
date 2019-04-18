from ID3_Class import ID3

# Creates ID3 object with a filepath to the config.json
id3_obj = ID3('config.json', local=False)

# prints the ID3 tree as `tree.png`
id3_obj.print_tree('tree')

# prints the list of all the variant names
print(id3_obj.api.variant_name_list)

# predicts ancestry of the person with the variant `22:50121766:50121767` and no other variant in `variant_name_list`
id3_obj.predict(['22:50121766:50121767'])