def build_module_label(dataset, module, name, columns, path, key):
    print("Building labels named {} for dataset {} for columns {} in {} on key={} based on module {}".format(
        name, dataset, columns, path, key, module))

def build_custom_label(dataset, name, columns, path, key):
    print("Building custom labels named {} for dataset {} for columns {} in {} on key={}".format(
        name, dataset, columns, path, key))
