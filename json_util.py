"""
Utility functions for IO opteration in JSON format
"""

import json


def write_json(data, json_file='data_input.json'):
    with open(json_file, 'w') as j_out:
        json.dump(data, j_out)
    print(f"write to {json_file}")


def read_json(json_file= 'data_input.json'):
    with open(json_file, 'r') as j_file:
        data = json.load(j_file)
    print(f"read from {json_file}")
    return data
