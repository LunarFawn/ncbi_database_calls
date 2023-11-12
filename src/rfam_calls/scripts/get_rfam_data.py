"""
File to sandbox getting rfam data for RNA hub
"""

import json
import requests


from typing import Dict, Any




def post_sequence_json():
    url = 'https://rfam.org/search/sequence'
    headers = {'Accept': 'application/json'}
    payload = {'seq': "GCGGCCGGTGATGAGAACTTCTCCCACTCACATTCGAGTTTCCCGACCATGAGATGACTCCACATGCACTACCATCTGAGGCCAC"}
    r = requests.post(url=url, headers=headers, data=payload)
    print (r.json())
    
def get_sequence_json():
    url = 'https://rfam.org/search/sequence/3F2368D2-81AA-11EE-95C6-1FD6CE29BD2E'
    headers = {'Accept': 'application/json'}
    #payload = {'seq': "AGTTACGGCCATACCTCAGAGAATATACCGTATCCCGTTCGATCTGCGAAGTTAAGCTCTGAAGGGCGTCGTCAGTACTATAGTGGGTGACCATATGGGAATACGACGTGCTGTAGCTT"}
    r = requests.get(url=url, headers=headers)
    # print (r.json()["hits"])
    
      
    hits_dict:Dict[Any,Any] = r.json()["hits"]
    
    for hit in hits_dict.keys():
        print(f'Hit = {hit}')
        new_thing = hits_dict[hit]
        for key, value in new_thing[0]['alignment'].items():
            print(f'{key} = {value}')
            # print (new_thing[0]['alignment']) 


get_sequence_json()

