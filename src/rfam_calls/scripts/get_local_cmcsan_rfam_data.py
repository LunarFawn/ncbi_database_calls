"""
File to get rfam data via local cmscan search
"""



import subprocess
import os
import gzip
from typing import List

def find_alignment_from_seed(seed_path:str):
    with gzip.open(seed_path,'rt', encoding='ISO-8859-1') as f:
        for line in f:
            print('got line', line)

def run_cmscan(db_path:str, fasta_path:str):
    command:str = f'cmscan --noali {db_path} {fasta_path}'
    result = os.popen(cmd=command)
    
    full_result:List[str] = result.read().split('\n')
    
    alignment_hits_modelnames:List[str] = []
    
    in_hits:bool = False
    exited_hits:bool = False
    start_hits:bool = False
    start_index:int = -1
    end_index:int = -1
    
    hit_count:int = -1
    
    for line in full_result:
        
        if exited_hits is True:

            if 'Total CM and HMM hits reported:' in line:
                number_string:str = line[31:].replace(' ', '')
                hit_count = int(number_string)
                print (f'hits = {alignment_hits_modelnames}')
                print(f'hit count = {hit_count}')
                exited_hits  = False
        
        if in_hits is True:
            if '----' not in line and line != '':
                temp_string:str = line[start_index:]
                end_index = start_index
                for char in temp_string:
                    if char == " ":
                        break
                    end_index += 1
                hit_name:str = line[start_index:end_index]
                alignment_hits_modelnames.append(hit_name)
            if line == '':
                exited_hits = True
                in_hits = False
                
                
        if start_hits is True:
            if "modelname" in line:
                start_index = line.index("modelname")
                in_hits = True
                start_hits = False
        
        if "Hit scores:" in line:
            start_hits = True
    
    if hit_count == len(alignment_hits_modelnames):
        return alignment_hits_modelnames
    else:
        raise Exception("number of hits does not match hits collected")
             
    # find_alignment_from_seed(seed_path='/home/rnahub/rnahub/db/rfam/Rfam.seed.gz')    

    
run_cmscan(db_path='/home/rnahub/rnahub/db/rfam/Rfam.cm',
           fasta_path='/home/pearljen/repo/ncbi_database_calls/src/rfam_calls/scripts/example/example.fa')
