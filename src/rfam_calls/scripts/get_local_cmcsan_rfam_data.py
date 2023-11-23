"""
File to get rfam data via local cmscan search
"""



import subprocess
import os
import gzip
from typing import List, Dict
from dataclasses import dataclass, field
import heapq
from attr import define


# @define(kw_only=True)
# class SeedAlignmentResult():
#     name:str
#     alignment:str

@define(kw_only=True)
class CmscanResult():
    model_name:str
    score:float
    comment:str
    is_error:bool = False
    seed_alignments:List[str] = []
    seed_alignemnt_dict:Dict[str, str] = {}
    

def run_cmscan(db_path:str, fasta_path:str)-> List[CmscanResult]:
    command:str = f'cmscan --noali {db_path} {fasta_path}'
    result = os.popen(cmd=command)
    raw_result:str = result.read()
    print(raw_result)
    
    full_result:List[str] =raw_result.split('\n')
    
    
    
    alignment_hits_modelnames:List[str] = []
    alignment_hits:Dict[str, CmscanResult] = {}
    
    in_hits:bool = False
    exited_hits:bool = False
    start_hits:bool = False
    
    score_word:str = 'score'
    score_start_index:int = -1
    score_end_index:int = -1
    
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
                is_error:bool = False
                comment:str = ''
                #first get modelname for alignment
                temp_string:str = line[start_index:]
                end_index = start_index
                for char in temp_string:
                    if char == " ":
                        break
                    end_index += 1
                hit_name:str = line[start_index:end_index]
                # alignment_hits_modelnames.append(hit_name)
                
                #now get score by working backwards until a space is found
                temp_string = line[:score_end_index]
                score_start_index = score_end_index
                for index in range(len(temp_string)-1, -1, -1):
                    if temp_string[index] == " ":
                        break
                    score_start_index -= 1
                score_string:str = line[score_start_index:score_end_index]
                score:float = 0
                try:
                    score = float(score_string)
                except:
                    is_error = True
                    comment = "unable to get float from score"
                    score = -1
                
                new_result:CmscanResult = CmscanResult(model_name=hit_name,
                                                       score=score,
                                                       comment=comment,
                                                       is_error=is_error)
                alignment_hits[hit_name] =new_result
                
            if line == '':
                exited_hits = True
                in_hits = False
                
                
        if start_hits is True:
            #get modelname index so that we know the start index of the name of
            #the hit that we will do alignment search against
            if "modelname" in line:
                start_index = line.index("modelname")
                in_hits = True
                start_hits = False
                
            #now get end-index of score as we will move backwards to get
            #start-index later
            if score_word in line:
                score_end_index = line.index(score_word) + len(score_word)
            
        if "Hit scores:" in line:
            start_hits = True
    
    if hit_count == len(alignment_hits):
        return alignment_hits
    else:
        raise Exception(f'number of hits {hit_count} does not match hits collected {alignment_hits}')


def get_alignment_data_from_seed():
    pass

def find_alignment_from_seed(cmscan_hits:Dict[str, CmscanResult], seed_path:str,):
    #make a heapq list que that can be poped as items are found
    hit_name_list:List[str] = list(cmscan_hits.keys())
   
    stockhold_1_0_id:str = '#=GF ID'
    stockhold_1_0_gap:str = '   '
    
    in_cmscan_hit_now:bool = False
    found_empty_lines:bool = False
    
    current_found_model_name:str = ''
    current_found_alignments:List[str] = []
    
    current_found_alignments_dict:Dict[str, str] = {}
    
    with gzip.open(seed_path,'rt', encoding='ISO-8859-1') as file:
        for line in file:
            
            if found_empty_lines is True:
                if line != '\n':
                    if line == '//\n':    
                        #this is teh end of the aligments so
                        #reset and save the data
                        found_empty_lines = False
                        current_hit:CmscanResult = cmscan_hits[current_found_model_name]
                        current_hit.seed_alignments = current_found_alignments
                        
                        current_hit.seed_alignemnt_dict = current_found_alignments_dict
                        cmscan_hits[current_found_model_name] = current_hit
                        #now reset the alignments
                        current_found_alignments = []
                        current_found_alignments_dict = {}
                        hit_name_list.remove(current_hit.model_name)
                        if len(hit_name_list) > 0:
                            continue
                        else:
                            break
                    
                    current_found_alignments.append(line)
                    
                    end_name_index:int = 0
                    start_alignment_index:int = -1
                    
                    temp_line:str = line
                    
                    special_string:str = "#=GC"

                    if special_string in line:
                        special_index = line.index(special_string) + len(special_string) + 1
                        temp_line = line[special_index:]
                    
                    for index in range(len(temp_line)):                       
                        
                        if temp_line[index] != " ":
                            end_name_index += 1
                        else:
                            start_alignment_index = index
                            break
                    
                    alignment_name:str = temp_line[:end_name_index]
                    if special_string in line:
                        alignment_name = f'#=GC {alignment_name}'
                    alignment_string:str = temp_line[start_alignment_index:].strip(' ').strip('\n')
                    current_found_alignments_dict[alignment_name] = alignment_string
                    
                    
                    
            
            if in_cmscan_hit_now is True:
                if line == '\n':
                    found_empty_lines = True
                    in_cmscan_hit_now = False
            
            
            if stockhold_1_0_id in line:
                #get name of this entry
                split_line:List[str] = line.split(stockhold_1_0_gap)
                model_name:str = split_line[1].strip('\n')
                if model_name in hit_name_list:   
                    current_found_model_name = model_name                 
                    print(current_found_model_name)
                    in_cmscan_hit_now = True
                # begin_seak_index = line.index(stockhold_1_0_id) + len(stockhold_1_0_id)
                # temp_string:str = line[begin_seak_index:]
                # name_start_index:int = -1
                # name_end_index:int = -1
                # for seak_index in range(len(temp_string)):
                    
            # print('got line', line)    
    return cmscan_hits     
    
results:CmscanResult = run_cmscan(db_path='/home/rnahub/rnahub/db/rfam/Rfam.cm',
                                fasta_path='/home/pearljen/repo/ncbi_database_calls/src/rfam_calls/scripts/example/example.fa')
print (results)

seed_resutls:Dict[str, CmscanResult] = find_alignment_from_seed(seed_path='/home/rnahub/rnahub/db/rfam/Rfam.seed.gz',
                         cmscan_hits=results)

print()
print('Entero_5_CRE')
# print(seed_resutls['Entero_5_CRE'].seed_alignemnt_dict)
print("Alignments")
print()
for item in seed_resutls['Entero_5_CRE'].seed_alignemnt_dict:
    print (f'name = {item}, Alignment = {seed_resutls["Entero_5_CRE"].seed_alignemnt_dict[item]}')
print()
print('MAT2A_D')
# print(seed_resutls['MAT2A_D'].seed_alignemnt_dict)
for item in seed_resutls['MAT2A_D'].seed_alignemnt_dict:
    print (f'name = {item}, Alignment = {seed_resutls["MAT2A_D"].seed_alignemnt_dict[item]}')
