"""
File to get rfam data via local cmscan search
"""



import subprocess
import os


def run_cmscan(db_path:str, fasta_path:str):
    command:str = f'cmscan --noali {db_path} {fasta_path}'
    result = os.popen(cmd=command)
    print(result.read())
    
def find_alignment_from_seed():
    pass
    
run_cmscan(db_path='/home/rnahub/rnahub/db/rfam/Rfam.cm',
           fasta_path='/home/pearljen/repo/ncbi_database_calls/src/rfam_calls/scripts/example/example.fa')
