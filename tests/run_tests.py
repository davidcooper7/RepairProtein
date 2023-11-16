"""
Run tests to make sure that RepairProtein class if functioning correctly.

Usage:
>python run_tests.py
"""


import os, sys, shutil
sys.path.append('../')
from RepairProtein import RepairProtein

# Input files 
working_dir = os.getcwd() + '/'
prot_dir = working_dir + 'test_input'
fasta_fn = working_dir + '4OBO.fasta'

# Make intermediate directory 
int_dir = working_dir + 'modeller_intermediates'
if not os.path.exists(int_dir):
    os.mkdir(int_dir)

# Make output directory
prot_out_dir = working_dir + 'test_output'
if not os.path.exists(prot_out_dir):
    os.mkdir(prot_out_dir)

# Iterate through input files
try:
    for pdb in os.listdir(prot_dir):
        rp = RepairProtein(pdb_fn=prot_dir + '/' + pdb,
                        fasta_fn=fasta_fn,
                        working_dir=int_dir)
        rp.run(pdb_out_fn=prot_out_dir + '/' + pdb)

    print('All tests passed')

# Check where RepairProtein failed
except:
    print(os.getcwd())
    print('\n!!!Tests failed!!!\n')
    print('Running diagnostics...')

    # Check for .pir files
    print('Testing for template .pir file...')
    if not os.path.exists(rp.temp_pir_fn):
        print('\tTest failed: ' + rp.temp_pir_fn + ' does not exist.' )
    else:
        print('\tTest passed.')

    print('Testing for target .pir file...')
    if not os.path.exists(rp.tar_pir_fn):
        print('\tTest failed: ' + rp.tar_pir_fn + ' does not exist.' )
    else:
        print('\tTest passed.')

    # Check for missing residues
    print('Testing for missing residues...')
    if len(rp.missing_residues) != 39:
        print('\tTest failed: number of missing residues (' + len(rp.missing_residues) + ') should equal 39.' )
    else:
        print('\tTest passed.')

    # Check for .ali files
    print('Testing for .ali file...')
    if not os.path.exists(rp.ali_fn):
        print('\tTest failed: ' + rp.ali_fn + ' does not exist.' )
    else:
        print('\tTest passed.')  

    # Check for output files
    print('Testing for output file...')
    if not os.path.exists(rp.pdb_out_fn):
        print('\tTest failed: ' + rp.pdb_out_fn + ' does not exist.' )
    else:
        print('\tTest passed.')  

    print('Traceback log\n\n')

    




