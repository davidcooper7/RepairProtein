
"""
TODO: Add chain options, thermodynamic mutations, unnatural residues?
"""
import shutil, os
from modeller import *
from modeller.automodel import *
from file_handling import parse_sequence, fasta_to_pir, pdb_to_pir
from SequenceWrapper import SequenceWrapper as SeqWrap
import mdtraj as md
import numpy as np


class RepairProtein():
    """
    Repair a broken protein 

    Attributes:
    -----------
        pdb_fn (str):
            String path to .pdb file to repair.

        fasta_fn (str):
            String path to .fasta file that contains sequence to use as a template to repair the protein .pdb.     

        working_dir (str):
            String path to working directory where all intermediate files made by UCSF modeller will be stored. Default is current working directory. 

        name (str):
            Name of pdb without .pdb extension

        pdb_out_fn (str):
            String path to write repaired .pdb file. 


    
    """

    def __init__(self, pdb_fn: str, fasta_fn: str, working_dir: str='./'):
        """
        Initialize RepairProtein object.

        Parameters:
        -----------
            pdb_fn (str):
                String path to .pdb file to repair.
            
            fasta_fn (str):
                String path to .fasta file that contains sequence to use as a template to repair the protein .pdb.     

            working_dir (str):
                String path to working directory where all intermediate files made by UCSF modeller will be stored. Default is current working directory. 

            temp_seq (str):
                String of template sequence.

            tar_seq (str):
                String of target sequence.

            missing_residues (np.array): 
                2-D Array of indices of missing sequence entries and what the missing entry is. e.g. [[0, 'A'], [1, 'B'], [4, 'E'], ....]

            term_residues (List):
                List of indices that represent the index N and C terminus of target sequence in template sequence. e.g. [4, 345]
        """

        # Initialize variables
        self.pdb_fn = pdb_fn
        self.fasta_fn = fasta_fn
        self.working_dir = working_dir
        self.name = self.pdb_fn.split('.pdb')[0]
        try:
            self.name = self.name.split('/')[-1]
        except:
            pass


    def run(self, pdb_out_fn: str, tails: bool=False):
        """
        Run the remodelling.

        Parameters:
        -----------
            pdb_out_fn (str):
                String path to write repaired .pdb file. 

            tails (bool):
                If True, add missing residues to N and C termini. Default is False.
        """
        # Parse template sequence from .fasta
        self.temp_pir_fn = self.working_dir + '/' + self.name + '.pir'
        fasta_to_pir(self.fasta_fn, self.temp_pir_fn)
        self.temp_seq = parse_sequence(self.temp_pir_fn)

        # Parse target sequence from .pdb
        self.tar_pir_fn = self.working_dir + '/' + self.name + '.pir'
        shutil.copy(self.pdb_fn, self.working_dir + '/' + self.pdb_fn.split('/')[-1])
        pdb_to_pir(self.name, self.working_dir)
        self.tar_seq = parse_sequence(self.tar_pir_fn)

        # Find missing residues
        sw = SeqWrap(self.temp_seq, self.tar_seq)
        sw.find_missing_residues(verbose=False, tails=tails)
        self.missing_residues = sw.missing_residues
        self.term_residues = sw.term_residues
        self.ali_fn = self.working_dir + '/' + self.name + '.ali'
        sw.write_alignment_file(self.ali_fn, self.temp_pir_fn)

        # Model 
        env = Environ()
        env.io.atom_files_directory = ['.', self.working_dir]

        class MyModel(AutoModel):
            def select_atoms(self):
                return Selection(self.residue_range('1:A', '1:A'))
            
        model = MyModel(env, alnfile = self.ali_fn, knowns=self.name, sequence=self.name + '_fill')
        model.starting_model = 1
        model.ending_model = 1

        cwd = os.getcwd()
        os.chdir(self.working_dir)
        model.make()
        os.chdir(cwd)

        # Move UCSF modeller output to desired location
        self.pdb_out_fn = pdb_out_fn
        modeller_pdb_out_fn = self.working_dir + '/' + model.outputs[0]['name']
        shutil.copy(modeller_pdb_out_fn, self.pdb_out_fn)

        # Delete tails if necessary
        if not tails:
            traj = md.load_pdb(self.pdb_out_fn)
            top = traj.topology
            print('!!!' + str(self.term_residues))
            resid_range = ' '.join(str(i) for i in range(self.term_residues[0], self.term_residues[1]))
            print('!!!' + resid_range)
            sele = top.select(f'resid {resid_range}')
            traj = traj.atom_slice(sele)
            traj.save_pdb(self.pdb_out_fn)





