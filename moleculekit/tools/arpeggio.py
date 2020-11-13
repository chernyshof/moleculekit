import logging
from collections import OrderedDict
from functools import reduce
import operator

import numpy as np

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import NeighborSearch
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Polypeptide import PPBuilder

#import openbabel as ob
from openbabel import openbabel as ob

from moleculekit.config_mol import ATOM_TYPES, METALS, HALOGENS, STD_RES, PROT_ATOM_TYPES




#pdb_filename = 'KXZ_ideal.pdb'#'protein.clean.pdb'


ph = 7.4 #'pH for hydrogen addition.'
headers = True # Write out column headers in output files.
use_ambiguities = True #Turn on abiguous definitions for ambiguous contacts.


###########
# CLASSES #
###########

class HydrogenError(Exception):

    def __init__(self):
        logging.error('Please remove all hydrogens from the structure then re-run.')

class OBBioMatchError(Exception):

    def __init__(self, serial=''):

        if not serial:
            logging.error('An OpenBabel atom could not be matched to a BioPython counterpart.')

        else:
            logging.error('OpenBabel OBAtom with PDB serial number {} could not be matched to a BioPython counterpart.'.format(serial))

class AtomSerialError(Exception):

    def __init__(self):
        logging.error('One or more atom serial numbers are duplicated.')

class SiftMatchError(Exception):

    def __init__(self):
        logging.error('Seeing is not believing.')

class SelectionError(Exception):

    def __init__(self, selection):
        logging.error('Invalid selector: {}'.format(selection))





# ADDRESS AMBIGUITIES
if not use_ambiguities:

    # REMOVE IF NOT USING THE AMBIGUITIES (DEFAULT)

    # REMOVE FROM SMARTS DEFINITIONS
    ATOM_TYPES['hbond acceptor'].pop('NH2 terminal amide', None)
    ATOM_TYPES['hbond donor'].pop('oxygen amide term', None)
    ATOM_TYPES['xbond acceptor'].pop('NH2 terminal amide', None)
    ATOM_TYPES['weak hbond acceptor'].pop('NH2 terminal amide', None)

    # REMOVE FROM PROTEIN ATOM DEFINITIONS
    PROT_ATOM_TYPES['hbond acceptor'] = [x for x in PROT_ATOM_TYPES['hbond acceptor'] if x not in ('ASNND2', 'GLNNE2', 'HISCE1', 'HISCD2')]
    PROT_ATOM_TYPES['hbond donor'] = [x for x in PROT_ATOM_TYPES['hbond donor'] if x not in ('ASNOD1', 'GLNOE1', 'HISCE1', 'HISCD2')]
    PROT_ATOM_TYPES['xbond acceptor'] = [x for x in PROT_ATOM_TYPES['xbond acceptor'] if x not in ('ASNND2', 'GLNNE2', 'HISCE1', 'HISCD2')]
    PROT_ATOM_TYPES['weak hbond acceptor'] = [x for x in PROT_ATOM_TYPES['weak hbond acceptor'] if x not in ('ASNND2', 'GLNNE2', 'HISCE1', 'HISCD2')]








def process_arpeggio(pdb_filename):

    # LOAD STRUCTURE (BIOPYTHON)
    pdb_parser = PDBParser()
    s = pdb_parser.get_structure('structure', pdb_filename)
    s_atoms = list(s.get_atoms())

    logging.info('Loaded PDB structure (BioPython)')

    # CHECK FOR HYDROGENS IN THE INPUT STRUCTURE
    input_has_hydrogens = False
    hydrogens = [x for x in s_atoms if x.element == 'H']

    if hydrogens:
        logging.info('Detected that the input structure contains hydrogens. Hydrogen addition will be skipped.')
        input_has_hydrogens = True




    # LOAD STRUCTURE (OPENBABEL)
    ob_conv = ob.OBConversion()
    ob_conv.SetInFormat('pdb')
    mol = ob.OBMol()
    ob_conv.ReadFile(mol, pdb_filename)



    # CHECK THAT EACH ATOM HAS A UNIQUE SERIAL NUMBER
    all_serials = [x.serial_number for x in s_atoms]

    if len(all_serials) > len(set(all_serials)):
        raise AtomSerialError









    # MAPPING OB ATOMS TO BIOPYTHON ATOMS AND VICE VERSA

    # FIRST MAP PDB SERIAL NUMBERS TO BIOPYTHON ATOMS FOR SPEED LATER
    # THIS AVOIDS LOOPING THROUGH `s_atoms` MANY TIMES
    serial_to_bio = {x.serial_number: x for x in s_atoms}

    # DICTIONARIES FOR CONVERSIONS
    ob_to_bio = {}
    bio_to_ob = {}

    for ob_atom in ob.OBMolAtomIter(mol):

        serial = ob_atom.GetResidue().GetSerialNum(ob_atom)

        # MATCH TO THE BIOPYTHON ATOM BY SERIAL NUMBER
        try:
            biopython_atom = serial_to_bio[serial]

        except KeyError:
            # ERRORWORTHY IF WE CAN'T MATCH AN OB ATOM TO A BIOPYTHON ONE
            raise OBBioMatchError(serial)

        # `Id` IS A UNIQUE AND STABLE ID IN OPENBABEL
        # CAN RECOVER THE ATOM WITH `mol.GetAtomById(id)`
        ob_to_bio[ob_atom.GetId()] = biopython_atom
        bio_to_ob[biopython_atom] = ob_atom.GetId()

    logging.info('Mapped OB to BioPython atoms and vice-versa.')

    # ADD EMPTY DATA STRUCTURES FOR TAGGED ATOM DATA
    # IN A SINGLE ITERATION
    for atom in s_atoms:

        # FOR ATOM TYPING VIA OPENBABEL
        atom.atom_types = set([])

        # LIST FOR EACH ATOM TO STORE EXPLICIT HYDROGEN COORDINATES
        atom.h_coords = []

        # DETECT METALS
        if atom.element.upper() in METALS:
            atom.is_metal = True
        else:
            atom.is_metal = False

        # DETECT HALOGENS
        if atom.element.upper() in HALOGENS:
            atom.is_halogen = True
        else:
            atom.is_halogen = False

    # ADD EXPLICIT HYDROGEN COORDS FOR H-BONDING INTERACTIONS
    # ADDING HYDROGENS DOESN'T SEEM TO INTERFERE WITH ATOM SERIALS (THEY GET ADDED AS 0)
    # SO WE CAN STILL GET BACK TO THE PERSISTENT BIOPYTHON ATOMS THIS WAY.
    if not input_has_hydrogens:
        mol.AddHydrogens(False, True, ph) # polaronly, correctForPH, pH

        logging.info('Added hydrogens.')






    # ATOM TYPING VIA OPENBABEL
    # ITERATE OVER ATOM TYPE SMARTS DEFINITIONS
    for atom_type, smartsdict in ATOM_TYPES.items():

        #logging.info('Typing: {}'.format(atom_type))

        # FOR EACH ATOM TYPE SMARTS STRING
        for smarts in smartsdict.values():

            #logging.info('Smarts: {}'.format(smarts))

            # GET OPENBABEL ATOM MATCHES TO THE SMARTS PATTERN
            ob_smart = ob.OBSmartsPattern()
            ob_smart.Init(str(smarts))

            #logging.info('Initialised for: {}'.format(smarts))

            ob_smart.Match(mol)

            #logging.info('Matched for: {}'.format(smarts))

            matches = [x for x in ob_smart.GetMapList()]

            #logging.info('List comp matches: {}'.format(smarts))

            if matches:

                # REDUCE TO A SINGLE LIST
                matches = set(reduce(operator.add, matches))

                #logging.info('Set reduce matches: {}'.format(smarts))

                for match in matches:

                    atom = mol.GetAtom(match)
                    ob_to_bio[atom.GetId()].atom_types.add(atom_type)

                #logging.info('Assigned types: {}'.format(smarts))


    # ALL WATER MOLECULES ARE HYDROGEN BOND DONORS AND ACCEPTORS
    for atom in (x for x in s_atoms if x.get_full_id()[3][0] == 'W'):
        atom.atom_types.add('hbond acceptor')
        atom.atom_types.add('hbond donor')

    # OVERRIDE PROTEIN ATOM TYPING FROM DICTIONARY
    for residue in s.get_residues():

        if residue.resname in STD_RES:

            for atom in residue.child_list:

                # REMOVE TYPES IF ALREADY ASSIGNED FROM SMARTS
                for atom_type in PROT_ATOM_TYPES.keys():
                    atom.atom_types.discard(atom_type)

                # ADD ATOM TYPES FROM DICTIONARY
                for atom_type, atom_ids in PROT_ATOM_TYPES.items():

                    atom_id = residue.resname.strip() + atom.name.strip()

                    if atom_id in atom_ids:
                        atom.atom_types.add(atom_type)





    def make_pymol_string(entity):
        '''
        Feed me a BioPython atom or BioPython residue.
        See `http://pymol.sourceforge.net/newman/user/S0220commands.html`.
        chain-identifier/resi-identifier/name-identifier
        chain-identifier/resi-identifier/
        '''

        if isinstance(entity, Atom):

            chain = entity.get_parent().get_parent()
            residue = entity.get_parent()
            atom_name = entity.name

        elif isinstance(entity, Residue):
            chain = entity.get_parent()
            residue = entity
            atom_name = ''

        else:
            raise TypeError('Cannot make a PyMOL string from a non-Atom or Residue object.')

        res_num = residue.id[1]

        # ADD INSERTION CODE IF NEED BE
        if residue.id[2] != ' ':
            res_num = str(res_num) + residue.id[2]

        macro =  '{}/{}/{}'.format(chain.id,
                                res_num,
                                atom_name)

        return macro

    '''

    with open(pdb_filename.replace('.pdb', '.atomtypes'), 'w') as fo:

        if headers:
            fo.write('{}\n'.format('\t'.join(
                ['atom', 'atom_types']
            )))

        for atom in s_atoms:
            fo.write('{}\n'.format('\t'.join([str(x) for x in [make_pymol_string(atom), sorted(tuple(atom.atom_types))]])))

    logging.info('Typed atoms.')


    '''

    return s_atoms



