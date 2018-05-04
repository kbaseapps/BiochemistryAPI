import csv
import os
from collections import defaultdict, namedtuple, OrderedDict

from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.DataStructs import FingerprintSimilarity

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
# TODO: add logging


def check_param(in_params, req_param, opt_param=list()):
    """
    Check if each of the params in the list are in the input params
    """
    for param in req_param:
        if param not in in_params:
            raise ValueError('{} parameter is required'.format(param))
    defined_param = set(req_param + opt_param)
    for param in in_params:
        if param not in defined_param:
            print("WARNING: received unexpected parameter {}".format(param))


def dict_from_file(path, key='id', dialect='excel-tab'):
    """
    Build a dictionary from an object array in a file
    :param path: local path to object
    :param key: what field should be used as the key
    :param dialect: excel-tab for TSV or excel for CSV
    :return:
    """
    if not os.path.exists(path):
        raise ValueError("File not found: {}".format(path))
    reader = csv.DictReader(open(path), dialect=dialect)
    return OrderedDict([(x[key], x) for x in reader])


def alias_dict_from_file(path, dialect='excel-tab'):
    """
    Build a dictionary from an object array in a file
    :param path: local path to object
    :param dialect: excel-tab for TSV or excel for CSV
    :return:
    """
    alias_mappings = defaultdict(list)
    with open(path) as infile:
        r = csv.DictReader(infile, dialect=dialect)
        for line in r:
            for seed_id in line['MS ID'].split('|'):
                if line['Source'] == 'Enzyme Class':
                    alias_mappings[seed_id].append(line['External ID'])
                else:
                    alias_mappings[seed_id].append('%s:%s' % (
                        line['Source'].strip(), line['External ID']))
    return alias_mappings


def make_mol_tuples(compound_dict, id_key="id", struct_key='structure', struct_type='inchi'):
    """
    Creates named tuples with (compound_id, RDKit Mol Object) from a dict with SMILES or InChI
    """
    MolTuple = namedtuple("MolTuple", ("id", "mol"))
    if struct_type.lower() == 'inchi':
        return [MolTuple(comp[id_key], AllChem.MolFromInchi(comp[struct_key]))
                for comp in compound_dict]
    if struct_type.lower() == 'smiles':
        return [MolTuple(comp[id_key], AllChem.MolFromSmiles(comp[struct_key]))
                for comp in compound_dict]


def substructure_search(query, structures):
    """Performs substructure search on 'query' in 'structures'"""
    pattern = AllChem.MolFromSmarts(query)
    return [x.id for x in structures if x.mol and x.mol.HasSubstructMatch(pattern)]


def similarity_search(query, structures, fp_type='maccs', min_similarity=0.8):
    """Perform return compound ids where tanimoto similarity of fingerprint 'fp_type is greater
    than 'min_similarity'"""
    if fp_type.lower() == 'maccs':
        fp1 = AllChem.GetMACCSKeysFingerprint(_get_mol(query))
        return [x.id for x in structures
                if x.mol and FingerprintSimilarity(fp1, AllChem.GetMACCSKeysFingerprint(x.mol))
                >= min_similarity]
    elif fp_type.lower() == 'rdkit':
        fp1 = AllChem.RDKFingerprint(_get_mol(query))
        return [x.id for x in structures
                if x.mol and FingerprintSimilarity(fp1, AllChem.RDKFingerprint(x.mol))
                >= min_similarity]
    else:
        fp_types = ", ".join(('maccs', 'rdkit'))
        raise ValueError('Invalid fingerprint type: choose one of {}'.format(fp_types))


def depict_compound(structure, size=(300, 300)):
    """
    Generate a SVG depiction of a chemical structure
    :param structure: SMILES or InChI string
    :param size: Tuple of (Height, Width)
    :return:
    """
    mol = _get_mol(structure)
    if not mol:
        return ""
    AllChem.Compute2DCoords(mol)
    dwr = rdMolDraw2D.MolDraw2DSVG(*size)
    dwr.DrawMolecule(mol)
    dwr.FinishDrawing()
    return dwr.GetDrawingText().replace('svg:', '')


def _get_mol(structure):
    if "InChI" in structure:
        mol = AllChem.MolFromInchi(structure)
    else:
        mol = AllChem.MolFromSmiles(structure)
    return mol
