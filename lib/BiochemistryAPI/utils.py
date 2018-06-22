import csv
import logging
import os
from collections import defaultdict, namedtuple, OrderedDict

from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.DataStructs import FingerprintSimilarity

rdk_lg = RDLogger.logger()
rdk_lg.setLevel(RDLogger.CRITICAL)
logging.basicConfig(level=logging.INFO)


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
            logging.warning("Received unexpected parameter {}".format(param))


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


def _get_mol(structure):
    if "InChI" in structure:
        mol = AllChem.MolFromInchi(structure)
    else:
        mol = AllChem.MolFromSmiles(structure)
    return mol


def make_mol_tuples(compound_dict, id_key="id", struct_key='structure', struct_type='inchi'):
    """
    Creates named tuples with (compound_id, RDKit Mol Object) from a dict with SMILES or InChI
    """
    MolTuple = namedtuple("MolTuple", "id mol maccs_fp rdkit_fp")
    tups = []
    for comp in compound_dict:
        mol = _get_mol(comp[struct_key])
        if mol:
            tups.append(MolTuple(comp[id_key],
                                 mol,
                                 AllChem.GetMACCSKeysFingerprint(mol),
                                 AllChem.RDKFingerprint(mol)))
    return tups


def substructure_search(query, structures):
    """Performs substructure search on 'query' in 'structures'"""
    pattern = AllChem.MolFromSmarts(query)
    return [x.id for x in structures if x.mol and x.mol.HasSubstructMatch(pattern)]


def similarity_search(query, structures, fp_type='maccs', min_similarity=0.8):
    """Perform return compound ids where tanimoto similarity of fingerprint 'fp_type is greater
    than 'min_similarity'"""
    if not isinstance(min_similarity, float) or not 0 <= min_similarity <= 1.0:
        raise ValueError('Invalid min_similarity. Value must be a float between 0 and 1.')
    if fp_type.lower() == 'maccs':
        fp1 = AllChem.GetMACCSKeysFingerprint(_get_mol(query))
        return [x.id for x in structures
                if FingerprintSimilarity(fp1, x.maccs_fp) >= min_similarity]
    elif fp_type.lower() == 'rdkit':
        fp1 = AllChem.RDKFingerprint(_get_mol(query))
        return [x.id for x in structures
                if FingerprintSimilarity(fp1, x.rdkit_fp) >= min_similarity]
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


def get_3d_mol(structure):
    mol = _get_mol(structure)
    if not mol:
        return ""
    AllChem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.RemoveHs(mol)
    return AllChem.MolToMolBlock(mol)
