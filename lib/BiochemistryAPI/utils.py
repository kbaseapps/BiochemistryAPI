from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


def depict_compound(structure, size=(300, 300)):
    if "InChI" in structure:
        mol = AllChem.MolFromInchi(structure)
    else:
        mol = AllChem.MolFromSmiles(structure)
    if not mol:
        return ""
    AllChem.Compute2DCoords(mol)
    dwr = rdMolDraw2D.MolDraw2DSVG(*size)
    dwr.DrawMolecule(mol)
    dwr.FinishDrawing()
    return dwr.GetDrawingText().replace('svg:', '')
