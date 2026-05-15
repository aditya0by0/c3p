"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: CHEBI:37664 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sterol backbone pattern (tetracyclic ring system with a 3-hydroxy group)
    sterol_pattern = Chem.MolFromSmarts("[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C")
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol backbone found"

    # Look for ester group (-O-C(=O)-) attached to the 3-hydroxy group
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester group found"

    # Check if the ester group is attached to the 3-hydroxy group of the sterol
    ester_attached_to_sterol = False
    for match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in ester_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                for neighbor_neighbor in neighbor.GetNeighbors():
                    if neighbor_neighbor.GetAtomicNum() == 8 and neighbor_neighbor.GetIdx() in [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetIdx() in [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetIdx() in [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]]]:
                        ester_attached_to_sterol = True
                        break
                if ester_attached_to_sterol:
                    break
        if ester_attached_to_sterol:
            break

    if not ester_attached_to_sterol:
        return False, "Ester group not attached to the 3-hydroxy group of the sterol"

    return True, "Contains sterol backbone with ester linkage at the 3-hydroxy group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35915',
                          'name': 'sterol ester',
                          'definition': 'A steroid ester obtained by formal '
                                        'condensation of the carboxy group of '
                                        'any carboxylic acid with the '
                                        '3-hydroxy group of a sterol.',
                          'parents': ['CHEBI:33308', 'CHEBI:47880'],
                          'xrefs': ['KEGG:C01958'],
                          'all_positive_examples': []},
    'config': None,
    'code_statistics': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O(C1C(O)C(OC(OC2=C(OC=3C(C2=O)=C(O)C=C(O)C3)C4=CC(O)=C(O)C=C4)C1O)CO)C5OC(C(O)C(O)C5O)CO',
                                     'name': 'Quercetin 3-beta-laminaribioside',
                                     'reason': 'No sterol backbone found'},
                                 {   'smiles': 'O=C(N[C@@H](CC1=CC=C(O)C=C1)C(O)=O)[C@H]2N(CCC2)C(=O)[C@@H](N)CC3=CC=CC=C3',
                                     'name': 'Phe-Pro-Tyr',
                                     'reason': 'No sterol backbone found'},
                                 {   'smiles': 'O[C@@H]1[C@@H](C(=C2C[C@](CO)(C)C[C@H]2[C@H](C1)C)CO)CO',
                                     'name': '4alpha,11,12,14-tetrahydroxy-1-tremulene',
                                     'reason': 'No sterol backbone found'},
                                 {   'smiles': 'O(CCCCCCCCCC\\C=C/C=C\\CC)C(=O)C',
                                     'name': '11Z,13Z-Hexadecadienyl acetate',
                                     'reason': 'No sterol backbone found'},
                                 {   'smiles': 'O=C(NC(C(=O)NC(CN(CCO)C)C)(C)C)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C1N(C(=O)C(CCCCCCCC)C)CCC1)CC(CC(O)CC(=O)CC)C)C)(C)C)C(CC)C)C(C)C)(C)C',
                                     'name': 'Roseoferin A1',
                                     'reason': 'No sterol backbone found'},
                                 {   'smiles': 'O=C1C(=C(CC2=CC=C(O)C=C2)C(C1(C3=CC=C(O)C=C3)C4(C5=CC=C(O)C=C5)C(=O)C(CC6=CC=C(O)C=C6)=C(C4=O)CC7=CC=C(O)C=C7)=O)CC8=CC=C(O)C=C8',
                                     'name': 'Nostotrebin 6',
                                     'reason': 'No sterol backbone found'},
                                 {   'smiles': 'ClC1=C(O)C=CC(=C1)NC2=NC=NC=3C2=CC(OCCCNCCCO)=C(C3)OC',
                                     'name': '2-chloro-4-[(6-{3-[(3-hydroxypropyl)amino]propoxy}-7-methoxyquinazolin-4-yl)amino]phenol',
                                     'reason': 'No sterol backbone found'},
                                 {   'smiles': 'C1[C@@H]2C([C@@H](N2)CN1CC3=COC=N3)C4=CC=C(C=C4)C5=CC=C(C=C5)C#N',
                                     'name': '4-[4-[(1S,5R)-3-(4-oxazolylmethyl)-3,6-diazabicyclo[3.1.1]heptan-7-yl]phenyl]benzonitrile',
                                     'reason': 'No sterol backbone found'},
                                 {   'smiles': 'S(CC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)C)[C@H](CC)C)C(O)=O)C',
                                     'name': 'Ala-Ile-Met',
                                     'reason': 'No sterol backbone found'},
                                 {   'smiles': 'O=C([C@H]1C(=C[C@@H](O)[C@@H]2[C@@H]1CC[C@@H](C2)C)C)CCOC(=O)C',
                                     'name': 'Pallidopenilline G',
                                     'reason': 'No sterol backbone found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]1(CC[C@@]2(C)C3=C(CC[C@]12C)[C@@]1(C)CC[C@H](OC(=O)CCCCCCCCCCCCCCC)C(C)(C)[C@]1([H])CC3)[C@H](C)CCC=C(C)C',
                                      'name': 'lanosteryl palmitate',
                                      'reason': 'No sterol backbone found'},
                                  {   'smiles': '[H][C@@]12CC=C3[C@]4([H])CC[C@]([H])([C@H](C)CCC(=C)C(C)C)[C@@]4(C)CC[C@]3([H])[C@@]1(C)CC[C@@H](C2)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                      'name': 'episteryl oleate',
                                      'reason': 'No sterol backbone found'},
                                  {   'smiles': 'C1=2[C@]3([C@](CC[C@@]1([C@@]4([C@](C(C2)=O)(C[C@H]([C@H](C4)O)O)[H])C)[H])([C@@](CC3)([H])[C@@H]([C@@H](CCC(C)(C)O)OP(=O)([O-])[O-])C)C)O',
                                      'name': 'ecdysone 22-phosphate(2-)',
                                      'reason': 'No sterol backbone found'},
                                  {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@H](C(C)=O)[C@@]4(C)CC[C@]3([H])[C@@]1(C)CC[C@@H](C2)OC(=O)CCC(O)=O',
                                      'name': '3beta-hydroxy-5beta-pregnan-20-one '
                                              'hemisuccinate',
                                      'reason': 'No sterol backbone found'},
                                  {   'smiles': '[H][C@@]12CCC3=C(CC[C@]4(C)[C@]([H])(CC[C@@]34[H])[C@H](C)CCC(=C)C(C)C)[C@@]1(C)CC[C@@H](C2)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                      'name': 'fecosteryl oleate',
                                      'reason': 'No sterol backbone found'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2(C)C3=C(CC[C@]12C)[C@@]1(C)CC[C@H](OC(=O)CCCCCCC\\C=C/CCCCCCCC)C(C)(C)[C@]1([H])CC3)[C@H](C)CCC=C(C)C',
                                      'name': 'lanosteryl oleate',
                                      'reason': 'No sterol backbone found'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCC\\C=C/CCCCCCCC)[C@H](C)\\C=C\\[C@H](C)C(C)C',
                                      'name': 'ergosteryl oleate',
                                      'reason': 'No sterol backbone found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])OC(=O)CCC1CCCC1',
                                      'name': 'testosterone cypionate',
                                      'reason': 'No sterol backbone found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])OC(C)=O',
                                      'name': 'testosterone acetate',
                                      'reason': 'No sterol backbone found'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCCCCCCCCCC)[C@H](C)\\C=C\\[C@H](C)C(C)C',
                                      'name': 'ergosteryl palmitate',
                                      'reason': 'No sterol backbone found'}],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as sterol esters, we need to '
                 'identify the presence of a sterol backbone and an ester '
                 'linkage. The sterol backbone typically has a characteristic '
                 'tetracyclic ring system, and the ester linkage is formed by '
                 'the condensation of a carboxylic acid with the 3-hydroxy '
                 'group of the sterol.\n'
                 '\n'
                 'Here’s the Python program that uses RDKit to classify sterol '
                 'esters:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 42,
    'num_false_positives': 80,
    'num_true_negatives': 142157,
    'num_false_negatives': 21,
    'num_negatives': None,
    'precision': 0.3442622950819672,
    'recall': 0.6666666666666666,
    'f1': 0.454054054054054,
    'accuracy': 0.9992902319044272,
    'negative_predictive_value': 0.9998522978238547}