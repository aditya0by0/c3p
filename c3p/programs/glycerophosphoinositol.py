"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: CHEBI:28874 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is a glycerophospholipid with inositol esterified to the phosphate group at the sn-3 position of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for phosphate group attached to the sn-3 position
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group found"

    # Look for inositol group (six-membered ring with multiple hydroxyls)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if len(inositol_matches) == 0:
        return False, "No inositol group found"

    # Check if the inositol is esterified to the phosphate group
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1OP(=O)(O)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "Inositol not esterified to phosphate group"

    # Check for fatty acid chains (long carbon chains attached to glycerol)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Missing fatty acid chains"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - glycerophosphoinositols typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycerophosphoinositol"

    return True, "Contains glycerol backbone with inositol esterified to phosphate and fatty acid chains"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36315',
                          'name': 'glycerophosphoinositol',
                          'definition': 'Any glycerophospholipid having the '
                                        'polar alcohol inositol esterified to '
                                        'the phosphate group at the sn-3 '
                                        'position of the glycerol backbone.',
                          'parents': ['CHEBI:37739'],
                          'xrefs': ['PMID:14706866'],
                          'all_positive_examples': []},
    'config': None,
    'code_statistics': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O(C1C(O)C(OC(OC2=C(OC=3C(C2=O)=C(O)C=C(O)C3)C4=CC(O)=C(O)C=C4)C1O)CO)C5OC(C(O)C(O)C5O)CO',
                                     'name': 'Quercetin 3-beta-laminaribioside',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C(N[C@@H](CC1=CC=C(O)C=C1)C(O)=O)[C@H]2N(CCC2)C(=O)[C@@H](N)CC3=CC=CC=C3',
                                     'name': 'Phe-Pro-Tyr',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O[C@@H]1[C@@H](C(=C2C[C@](CO)(C)C[C@H]2[C@H](C1)C)CO)CO',
                                     'name': '4alpha,11,12,14-tetrahydroxy-1-tremulene',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O(CCCCCCCCCC\\C=C/C=C\\CC)C(=O)C',
                                     'name': '11Z,13Z-Hexadecadienyl acetate',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C(NC(C(=O)NC(CN(CCO)C)C)(C)C)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C1N(C(=O)C(CCCCCCCC)C)CCC1)CC(CC(O)CC(=O)CC)C)C)(C)C)C(CC)C)C(C)C)(C)C',
                                     'name': 'Roseoferin A1',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'O=C1C(=C(CC2=CC=C(O)C=C2)C(C1(C3=CC=C(O)C=C3)C4(C5=CC=C(O)C=C5)C(=O)C(CC6=CC=C(O)C=C6)=C(C4=O)CC7=CC=C(O)C=C7)=O)CC8=CC=C(O)C=C8',
                                     'name': 'Nostotrebin 6',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'ClC1=C(O)C=CC(=C1)NC2=NC=NC=3C2=CC(OCCCNCCCO)=C(C3)OC',
                                     'name': '2-chloro-4-[(6-{3-[(3-hydroxypropyl)amino]propoxy}-7-methoxyquinazolin-4-yl)amino]phenol',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C1[C@@H]2C([C@@H](N2)CN1CC3=COC=N3)C4=CC=C(C=C4)C5=CC=C(C=C5)C#N',
                                     'name': '4-[4-[(1S,5R)-3-(4-oxazolylmethyl)-3,6-diazabicyclo[3.1.1]heptan-7-yl]phenyl]benzonitrile',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'S(CC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)C)[C@H](CC)C)C(O)=O)C',
                                     'name': 'Ala-Ile-Met',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C([C@H]1C(=C[C@@H](O)[C@@H]2[C@@H]1CC[C@@H](C2)C)C)CCOC(=O)C',
                                     'name': 'Pallidopenilline G',
                                     'reason': 'No phosphate group found'}],
    'sample_false_negatives': [   {   'smiles': 'CCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCC',
                                      'name': "1,2-dibutyryl-sn-glycero-3-phospho-(1'D-myo-inositol)",
                                      'reason': 'Molecular weight too low for '
                                                'glycerophosphoinositol'},
                                  {   'smiles': '[C@@H]1([C@@H]([C@@H](O)[C@H]([C@H]([C@@H]1O)OP(=O)(OCOC(=O)CCC)OCOC(=O)CCC)O)OP(=O)(OCOC(=O)CCC)OCOC(=O)CCC)OP(=O)(OCOC(=O)CCC)OCOC(=O)CCC',
                                      'name': 'D-myo-Ins(1,4,5)P3 '
                                              'hexakis(butyryloxymethyl) ester',
                                      'reason': 'No glycerol backbone found'},
                                  {   'smiles': 'P(OC1[C@H](O)[C@@H](O)C(O)[C@@H](O)[C@H]1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)([O-])=O',
                                      'name': '18:1-18:3-PI',
                                      'reason': 'No phosphate group found'},
                                  {   'smiles': 'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCC)COC(=O)CCCCCCCCCCCCCCCC)([O-])=O.[NH4+]',
                                      'name': "1-heptadecanoyl-2-(9Z-tetradecenoyl)-sn-glycero-3-phospho-(1'-myo-inositol)(ammonium "
                                              'salt)',
                                      'reason': 'No phosphate group found'},
                                  {   'smiles': 'P(OC1[C@H](O)[C@@H](O)C(O)[C@@H](O)[C@H]1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCCCCCCCCCC)([O-])=O',
                                      'name': '18:0-18:3-PI',
                                      'reason': 'No phosphate group found'},
                                  {   'smiles': 'P(OC1[C@H](O)[C@H](O)C(O)[C@H](O)[C@H]1O)(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)([O-])=O',
                                      'name': '18:2-18:2-PI',
                                      'reason': 'No phosphate group found'}],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as glycerophosphoinositols, we '
                 'need to identify specific structural features in the '
                 'molecule. A glycerophosphoinositol is a glycerophospholipid '
                 'with inositol esterified to the phosphate group at the sn-3 '
                 'position of the glycerol backbone. The key features to look '
                 'for are:\n'
                 '\n'
                 '1. **Glycerol backbone**: A three-carbon chain with hydroxyl '
                 'groups.\n'
                 '2. **Phosphate group**: Attached to the sn-3 position of the '
                 'glycerol backbone.\n'
                 '3. **Inositol group**: A six-membered ring with multiple '
                 'hydroxyl groups, esterified to the phosphate group.\n'
                 '4. **Fatty acid chains**: Typically attached to the sn-1 and '
                 'sn-2 positions of the glycerol backbone.\n'
                 '\n'
                 'The function `is_glycerophosphoinositol` will check for '
                 "these features using RDKit's substructure matching "
                 'capabilities.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 418,
    'num_false_positives': 8,
    'num_true_negatives': 141868,
    'num_false_negatives': 6,
    'num_negatives': None,
    'precision': 0.9812206572769953,
    'recall': 0.9858490566037735,
    'f1': 0.9835294117647059,
    'accuracy': 0.999901616303584,
    'negative_predictive_value': 0.9999577089530146}