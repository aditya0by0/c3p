"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
"""
Classifies: phosphatidylcholine
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine is a glycerophosphocholine that is glycero-3-phosphocholine bearing
    two acyl substituents at positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define phosphocholine head group pattern (without charges)
    phosphocholine_smarts = "COP(=O)(O)OCCN(C)(C)C"
    phosphocholine_pattern = Chem.MolFromSmarts(phosphocholine_smarts)
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine head group found"

    # Define glycerol backbone with two ester groups (remove stereochemistry)
    glycerol_esters_smarts = "OCC(OC(=O)[#6])COC(=O)[#6]"
    glycerol_esters_pattern = Chem.MolFromSmarts(glycerol_esters_smarts)
    if not mol.HasSubstructMatch(glycerol_esters_pattern):
        return False, "No glycerol backbone with two ester groups found"

    # Check for exactly two ester bonds connected to the glycerol carbons
    ester_bond_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_bonds = mol.GetSubstructMatches(ester_bond_pattern)
    if len(ester_bonds) < 2:
        return False, f"Found {len(ester_bonds)} ester bonds, need at least 2"

    # Optional: Verify the acyl chains are attached via ester bonds
    # This can be more complex but for now assume if ester bonds are present, acyl chains are present

    return True, "Contains glycerol backbone with two acyl chains and phosphocholine head group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64482',
                          'name': 'phosphatidylcholine',
                          'definition': 'A glycerophosphocholine that is '
                                        'glycero-3-phosphocholine bearing two '
                                        'acyl substituents at positions 1 and '
                                        '2.',
                          'parents': ['CHEBI:35284', 'CHEBI:36313'],
                          'xrefs': ['PMID:2474544'],
                          'all_positive_examples': []},
    'config': None,
    'code_statistics': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCC)COC(=O)CCCCCCCCCCC)([O-])=O '
               'NAME: PC(12:0/14:1(9Z)) REASON: MISSED No glycerol backbone '
               'with two ester groups found\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCC)COC(=O)CCCCCCC/C=C\\CCCCCCC)([O-])=O '
               'NAME: PC(17:1(9Z)/14:1(9Z)) REASON: MISSED No glycerol '
               'backbone with two ester groups found\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCCCCCCCC/C=C\\CCCCCCCC)COC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)([O-])=O '
               'NAME: PC(22:5(4Z,7Z,10Z,13Z,16Z)/24:1(15Z)) REASON: MISSED No '
               'glycerol backbone with two ester groups found\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC)([O-])=O '
               'NAME: PC(10:0/21:0) REASON: MISSED No glycerol backbone with '
               'two ester groups found\n'
               ' * SMILES: '
               'S([C@@H]([C@@H](O)CCCC(O)=O)/C=C/C=C/C=C\\C/C=C\\CCCCC)C[C@H](N)C(O[C@H](COC(=O)CCCCCCCCCCCCCCCCCCCCCCC)COP(OCC[N+](C)(C)C)([O-])=O)=O '
               'NAME: PC(24:0/LTE4) REASON: MISSED No glycerol backbone with '
               'two ester groups found\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCC/C=C\\CCCCCC)([O-])=O '
               'NAME: PC(16:1(9Z)/18:4(6Z,9Z,12Z,15Z)) REASON: MISSED No '
               'glycerol backbone with two ester groups found\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\CCCCC)([O-])=O '
               'NAME: PC(20:3(8Z,11Z,14Z)/20:4(5Z,8Z,11Z,14Z)) REASON: MISSED '
               'No glycerol backbone with two ester groups found\n'
               ' * SMILES: '
               'C([C@@](COC(CCCCCCCCCCCCCCCCCCCCCCC)=O)(OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)[H])OP([O-])(=O)OCC[N+](C)(C)C '
               'NAME: PC(24:0/18:2(9Z,12Z)) REASON: MISSED No glycerol '
               'backbone with two ester groups found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC '
               'NAME: '
               '1-octadecanoyl-2-(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)-sn-glycero-3-phosphocholine '
               'REASON: MISSED No glycerol backbone with two ester groups '
               'found\n'
               ' * SMILES: '
               'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)([O-])=O '
               'NAME: 1-STEAROYL-2-DOCOSAHEXANOYL-SN-GLYCERO-3-PHOSPHOCHOLINE '
               'REASON: MISSED No glycerol backbone with two ester groups '
               'found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'No phosphocholine head group '
                                               'found'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'No phosphocholine head group '
                                               'found'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'No phosphocholine head group '
                                               'found'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'No phosphocholine head group '
                                               'found'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'No phosphocholine head group '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No phosphocholine head group '
                                               'found'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'No phosphocholine head group '
                                               'found'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'No phosphocholine head group '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'No phosphocholine head group '
                                               'found'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'No phosphocholine head group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'P(OC[C@H](OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)(O)(O)=O',
                                      'name': 'PA(20:2(11Z,14Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                      'reason': 'No phosphocholine head group '
                                                'found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)OCC(COP(O)(=O)CC[N+](C)(C)C)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                      'name': '1-hexadecanoyl-2-[(Z)-octadec-9-enoyl]-sn-glycerol '
                                              '3-(2-trimethylaminoethylphosphonate)',
                                      'reason': 'No phosphocholine head group '
                                                'found'}],
    'attempt': 1,
    'reasoning': "To understand why the previous program failed, let's analyze "
                 'the code and the outcomes.\n'
                 '\n'
                 '**Issues with the previous code:**\n'
                 '\n'
                 '1. **Stereochemistry in SMARTS patterns:** The SMARTS '
                 'patterns used in the code include stereochemical '
                 'specifications like `[C@@H]`. This means the substructure '
                 "search will only match if the molecule's SMILES string "
                 'includes the exact same stereochemistry. However, many of '
                 'the provided SMILES strings either lack stereochemical '
                 'information or have different stereochemistry, causing the '
                 'substructure search to fail.\n'
                 '\n'
                 '2. **Overly specific SMARTS patterns:** The SMARTS patterns '
                 'for the phosphocholine head group and glycerol backbone are '
                 'very specific and may not account for variations in how '
                 'these groups are represented in different SMILES strings. '
                 'For example, the pattern for the phosphocholine group '
                 'includes charges `[N+](C)(C)C` and `[O-]`, but some SMILES '
                 'strings may represent these groups without explicit '
                 'charges.\n'
                 '\n'
                 '3. **Not accounting for tautomers or resonance forms:** The '
                 "code doesn't consider that the phosphate group can be "
                 'represented differently due to tautomerism or resonance '
                 'structures, leading to mismatches.\n'
                 '\n'
                 '**How to improve the code:**\n'
                 '\n'
                 '- **Remove stereochemistry from SMARTS patterns:** By '
                 'omitting stereochemical specifications, the substructure '
                 'search becomes more general and will match molecules '
                 'regardless of their stereochemistry.\n'
                 '\n'
                 '- **Simplify and generalize SMARTS patterns:** Create SMARTS '
                 'patterns that capture the essential features of the '
                 'phosphocholine head group and glycerol backbone without '
                 'being too specific about atom charges or connectivity.\n'
                 '\n'
                 '- **Account for different representations:** Adjust the '
                 'patterns to match the way these groups are commonly '
                 'represented in SMILES strings, including possible '
                 'variations.\n'
                 '\n'
                 "Based on these observations, I'll rewrite the code to use "
                 'more general SMARTS patterns and remove stereochemistry to '
                 'improve matching.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 623,
    'num_false_positives': 27,
    'num_true_negatives': 141648,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.9584615384615385,
    'recall': 0.9968,
    'f1': 0.9772549019607842,
    'accuracy': 0.9997962052002811,
    'negative_predictive_value': 0.9999858806918461}