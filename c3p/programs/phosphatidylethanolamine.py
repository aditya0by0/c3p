"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine is a glycerophospholipid with a phosphatidyl group esterified to ethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with two ester groups and one phosphate group)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Look for phosphate group (P(=O)(O))
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for ethanolamine group (OCCN)
    ethanolamine_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][NX3]")
    if not mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "No ethanolamine group found"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - phosphatidylethanolamines typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidylethanolamine"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for phosphatidylethanolamine"
    if o_count < 6:
        return False, "Too few oxygens for phosphatidylethanolamine"

    return True, "Contains glycerol backbone with 2 fatty acid chains, a phosphate group, and an ethanolamine group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16038',
                          'name': 'phosphatidylethanolamine',
                          'definition': 'A class of glycerophospholipids in '
                                        'which a phosphatidyl group is '
                                        'esterified to the hydroxy group of '
                                        'ethanolamine.',
                          'parents': ['CHEBI:36314'],
                          'xrefs': [   'DrugBank:DB04327',
                                       'HMDB:HMDB0060501',
                                       'KEGG:C00350',
                                       'LIPID_MAPS_instance:LMGP02010000',
                                       'PMID:10540156',
                                       'PMID:11042504',
                                       'PMID:11159918',
                                       'PMID:11829744',
                                       'PMID:12139474',
                                       'PMID:15653902',
                                       'PMID:16037249',
                                       'PMID:16303767',
                                       'PMID:16620109',
                                       'PMID:18034796',
                                       'PMID:18259190',
                                       'PMID:18398168',
                                       'PMID:18462396',
                                       'PMID:18570887',
                                       'PMID:18957134',
                                       'PMID:19393163',
                                       'PMID:23354482',
                                       'PMID:23369752',
                                       'PMID:23543734',
                                       'PMID:3196084',
                                       'PMID:7980848',
                                       'Wikipedia:Phosphatidylethanolamine'],
                          'all_positive_examples': []},
    'config': None,
    'code_statistics': {   'lines_of_code': 51,
                           'log_lines_of_code': 3.9318256327243257,
                           'indent_by_line': [   1,
                                                 1,
                                                 1,
                                                 0,
                                                 1,
                                                 2,
                                                 0,
                                                 1,
                                                 2,
                                                 2,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 0,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 2,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 0,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 0,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 0,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 0,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 0,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 0,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 1,
                                                 2,
                                                 0,
                                                 1,
                                                 0,
                                                 0],
                           'max_indent': 2,
                           'imports': [   'from rdkit import Chem',
                                          'from rdkit.Chem import AllChem',
                                          'from rdkit.Chem import '
                                          'rdMolDescriptors'],
                           'imports_count': 3,
                           'methods_called': [   'CalcExactMolWt',
                                                 'GetAtoms',
                                                 'GetSubstructMatches',
                                                 'HasSubstructMatch',
                                                 'GetAtomicNum',
                                                 'MolFromSmiles',
                                                 'MolFromSmarts',
                                                 'CalcNumRotatableBonds'],
                           'methods_called_count': 8,
                           'smarts_strings': [   '[PX4](=[OX1])([OX2])[OX2]',
                                                 '[OX2][CX3](=[OX1])',
                                                 '[OX2][CX4][CX4][NX3]',
                                                 '[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]',
                                                 '[CH2X4][CHX4][CH2X4]'],
                           'smarts_strings_count': 5,
                           'defs': [   'is_phosphatidylethanolamine(smiles: '
                                       'str):'],
                           'defs_count': 1,
                           'returns': [   'False, "Invalid SMILES string"',
                                          'False, "No glycerol backbone found"',
                                          'False, f"Found {len(ester_matches)} '
                                          'ester groups, need exactly 2"',
                                          'False, "No phosphate group found"',
                                          'False, "No ethanolamine group '
                                          'found"',
                                          'False, f"Missing fatty acid chains, '
                                          'got {len(fatty_acid_matches)}"',
                                          'False, "Chains too short to be '
                                          'fatty acids"',
                                          'False, "Molecular weight too low '
                                          'for phosphatidylethanolamine"',
                                          'False, "Too few carbons for '
                                          'phosphatidylethanolamine"',
                                          'False, "Too few oxygens for '
                                          'phosphatidylethanolamine"',
                                          'True, "Contains glycerol backbone '
                                          'with 2 fatty acid chains, a '
                                          'phosphate group, and an '
                                          'ethanolamine group"'],
                           'returns_count': 11,
                           'complexity': 5.186365126544866},
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C1COCCN1C2=C(C=C(C=C2)C(F)(F)F)NC(=O)C3=CC(=CN=C3)Br',
                                     'name': '5-bromo-N-[2-(4-morpholinyl)-5-(trifluoromethyl)phenyl]-3-pyridinecarboxamide',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C[C@H]1C[C@@H]2CC[C@@H](O2)[C@@H](C)C(=O)O[C@H](C)C[C@H]2CC[C@H](O2)[C@H](C)C(=O)O[C@@H](C)C[C@@H]2CC[C@@H](O2)[C@@H](C)C(=O)O[C@H](C)C[C@H]2CC[C@H](O2)[C@H](C)C(=O)O1',
                                     'name': 'Nonactin',
                                     'reason': 'Found 4 ester groups, need '
                                               'exactly 2'},
                                 {   'smiles': 'O1OC23C(C14C(C5C(C(CC5)C(CCC(C(C)C)C)C)(CC4)C)=CC2=O)(CCC(O)C3)C',
                                     'name': '5,9-Epidioxy-3-hydroxyergost-7-en-6-one',
                                     'reason': 'Found 0 ester groups, need '
                                               'exactly 2'},
                                 {   'smiles': 'CC1=C(SC=C1)C(=O)N2CCCC(C2)CNS(=O)(=O)C3=CC=C(C=C3)OC',
                                     'name': '4-methoxy-N-[[1-[(3-methyl-2-thiophenyl)-oxomethyl]-3-piperidinyl]methyl]benzenesulfonamide',
                                     'reason': 'Found 0 ester groups, need '
                                               'exactly 2'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCCCCCC)C[C@@H](O)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC',
                                     'name': 'DG(16:0/0:0/20:1n9)',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'Oc1cc2CC3(O)COc4c(O)c(O)ccc4C3c2cc1O',
                                     'name': 'haematoxylin',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O([C@@H]1[C@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO)O)O[C@@H]([C@H]1O)CO[C@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O)CO)[C@@H]4O)CO[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O)CO)CO)[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO',
                                     'name': 'N-[(3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4S,5R,6R)-4-[(2R,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5R,6R)-6-[[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O([C@@H]1[C@H](O)[C@@H](O)[C@@H](O[C@@H]1CO)O)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO',
                                     'name': 'beta-D-Glcp-(1->4)-beta-D-Galp',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C(N)C(/C=C/[N+]([O-])=NC(C(O)C)C)CCC',
                                     'name': 'Maniwamycin F',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O1C2=C(C(=O)C(C3=C(O)C=C(O)C=C3)=C1)C(OC)=CC(O)=C2',
                                     'name': 'Barpisoflavone A',
                                     'reason': 'No glycerol backbone found'}],
    'sample_false_negatives': [   {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCCCCCCCCCCCCCC)(OCC[NH3+])([O-])=O',
                                      'name': '24:0-18:2-PE',
                                      'reason': 'No ethanolamine group found'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCC(O)=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                      'name': 'OA-PE',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCC(=O)/C=C/C(O)=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                      'name': 'OKHdiA-PE',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCC(O)=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                      'name': 'OG-PE',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCC)COC(=O)CCCCCCC)(OCCN)(O)=O',
                                      'name': 'PE(8:0/8:0)',
                                      'reason': 'Molecular weight too low for '
                                                'phosphatidylethanolamine'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCC(=O)/C=C/C(O)=O)COC(=O)CCCCCCCCCCCCCCC)(OCCN)(O)=O',
                                      'name': 'PKODiA-PE',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'C(CN)OP(OC[C@@H](COC(CCCCC)=O)OC(=O)CCCCC)(=O)O',
                                      'name': '1,2-dihexanoyl-sn-glycero-3-phosphoethanolamine',
                                      'reason': 'Molecular weight too low for '
                                                'phosphatidylethanolamine'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCC=O)COC(=O)CCCCCCCCCCCCCCC)(OCCN)(O)=O',
                                      'name': 'POB-PE',
                                      'reason': 'Molecular weight too low for '
                                                'phosphatidylethanolamine'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCC(O)=O)COC(=O)CCCCCCCCCCCCCCC)(OCCN)(O)=O',
                                      'name': 'PG-PE',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OC[C@H](OC(=O)CCCC(=O)/C=C/C(O)=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                      'name': 'OKODiA-PE',
                                      'reason': 'Found 3 ester groups, need '
                                                'exactly 2'}],
    'attempt': 0,
    'reasoning': 'To classify phosphatidylethanolamines, we need to identify '
                 'the key structural features of this class of molecules. '
                 'Phosphatidylethanolamines are a type of glycerophospholipid '
                 'characterized by a glycerol backbone, two fatty acid chains '
                 'attached via ester bonds, and a phosphate group esterified '
                 'to ethanolamine. The general structure can be represented '
                 'as:\n'
                 '\n'
                 '- Glycerol backbone: `C-C-C` with two ester-linked fatty '
                 'acids and one phosphate group.\n'
                 '- Phosphate group: `P(=O)(O)`\n'
                 '- Ethanolamine group: `OCCN`\n'
                 '\n'
                 'The program will check for these features in the SMILES '
                 "string. Here's the Python code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 792,
    'num_false_positives': 96,
    'num_true_negatives': 141384,
    'num_false_negatives': 28,
    'num_negatives': None,
    'precision': 0.8918918918918919,
    'recall': 0.9658536585365853,
    'f1': 0.927400468384075,
    'accuracy': 0.9991286015460296,
    'negative_predictive_value': 0.9998019970016688}