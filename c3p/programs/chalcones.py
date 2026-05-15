"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:23003 chalcone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone that is 1,3-diphenylpropenone (benzylideneacetophenone), ArCH=CH(=O)Ar, and its derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core chalcone structure: Ar-CH=CH-C(=O)-Ar or Ar-CH2-CH2-C(=O)-Ar (dihydrochalcone)
    chalcone_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CH]=[CH]-[C](=O)-[c]2[c][c][c][c][c]2")
    dihydrochalcone_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CH2]-[CH2]-[C](=O)-[c]2[c][c][c][c][c]2")
    
    if not (mol.HasSubstructMatch(chalcone_pattern) or mol.HasSubstructMatch(dihydrochalcone_pattern)):
        return False, "No chalcone core structure found (Ar-CH=CH-C(=O)-Ar or Ar-CH2-CH2-C(=O)-Ar)"

    # Check for the presence of two aromatic rings (Ar)
    aromatic_rings = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"

    # Check for the presence of the propenone chain (CH=CH-C=O or CH2-CH2-C=O)
    propenone_chain = Chem.MolFromSmarts("[CH]=[CH]-[C](=O)")
    dihydro_propenone_chain = Chem.MolFromSmarts("[CH2]-[CH2]-[C](=O)")
    
    if not (mol.HasSubstructMatch(propenone_chain) or mol.HasSubstructMatch(dihydro_propenone_chain)):
        return False, "No propenone chain (CH=CH-C=O or CH2-CH2-C=O) found"

    # Check for the presence of a ketone group (C=O)
    ketone_group = Chem.MolFromSmarts("[C](=O)")
    if not mol.HasSubstructMatch(ketone_group):
        return False, "No ketone group (C=O) found"

    return True, "Contains the core chalcone structure (Ar-CH=CH-C(=O)-Ar or Ar-CH2-CH2-C(=O)-Ar) with possible substitutions on the aromatic rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23086',
                          'name': 'chalcones',
                          'definition': 'A ketone that is '
                                        '1,3-diphenylpropenone '
                                        '(benzylideneacetophenone), '
                                        'ArCH=CH(=O)Ar, and its derivatives '
                                        'formed by substitution.',
                          'parents': [   'CHEBI:51689',
                                         'CHEBI:72544',
                                         'CHEBI:76224'],
                          'xrefs': ['KEGG:C15589'],
                          'all_positive_examples': []},
    'config': None,
    'code_statistics': {   'lines_of_code': 32,
                           'log_lines_of_code': 3.4657359027997265,
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
                                                 0,
                                                 1,
                                                 1,
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
                           'methods_called': [   'MolFromSmarts',
                                                 'HasSubstructMatch',
                                                 'GetIsAromatic',
                                                 'GetAtoms',
                                                 'GetIdx',
                                                 'MolFromSmiles'],
                           'methods_called_count': 6,
                           'smarts_strings': [   '[c]1[c][c][c][c][c]1-[CH]=[CH]-[C](=O)-[c]2[c][c][c][c][c]2',
                                                 '[CH]=[CH]-[C](=O)',
                                                 '[CH2]-[CH2]-[C](=O)',
                                                 '[c]1[c][c][c][c][c]1-[CH2]-[CH2]-[C](=O)-[c]2[c][c][c][c][c]2',
                                                 '[C](=O)'],
                           'smarts_strings_count': 5,
                           'defs': ['is_chalcones(smiles: str):'],
                           'defs_count': 1,
                           'returns': [   'False, "Invalid SMILES string"',
                                          'False, "No chalcone core structure '
                                          'found (Ar-CH=CH-C(=O)-Ar or '
                                          'Ar-CH2-CH2-C(=O)-Ar)"',
                                          'False, "Less than two aromatic '
                                          'rings found"',
                                          'False, "No propenone chain '
                                          '(CH=CH-C=O or CH2-CH2-C=O) found"',
                                          'False, "No ketone group (C=O) '
                                          'found"',
                                          'True, "Contains the core chalcone '
                                          'structure (Ar-CH=CH-C(=O)-Ar or '
                                          'Ar-CH2-CH2-C(=O)-Ar) with possible '
                                          'substitutions on the aromatic '
                                          'rings"'],
                           'returns_count': 6,
                           'complexity': 3.693147180559945},
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'CC(C)=CCC[C@]1(C)Oc2c(C=C1)c(O)c(CC=C(C)C)c(O)c2C(=O)\\C=C\\c1ccc(O)c(O)c1 '
               'NAME: (S)-mallotophilippen E REASON: MISSED No chalcone core '
               'structure found (Ar-CH=CH-C(=O)-Ar)\n'
               ' * SMILES: O(C1=C(C(O)=C(C(=O)CCC2=CC=CC=C2)C(O)=C1)CC=C(C)C)C '
               "NAME: 2',6'-Dihydroxy-4'-methoxy-3'-prenyldihydrochalcone "
               'REASON: MISSED No chalcone core structure found '
               '(Ar-CH=CH-C(=O)-Ar)\n'
               ' * SMILES: [H]C(=C([H])c1ccccc1)C(=O)c1ccccc1 NAME: chalcone '
               'REASON: MISSED No chalcone core structure found '
               '(Ar-CH=CH-C(=O)-Ar)\n'
               ' * SMILES: '
               'O=C(O)C1=C(O)C(C(=O)[C@@H]([C@H](O)C2=C(OC)C=CC(=C2)OC)C)=CC=C1OC '
               'NAME: Chresdihydrochalcone REASON: MISSED No chalcone core '
               'structure found (Ar-CH=CH-C(=O)-Ar)\n'
               ' * SMILES: C1=CC(=C(C=C1C=CC(=O)C2=C(C=C(C=C2)O)O)O)O NAME: '
               '1-(2,4-dihydroxyphenyl)-3-(3,4-dihydroxyphenyl)-2-propen-1-one '
               'REASON: MISSED No chalcone core structure found '
               '(Ar-CH=CH-C(=O)-Ar)\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\Cc1c(O)c2C=CC(C)(C)Oc2c(C(=O)\\C=C\\c2ccc(O)c(O)c2)c1O '
               'NAME: mallotophilippen D REASON: MISSED No chalcone core '
               'structure found (Ar-CH=CH-C(=O)-Ar)\n'
               ' * SMILES: '
               'C1(=C(CC=C(C)C)C(=C(C(=O)/C=C/C2=CC=C(C=C2)O)C(=C1)OC)O)OC '
               "NAME: 4'-O-methylxanthohumol REASON: MISSED No chalcone core "
               'structure found (Ar-CH=CH-C(=O)-Ar)\n'
               ' * SMILES: '
               'O1C(C=CC2=C1C=CC(=C2O)C(=O)/C=C/C3=CC=C(OC)C=C3)(C)C NAME: '
               '4-Methoxylonchocarpin REASON: MISSED No chalcone core '
               'structure found (Ar-CH=CH-C(=O)-Ar)\n'
               ' * SMILES: O(C1=C(C=C(OC)C=C1)/C=C/C(=O)C=2C(OC)=CC=CC2O)C '
               "NAME: 2'-hydroxy-2,5,6'-trimethoxychalcone REASON: MISSED No "
               'chalcone core structure found (Ar-CH=CH-C(=O)-Ar)\n'
               ' * SMILES: S(OC=1C=C(CCC(=O)C2=C(O)C=CC(O)=C2)C=CC1)(O)(=O)=O '
               'NAME: '
               '{3-[3-(2,5-dihydroxyphenyl)-3-oxopropyl]phenyl}oxidanesulfonic '
               'acid REASON: MISSED No chalcone core structure found '
               '(Ar-CH=CH-C(=O)-Ar)\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C1N(C(=O)[C@@H](NC(=O)[C@H](NC(=O)C/C=C/[C@H]2O[C@H]2C(=O)CC(C)C)CN(C)C)C(C)C)[C@@H](CC(C)C)C(=C1C)O',
                                     'name': 'Sandarazol B',
                                     'reason': 'No chalcone core structure '
                                               'found (Ar-CH=CH-C(=O)-Ar or '
                                               'Ar-CH2-CH2-C(=O)-Ar)'},
                                 {   'smiles': 'ClC1=CC=C(C=2OC(COC(C)(C)C(O)=O)=C(N2)C)C=C1',
                                     'name': 'romazarit',
                                     'reason': 'No chalcone core structure '
                                               'found (Ar-CH=CH-C(=O)-Ar or '
                                               'Ar-CH2-CH2-C(=O)-Ar)'},
                                 {   'smiles': 'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C1N=CN=C2N',
                                     'name': '(5Z,8Z,11Z,14Z,17Z)-icosapentaenoyl-CoA(4-)',
                                     'reason': 'No chalcone core structure '
                                               'found (Ar-CH=CH-C(=O)-Ar or '
                                               'Ar-CH2-CH2-C(=O)-Ar)'},
                                 {   'smiles': 'C[C@H]1CN([C@@H](COC2=C(C=C(C=C2)NC(=O)NC3CCCCC3)C(=O)N(C[C@@H]1OC)C)C)C',
                                     'name': '1-cyclohexyl-3-[(4R,7S,8R)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'No chalcone core structure '
                                               'found (Ar-CH=CH-C(=O)-Ar or '
                                               'Ar-CH2-CH2-C(=O)-Ar)'},
                                 {   'smiles': 'C[C@H]1CCCCO[C@@H]([C@H](CN(C(=O)C2=C(O1)C=CC(=C2)NS(=O)(=O)C3=CC=CS3)[C@H](C)CO)C)CN(C)C(=O)NC4=CC=CC5=CC=CC=C54',
                                     'name': '1-[[(3S,9S,10S)-12-[(2R)-1-hydroxypropan-2-yl]-3,10-dimethyl-13-oxo-16-(thiophen-2-ylsulfonylamino)-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-9-yl]methyl]-1-methyl-3-(1-naphthalenyl)urea',
                                     'reason': 'No chalcone core structure '
                                               'found (Ar-CH=CH-C(=O)-Ar or '
                                               'Ar-CH2-CH2-C(=O)-Ar)'},
                                 {   'smiles': 'P(OC[C@@H](COC(CCCCCCCCCCCCCCCC)=O)OC(CCCCCCCCCCCCCCCC)=O)(=O)(OCC[NH3+])[O-]',
                                     'name': '1,2-diheptadecanoyl-sn-glycero-3-phosphoethanolamine '
                                             'zwitterion',
                                     'reason': 'No chalcone core structure '
                                               'found (Ar-CH=CH-C(=O)-Ar or '
                                               'Ar-CH2-CH2-C(=O)-Ar)'},
                                 {   'smiles': 'O=C1N(C(CC1)C(O)=O)CCC(N)C(O)=O',
                                     'name': "(2S,3'S)-alpha-Amino-2-carboxy-5-oxo-1-pyrrolidinebutanoic "
                                             'acid',
                                     'reason': 'No chalcone core structure '
                                               'found (Ar-CH=CH-C(=O)-Ar or '
                                               'Ar-CH2-CH2-C(=O)-Ar)'},
                                 {   'smiles': 'CC(=O)[C@@H](O)[C@H](O)CC(=O)C([O-])=O',
                                     'name': '3,7-dideoxy-D-threo-hepto-2,6-diuolosonate',
                                     'reason': 'No chalcone core structure '
                                               'found (Ar-CH=CH-C(=O)-Ar or '
                                               'Ar-CH2-CH2-C(=O)-Ar)'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO)[C@H](O[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)[C@@H](NC(=O)C)[C@@H]1OC[C@H]4O[C@H](O)[C@H](NC(=O)C)[C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5O)CO)[C@H]4O)CO',
                                     'name': 'N-[(2S,3R,4R,5R,6R)-6-[[(2R,3R,4R,5S,6R)-3-Acetamido-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxymethyl]-2,5-dihydroxy-4-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-3-yl]acetamide',
                                     'reason': 'No chalcone core structure '
                                               'found (Ar-CH=CH-C(=O)-Ar or '
                                               'Ar-CH2-CH2-C(=O)-Ar)'},
                                 {   'smiles': 'COC(=O)C1=CC=CC=C1NC(=O)NC2=NC=C(C=C2)Cl',
                                     'name': '2-[[[(5-chloro-2-pyridinyl)amino]-oxomethyl]amino]benzoic '
                                             'acid methyl ester',
                                     'reason': 'No chalcone core structure '
                                               'found (Ar-CH=CH-C(=O)-Ar or '
                                               'Ar-CH2-CH2-C(=O)-Ar)'}],
    'sample_false_negatives': [   {   'smiles': 'O=C(O)C1=C(O)C(C(=O)[C@@H]([C@H](O)C2=C(OC)C=CC(=C2)OC)C)=CC=C1OC',
                                      'name': 'Chresdihydrochalcone',
                                      'reason': 'No chalcone core structure '
                                                'found (Ar-CH=CH-C(=O)-Ar or '
                                                'Ar-CH2-CH2-C(=O)-Ar)'},
                                  {   'smiles': 'OC(CC1=CC=C(O)C=C1)C(=O)C2=C(O)C=C(OC\\C=C(\\CC=C(C)C)/C)C=C2',
                                      'name': "alpha,4,2'-Trihydroxy-4-O-geranyldihydrochalcone",
                                      'reason': 'No chalcone core structure '
                                                'found (Ar-CH=CH-C(=O)-Ar or '
                                                'Ar-CH2-CH2-C(=O)-Ar)'},
                                  {   'smiles': 'CO\\C(=C/C(=O)c1c(OC)cc2OC(C)(C)C=Cc2c1OC)c1ccccc1',
                                      'name': '7-methoxypraecansone B',
                                      'reason': 'No chalcone core structure '
                                                'found (Ar-CH=CH-C(=O)-Ar or '
                                                'Ar-CH2-CH2-C(=O)-Ar)'},
                                  {   'smiles': 'O(C1=C(C(O)=C(C(OC)=C1)C)C(=O)/C=C(\\O)/C2=CC=CC=C2)C',
                                      'name': 'beta,2-Dihydroxy-4,6-dimethoxy-3-methylchalcone',
                                      'reason': 'No chalcone core structure '
                                                'found (Ar-CH=CH-C(=O)-Ar or '
                                                'Ar-CH2-CH2-C(=O)-Ar)'},
                                  {   'smiles': 'O(C=1C(CC=C(C)C)=CC(=C(OC)C1)C(=O)/C=C(\\O)/C2=CC=CC=C2)C',
                                      'name': 'Pongagallone A',
                                      'reason': 'No chalcone core structure '
                                                'found (Ar-CH=CH-C(=O)-Ar or '
                                                'Ar-CH2-CH2-C(=O)-Ar)'},
                                  {   'smiles': 'S(OC(CC(=O)C1=C(O)C=CC(O)=C1)C2=CC=CC=C2)(O)(=O)=O',
                                      'name': '[3-(2,5-dihydroxyphenyl)-3-oxo-1-phenylpropoxy]sulfonic '
                                              'acid',
                                      'reason': 'No chalcone core structure '
                                                'found (Ar-CH=CH-C(=O)-Ar or '
                                                'Ar-CH2-CH2-C(=O)-Ar)'},
                                  {   'smiles': 'OC(CC1=CC(CC=C(C)C)=C(O)C=C1)C(=O)C=2C=C(CC=C(C)C)C(O)=CC2O',
                                      'name': '(R)-Kanzonol Y',
                                      'reason': 'No chalcone core structure '
                                                'found (Ar-CH=CH-C(=O)-Ar or '
                                                'Ar-CH2-CH2-C(=O)-Ar)'},
                                  {   'smiles': 'COC1=C(C(=CC(=C1)C=C(C#N)C(=O)C2=CC(=C(C=C2)O)O)I)O',
                                      'name': '2-[(3,4-dihydroxyphenyl)-oxomethyl]-3-(4-hydroxy-3-iodo-5-methoxyphenyl)-2-propenenitrile',
                                      'reason': 'No chalcone core structure '
                                                'found (Ar-CH=CH-C(=O)-Ar or '
                                                'Ar-CH2-CH2-C(=O)-Ar)'},
                                  {   'smiles': 'S(OC=1C(OC)=CC(CC(=O)C(=O)C2=C(O)C(OC)=C(O)C=C2O)=CC1O)(O)(=O)=O',
                                      'name': '{4-[2,3-dioxo-3-(2,4,6-trihydroxy-3-methoxyphenyl)propyl]-2-hydroxy-6-methoxyphenyl}oxidanesulfonic '
                                              'acid',
                                      'reason': 'No chalcone core structure '
                                                'found (Ar-CH=CH-C(=O)-Ar or '
                                                'Ar-CH2-CH2-C(=O)-Ar)'},
                                  {   'smiles': 'OC(CC1=CC(CC=C(C)C)=C(O)C(OC)=C1)C(=O)C=2C=C(CC=C(C)C)C(O)=CC2O',
                                      'name': '1-[2,4-dihydroxy-5-(3-methylbut-2-en-1-yl)phenyl]-2-hydroxy-3-[4-hydroxy-3-methoxy-5-(3-methylbut-2-en-1-yl)phenyl]propan-1-one',
                                      'reason': 'No chalcone core structure '
                                                'found (Ar-CH=CH-C(=O)-Ar or '
                                                'Ar-CH2-CH2-C(=O)-Ar)'}],
    'attempt': 1,
    'reasoning': 'The previous code failed to correctly identify chalcones '
                 'because the SMARTS pattern used to detect the core chalcone '
                 'structure was too restrictive. Specifically, the pattern '
                 '`[c;H1]:[c;H1]-[CH]=[CH]-[C](=O)-[c;H1]:[c;H1]` requires '
                 'that the aromatic rings have exactly one hydrogen atom, '
                 'which is not always the case due to substitutions. '
                 'Additionally, the pattern does not account for the '
                 'possibility of dihydrochalcones, where the double bond in '
                 'the propenone chain is reduced to a single bond.\n'
                 '\n'
                 'To improve the code, we need to:\n'
                 '1. Use a more flexible SMARTS pattern that allows for '
                 'substitutions on the aromatic rings.\n'
                 '2. Account for both chalcones and dihydrochalcones by '
                 'allowing for either a double or single bond in the propenone '
                 'chain.\n'
                 '3. Ensure that the pattern matches the core structure of '
                 'chalcones, which consists of two aromatic rings connected by '
                 'a propenone chain.\n'
                 '\n'
                 'Here is the revised code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 205,
    'num_false_positives': 100,
    'num_true_negatives': 96977,
    'num_false_negatives': 18,
    'num_negatives': None,
    'precision': 0.6721311475409836,
    'recall': 0.9192825112107623,
    'f1': 0.7765151515151516,
    'accuracy': 0.9987872559095581,
    'negative_predictive_value': 0.9998144234238878}