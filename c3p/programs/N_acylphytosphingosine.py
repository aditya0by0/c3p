"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:35801 N-acylphytosphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phytosphingosine backbone pattern
    phytosphingosine_pattern = Chem.MolFromSmarts("[C@H](CO)[C@@H](O)[C@@H](O)CCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"
    
    # Look for acyl group attached to nitrogen
    acyl_pattern = Chem.MolFromSmarts("[N;H1][C;X3](=[O;X1])[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group attached to nitrogen found"
    
    # Look for galactose or glucose moiety
    sugar_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)O)O)O)O")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No galactose or glucose moiety found"
    
    return True, "Contains phytosphingosine backbone with acyl group attached to nitrogen and galactose/glucose moiety"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:31998',
                          'name': 'N-acylphytosphingosine',
                          'definition': 'A ceramide that is phytosphingosine '
                                        'having a fatty acyl group attached to '
                                        'the nitrogen.',
                          'parents': ['CHEBI:139051'],
                          'xrefs': [   'KEGG:C12145',
                                       'LIPID_MAPS_instance:LMSP01030000'],
                          'all_positive_examples': []},
    'config': None,
    'code_statistics': {   'lines_of_code': 26,
                           'log_lines_of_code': 3.258096538021482,
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
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 1,
                                                 1,
                                                 1,
                                                 1,
                                                 2,
                                                 1,
                                                 1,
                                                 0],
                           'max_indent': 2,
                           'imports': [   'from rdkit import Chem',
                                          'from rdkit.Chem import AllChem'],
                           'imports_count': 2,
                           'methods_called': [   'MolFromSmiles',
                                                 'MolFromSmarts',
                                                 'HasSubstructMatch'],
                           'methods_called_count': 3,
                           'smarts_strings': [   '[C@H](CO)[C@@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                                 '[N;H1][C;X3](=[O;X1])[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]~[C;X4]',
                                                 '[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)O)O)O)O'],
                           'smarts_strings_count': 3,
                           'defs': ['is_N_acylphytosphingosine(smiles: str):'],
                           'defs_count': 1,
                           'returns': [   'False, "Invalid SMILES string"',
                                          'False, "No phytosphingosine '
                                          'backbone found"',
                                          'False, "No acyl group attached to '
                                          'nitrogen found"',
                                          'False, "No galactose or glucose '
                                          'moiety found"',
                                          'True, "Contains phytosphingosine '
                                          'backbone with acyl group attached '
                                          'to nitrogen and galactose/glucose '
                                          'moiety"'],
                           'returns_count': 5,
                           'complexity': 2.851619307604296},
    'message': '\n'
               'Error: unterminated string literal (detected at line 19) '
               '(<string>, line 19)\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n'
               'IMPORTANT NOTE: I do not have 100% confidence in the benchmark '
               'I am using. There may be occasional and systematic mistakes. '
               'Use your best judgment, and if you think the classifications '
               'your program are consistent with your understanding if the '
               'meaning of the chemical class, then you can ignore outliers, '
               'but explain your reasoning in doing so. I have great '
               'confidence in your broad understanding of chemistry and your '
               'ability to translate this into code.',
    'sample_true_negatives': [   {   'smiles': 'C1CCC(CC1)C#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CNCC(=O)N4[C@H]3CO',
                                     'name': '(6R,7R,8R)-7-[4-(2-cyclohexylethynyl)phenyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one',
                                     'reason': 'No phytosphingosine backbone '
                                               'found'},
                                 {   'smiles': 'CCCNC(=O)[C@@H]1[C@H]([C@@H]2CN3C(=CC=C(C3=O)C=CC)[C@H]1N2CC4CCCC4)CO',
                                     'name': 'LSM-14611',
                                     'reason': 'No phytosphingosine backbone '
                                               'found'},
                                 {   'smiles': 'c1ccc(nc1)-c1cccc(n1)-c1ccccn1',
                                     'name': "2,2':6',2''-terpyridine",
                                     'reason': 'No phytosphingosine backbone '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=CC=CC=C2C3=CC=CC=C3CO[C@@H]1CN(C)C(=O)CC4=CC=CC=C4)[C@H](C)CO',
                                     'name': 'LSM-12861',
                                     'reason': 'No phytosphingosine backbone '
                                               'found'},
                                 {   'smiles': 'C1(=CC(=CC(=C1C(=O)CCC2=CC(=C(C=C2)O)O)O)O)O[C@H]3[C@@H]([C@H]([C@@H]([C@@H](CO[C@H]4[C@@H]([C@H]([C@@H](CO4)O)O)O)O3)O)O)O',
                                     'name': '3-hydroxyphloretin '
                                             "2'-O-xylosylglucoside",
                                     'reason': 'No phytosphingosine backbone '
                                               'found'},
                                 {   'smiles': 'C1[C@]2([C@@]3(CC[C@H]4C(=C[C@]3(CC[C@@]2(C(C[C@@H]1O)(C)C)[H])C4)C)[H])C',
                                     'name': '2alpha-hydroxy-ent-isokaurene',
                                     'reason': 'No phytosphingosine backbone '
                                               'found'},
                                 {   'smiles': 'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]3O)CO)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]8[C@H](O)[C@H](O[C@@H](O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO)O[C@H]%10[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]%10CO[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)O)[C@H]8O)CO[C@H]%12O[C@@H]([C@@H](O)[C@H](O)[C@@H]%12O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O[C@]%15(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%15)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%14O)CO)[C@H](O)[C@H]%13NC(=O)C)CO)CO[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO[C@]%18(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%18)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%16NC(=O)C)CO',
                                     'name': 'CID 91849359',
                                     'reason': 'No phytosphingosine backbone '
                                               'found'},
                                 {   'smiles': 'O=C1OCC=2[C@@]1(O)[C@@]3([C@H](C(CCC3)(C)C)[C@@H](C2)OC(=O)CCCCCCC)C',
                                     'name': 'Nanangenine E',
                                     'reason': 'No phytosphingosine backbone '
                                               'found'},
                                 {   'smiles': 'Nc1ccn([C@H]2C[C@H](O)[C@@H](COP([O-])([O-])=O)O2)c(=O)n1',
                                     'name': "2'-deoxycytosine "
                                             "5'-monophosphate(2-)",
                                     'reason': 'No phytosphingosine backbone '
                                               'found'},
                                 {   'smiles': 'O=C1OC(CC(=O)N[C@H](C(=O)N[C@H](C(=O)OC(CC(N[C@H](C(N[C@H]1CO)=O)CC(C)C)=O)CCCCCCC)CC(C)C)CO)CCCCC',
                                     'name': 'Icosalide A2',
                                     'reason': 'No phytosphingosine backbone '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': '[C@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)OC[C@@H]([C@@H]([C@@H](CCCCCCCCCCCCCC)O)O)NC(CCCCCCCC2=CC=C(C=C2)C(F)(F)F)=O',
                                      'name': '1-O-(alpha-D-galactopyranosyl)-N-[8-(4-trifluoromethyl)phenyloctanoyl]phytosphingosine',
                                      'reason': 'No acyl group attached to '
                                                'nitrogen found'},
                                  {   'smiles': 'C(CCCCCCCCCC)CCC[C@H]([C@@H](O)[C@@H](NC(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)CO)O',
                                      'name': 'N-(13Z-docosenoyl)-phytosphingosine',
                                      'reason': 'No galactose or glucose '
                                                'moiety found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC',
                                      'name': 'N-octacosanoylphytosphingosine',
                                      'reason': 'No galactose or glucose '
                                                'moiety found'},
                                  {   'smiles': '[C@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)OC[C@@H]([C@@H]([C@@H](CCCCCCCCCCCCCC)O)O)NC(CCCCCC2=CC=C(C=C2)OC)=O',
                                      'name': '1-O-(alpha-D-galactopyranosyl)-N-[6-(4-methoxyphenyl)hexanoyl]phytosphingosine',
                                      'reason': 'No acyl group attached to '
                                                'nitrogen found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC',
                                      'name': 'N-hexadecanoylphytosphingosine',
                                      'reason': 'No galactose or glucose '
                                                'moiety found'},
                                  {   'smiles': '[C@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)OC[C@@H]([C@@H]([C@@H](CCCCCCCCCCCCCC)O)O)NC(CCCCCC2=CC=C(C=C2)C(F)(F)F)=O',
                                      'name': '1-O-(alpha-D-galactopyranosyl)-N-{6-[4-(trifluoromethyl)phenyl]hexanoyl}phytosphingosine',
                                      'reason': 'No acyl group attached to '
                                                'nitrogen found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC',
                                      'name': 'N-tetracosanoylphytosphingosine',
                                      'reason': 'No galactose or glucose '
                                                'moiety found'},
                                  {   'smiles': 'CCCCCCCCCCCCCC[C@@H](O)[C@@H](O)[C@H](CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCc1ccccc1',
                                      'name': '1-O-(alpha-D-galactopyranosyl)-N-(8-phenyloctanoyl)phytosphingosine',
                                      'reason': 'No acyl group attached to '
                                                'nitrogen found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC',
                                      'name': 'N-eicosanoyl-4-hydroxysphinganine',
                                      'reason': 'No galactose or glucose '
                                                'moiety found'},
                                  {   'smiles': 'O[C@@H](CCCCCC/C=C/[C@H](O)[C@@H](O)[C@H](O)[C@H](NC(=O)C)C(O)=O)CCCCCC',
                                      'name': 'Sphingofungin D',
                                      'reason': 'No phytosphingosine backbone '
                                                'found'}],
    'attempt': 2,
    'reasoning': "Here's my attempt at a program to classify "
                 'N-acylphytosphingosines:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 49,
    'num_false_positives': 32,
    'num_true_negatives': 142174,
    'num_false_negatives': 45,
    'num_negatives': None,
    'precision': 0.6049382716049383,
    'recall': 0.5212765957446809,
    'f1': 0.56,
    'accuracy': 0.9994588896697119,
    'negative_predictive_value': 0.9996835865812584}