"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.

    A 3-hydroxy fatty acyl-CoA(4-) is an acyl-CoA molecule where the fatty acyl chain has a hydroxyl group at the 3-position,
    and the molecule is deprotonated at the phosphate and diphosphate groups, resulting in a net charge of -4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check net charge
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Total charge is {total_charge}, expected -4"

    # Check for Coenzyme A moiety
    # Simplified CoA pattern
    coenzyme_a_smarts = 'NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc23)[C@H](O)[C@H]1OP([O-])([O-])=O'
    coenzyme_a = Chem.MolFromSmarts(coenzyme_a_smarts)
    if not mol.HasSubstructMatch(coenzyme_a):
        return False, "Coenzyme A moiety not found"

    # Check for thioester linkage
    thioester_smarts = 'C(=O)SCCNC(=O)'
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester):
        return False, "Thioester linkage not found"

    # Check for 3-hydroxy fatty acyl chain
    # The chain attached to the thioester should have a hydroxyl at the 3-position
    hydroxy_chain_smarts = '[C;!R][C@H](O)[C;!R][CX3](=O)SCCNC(=O)'
    hydroxy_chain = Chem.MolFromSmarts(hydroxy_chain_smarts)
    if not mol.HasSubstructMatch(hydroxy_chain):
        return False, "3-hydroxy fatty acyl chain not found"

    return True, "Molecule is a 3-hydroxy fatty acyl-CoA(4-)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:65102',
                          'name': '3-hydroxy fatty acyl-CoA(4-)',
                          'definition': 'An acyl-CoA(4-) oxoanion arising from '
                                        'deprotonation of the phosphate and '
                                        'diphosphate OH groups of any '
                                        '3-hydroxy fatty acyl-CoA; major '
                                        'species at pH 7.3.',
                          'parents': ['CHEBI:77636'],
                          'xrefs': ['PMID:17719544'],
                          'all_positive_examples': []},
    'config': None,
    'code_statistics': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'Total charge is 0, expected '
                                               '-4'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'Total charge is 0, expected '
                                               '-4'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'Total charge is 0, expected '
                                               '-4'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'Total charge is 0, expected '
                                               '-4'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'Total charge is 0, expected '
                                               '-4'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Total charge is 0, expected '
                                               '-4'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'Total charge is 0, expected '
                                               '-4'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'Total charge is 0, expected '
                                               '-4'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'Total charge is 0, expected '
                                               '-4'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'Total charge is 0, expected '
                                               '-4'}],
    'sample_false_negatives': [   {   'smiles': '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)C[C@@H](O)CCCCCCCCCCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O',
                                      'name': '(3S)-3-hydroxytetracosanoyl-CoA',
                                      'reason': 'Total charge is 0, expected '
                                                '-4'},
                                  {   'smiles': '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(C[C@H](CC)O)=O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O',
                                      'name': '(S)-3-hydroxypentanoyl-CoA',
                                      'reason': 'Total charge is 0, expected '
                                                '-4'},
                                  {   'smiles': '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)C[C@@H](O)CCC/C=C\\CCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O',
                                      'name': '(3S,7Z)-3-hydroxytetradec-7-enoyl-CoA(4-)',
                                      'reason': 'Total charge is 0, expected '
                                                '-4'}],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a **3-hydroxy fatty '
                 'acyl-CoA(4-)**, we need to check for the following '
                 'features:\n'
                 '\n'
                 '1. **Coenzyme A moiety**: The molecule must contain the CoA '
                 'structure.\n'
                 '2. **Thioester linkage**: There should be a thioester bond '
                 'between the fatty acyl chain and CoA.\n'
                 '3. **3-hydroxy fatty acyl chain**: The fatty acyl chain '
                 'attached to CoA must have a hydroxyl group at the '
                 '3-position.\n'
                 '4. **Net charge of -4**: Due to deprotonation of phosphate '
                 'groups, the molecule should have an overall charge of -4.\n'
                 '\n'
                 'Here is the code that implements these checks:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 61,
    'num_false_positives': 35,
    'num_true_negatives': 142201,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.6354166666666666,
    'recall': 0.953125,
    'f1': 0.7625000000000001,
    'accuracy': 0.9997329585382994,
    'negative_predictive_value': 0.9999789035470169}