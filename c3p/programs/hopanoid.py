"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:36243 hopanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the hopane skeleton pattern
    hopane_pattern = Chem.MolFromSmarts(
        "[C@H]1CC[C@@]2(C)[C@H]1CC[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3CC[C@@]12C"
    )
    
    # Check if the molecule contains the hopane skeleton
    if not mol.HasSubstructMatch(hopane_pattern):
        return False, "No hopane skeleton found"

    # Check the number of rings to ensure it's a pentacyclic structure
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 5:
        return False, f"Found {n_rings} rings, need at least 5 for hopanoid"

    # Check molecular weight - hopanoids typically have a high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for hopanoid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 30:
        return False, "Too few carbons for hopanoid"
    if o_count > 10:
        return False, "Too many oxygens for typical hopanoid"

    return True, "Contains hopane skeleton with appropriate ring structure and molecular weight"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51963',
                          'name': 'hopanoid',
                          'definition': 'A triterpenoid based on a  hopane '
                                        'skeleton.',
                          'parents': ['CHEBI:36615'],
                          'xrefs': ['KEGG:C06084', 'Wikipedia:Hopanoids'],
                          'all_positive_examples': []},
    'config': None,
    'code_statistics': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cl.C=1C=C(OC)C(OC)=C(C1CN2CCNCC2)OC.Cl',
                                     'name': 'Trimetazidine hydrochloride',
                                     'reason': 'No hopane skeleton found'},
                                 {   'smiles': 'C1[C@@]2([C@]3(C[C@@H]([C@]4([C@]([C@@]3([C@H](C[C@@]2(C[C@@H](C1)O)[H])O)[H])(CC[C@]4([H])[C@@H](CCC(NCCS([O-])(=O)=O)=O)C)[H])C)O)[H])C',
                                     'name': 'tauroursocholate',
                                     'reason': 'No hopane skeleton found'},
                                 {   'smiles': '[H]S(=O)(=O)C[C@H](N)C(O)=O',
                                     'name': 'L-cysteine-S-dioxide',
                                     'reason': 'No hopane skeleton found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)CCCN2C=C(CO[C@@H]1CN(C)CC3=CC=CC=C3OC)N=N2)[C@H](C)CO',
                                     'name': '(8R,9S)-6-[(2R)-1-hydroxypropan-2-yl]-9-[[(2-methoxyphenyl)methyl-methylamino]methyl]-8-methyl-10-oxa-1,6,13,14-tetrazabicyclo[10.2.1]pentadeca-12(15),13-dien-5-one',
                                     'reason': 'No hopane skeleton found'},
                                 {   'smiles': 'Cc1ccc(cc1)N=C=Nc1ccc(C)cc1',
                                     'name': '1,3-di(p-tolyl)carbodiimide',
                                     'reason': 'No hopane skeleton found'},
                                 {   'smiles': 'O=C1[C@H]([C@H]2[C@](C3=C([C@]4([C@]([C@@H]([C@@H]([C@H](O)CC(=C)C(C)C)C)CC4)(C)CC3)C)CC2)(C)CC1)C',
                                     'name': 'Gilvsin C',
                                     'reason': 'No hopane skeleton found'},
                                 {   'smiles': 'C[C@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=CC=C3)C(=O)N(C[C@H]1OC)C)C)C(=O)C',
                                     'name': '1-[(4S,7S,8S)-5-acetyl-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]-3-phenylurea',
                                     'reason': 'No hopane skeleton found'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@@H]1OC)C)C)C(=O)C3=NC=CN=C3',
                                     'name': 'N-[(4S,7R,8R)-8-methoxy-4,7,10-trimethyl-11-oxo-5-[oxo(2-pyrazinyl)methyl]-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No hopane skeleton found'},
                                 {   'smiles': 'COc1ccc(cc1)\\C=C(NC=O)\\C(NC=O)=C\\c1ccc(OS(O)(=O)=O)cc1',
                                     'name': 'Fumiformamide',
                                     'reason': 'No hopane skeleton found'},
                                 {   'smiles': 'CCCCC(CC)COCCCN',
                                     'name': '3-(2-ethylhexoxy)propan-1-amine',
                                     'reason': 'No hopane skeleton found'}],
    'sample_false_negatives': [   {   'smiles': 'O(C1[C@@](O)([C@@H](O)[C@@H]([C@H]1N)O)CO)C[C@H](O)[C@H](O)[C@H](O)CCC([C@@H]2C3[C@](C4[C@@]([C@@]5(C=CC6C(CCC[C@@]6(C5CC4)C)(C)C)C)(C)CC3)(C)CC2)C',
                                      'name': 'Bacteriohop-6-enetetrol '
                                              'carbapseudopentose ether',
                                      'reason': 'No hopane skeleton found'},
                                  {   'smiles': 'O=C1C=CC(C)(C)C2[C@]1([C@@H]3[C@]([C@]4([C@@H]([C@@]5([C@H]([C@@H](C(O)(C)C)CC5)CC4)C)CC3)C)(CC2)C)C',
                                      'name': '22-Hydroxy-2-hopen-1-one',
                                      'reason': 'No hopane skeleton found'},
                                  {   'smiles': 'O[C@@H]1[C@@H]2[C@@]([C@H]3C[C@@H](O)[C@H]4[C@]([C@@]3(C1)C)(CCC=5[C@@]4(CCC5C(C)C)C)C)(CCCC2(C)C)C',
                                      'name': '17(21)-hopene-6alpha,12beta-diol',
                                      'reason': 'No hopane skeleton found'},
                                  {   'smiles': 'O=C(O[C@@H]1[C@](O)([C@H](OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)CCC([C@@H]2[C@H]3[C@]([C@@H]4[C@@]([C@]5([C@@H]([C@@]6([C@H](C(CCC6)(C)C)CC5)C)CC4)C)(C)CC3)(C)CC2)C)[C@@H]([C@H]1OC(=O)C)NC(=O)C)COC(=O)C)C',
                                      'name': '[(1R,2R,3R,4R,5S)-2-[7-[(3R,3aS,5aR,5bR,7aS,11aS,11bR,13aR,13bS)-5a,5b,8,8,11a,13b-hexamethyl-1,2,3,3a,4,5,6,7,7a,9,10,11,11b,12,13,13a-hexadecahydrocyclopenta[a]chrysen-3-yl]-2,3,4-triacetyloxyoctoxy]-3-acetamido-4,5-diacetyloxy-1-hydroxycyclopentyl]methyl '
                                              'acetate',
                                      'reason': 'Too many oxygens for typical '
                                                'hopanoid'},
                                  {   'smiles': 'O=C1C(C2[C@]([C@@H]3[C@@]([C@]4([C@@H]([C@@]5(C(=C(C(C)C)CC5)CC4)C)CC3)C)(C)CC2)(C)CC1)(C)C',
                                      'name': '17(21)-Hopen-3-one',
                                      'reason': 'No hopane skeleton found'},
                                  {   'smiles': 'OCCCC([C@@H]1[C@H]2[C@]([C@@H]3[C@@]([C@@]4(C=C[C@H]5C(CCC[C@@]5([C@H]4CC3)C)(C)C)C)(C)CC2)(C)CC1)C',
                                      'name': 'Hopene',
                                      'reason': 'No hopane skeleton found'},
                                  {   'smiles': 'O=C(O[C@@H]1[C@](O)([C@H](OCC(OC(=O)C)C(OC(=O)C)C(OC(=O)C)CCC([C@@H]2[C@H]3[C@]([C@@H]4[C@@]([C@]5([C@@H]([C@@]6([C@H](C(C[C@H](C6)C)(C)C)CC5)C)CC4)C)(C)CC3)(C)CC2)C)[C@@H]([C@H]1OC(=O)C)NC(=O)C)COC(=O)C)C',
                                      'name': '[(1R,2R,3R,4R,5S)-2-[7-[(3R,3aS,5aR,5bR,7aS,10R,11aS,11bR,13aR,13bS)-5a,5b,8,8,10,11a,13b-heptamethyl-1,2,3,3a,4,5,6,7,7a,9,10,11,11b,12,13,13a-hexadecahydrocyclopenta[a]chrysen-3-yl]-2,3,4-triacetyloxyoctoxy]-3-acetamido-4,5-diacetyloxy-1-hydroxycyclopentyl]methyl '
                                              'acetate',
                                      'reason': 'Too many oxygens for typical '
                                                'hopanoid'},
                                  {   'smiles': 'O(C1C(O)(C(O)C(C1N)O)CO)CC(O)C(O)C(O)C(O)CC([C@@H]2C3[C@](C4C=CC5[C@@]6(C(C([C@@H](C)CC6)(C)C)CC[C@]5([C@]4(C)CC3)C)C)(C)CC2)C',
                                      'name': '7-[(3R,5aR,5bR,9S,11aS,13bS)-5a,5b,8,8,9,11a,13b-heptamethyl-1,2,3,3a,4,5,6,7,7a,9,10,11,11b,13a-tetradecahydrocyclopenta[a]chrysen-3-yl]-1-[5-amino-2,3,4-trihydroxy-2-(hydroxymethyl)cyclopentyl]oxyoctane-2,3,4,5-tetrol',
                                      'reason': 'No hopane skeleton found'},
                                  {   'smiles': 'O(C1C(O)(C(O)C(C1N)O)CO)CC(O)C(O)C(O)C(CC([C@@H]2C3[C@](C4C=CC5[C@@]6(C(C([C@@H](C)C(C6)C)(C)C)CC[C@]5([C@]4(C)CC3)C)C)(C)CC2)C)C',
                                      'name': '5-[7-[(3R,5aR,5bR,9S,11aS,13bS)-5a,5b,8,8,9,10,11a,13b-octamethyl-1,2,3,3a,4,5,6,7,7a,9,10,11,11b,13a-tetradecahydrocyclopenta[a]chrysen-3-yl]-2,3,4-trihydroxy-5-methyloctoxy]-4-amino-1-(hydroxymethyl)cyclopentane-1,2,3-triol',
                                      'reason': 'No hopane skeleton found'},
                                  {   'smiles': 'O=C(O[C@@H]1C2[C@@]([C@H]3CC[C@H]4[C@]([C@@]3(C1)C)(CCC=5[C@@]4(CCC5C(C)C)C)C)(CCCC2(C)C)C)C',
                                      'name': '6a-Acetoxy-17(21)-hopene',
                                      'reason': 'No hopane skeleton found'}],
    'attempt': 0,
    'reasoning': 'To classify hopanoids, we need to identify the '
                 'characteristic hopane skeleton, which is a pentacyclic '
                 'triterpenoid structure. The hopane skeleton consists of five '
                 'fused rings (four six-membered and one five-membered) with '
                 'specific stereochemistry. We can use RDKit to search for '
                 'this core structure in the SMILES string.\n'
                 '\n'
                 'Here is the Python code to classify hopanoids:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 61,
    'num_false_positives': 5,
    'num_true_negatives': 142220,
    'num_false_negatives': 14,
    'num_negatives': None,
    'precision': 0.9242424242424242,
    'recall': 0.8133333333333334,
    'f1': 0.8652482269503545,
    'accuracy': 0.9998664792691497,
    'negative_predictive_value': 0.9999015706511805}