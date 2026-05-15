"""
Classifies: CHEBI:33303 chalcogen
"""
"""
Classifies: chalcogen
A chalcogen is any p-block element belonging to the group 16 family of the periodic table.
"""
from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a chemical entity is a chalcogen based on its SMILES string.
    Only identifies pure chalcogen atoms, not compounds or ions containing chalcogens.
    
    Args:
        smiles (str): SMILES string of the chemical entity
        
    Returns:
        tuple: (bool, str) indicating if it's a chalcogen and the reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if it's a single atom
    if mol.GetNumAtoms() != 1:
        return False, "Not a single atom"
    
    # Get the atom
    atom = mol.GetAtomWithIdx(0)
    
    # List of chalcogen atomic numbers
    chalcogens = {
        8: "oxygen",
        16: "sulfur",
        34: "selenium",
        52: "tellurium",
        84: "polonium"
    }
    
    atomic_num = atom.GetAtomicNum()
    
    # Check if it's a chalcogen
    if atomic_num not in chalcogens:
        return False, f"Not a chalcogen (atomic number {atomic_num})"
    
    # Check for formal charge - pure atoms should have no charge
    if atom.GetFormalCharge() != 0:
        return False, "Not a pure atom (has formal charge)"
        
    # Check for explicit hydrogens
    if atom.GetNumExplicitHs() > 0:
        return False, "Not a pure atom (has explicit hydrogens)"
        
    # Check for implicit hydrogens
    if atom.GetNumImplicitHs() > 0:
        return False, "Not a pure atom (has implicit hydrogens)"
        
    # Check for bonds
    if len(atom.GetBonds()) > 0:
        return False, "Not a pure atom (has bonds)"
    
    # Get isotope information if present
    isotope = atom.GetIsotope()
    element_name = chalcogens[atomic_num]
    
    if isotope:
        return True, f"Chalcogen: {element_name}-{isotope}"
    else:
        return True, f"Chalcogen: {element_name}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33303',
                          'name': 'chalcogen',
                          'definition': 'Any p-block element belonging to the '
                                        'group 16 family of the periodic '
                                        'table.',
                          'parents': ['CHEBI:33560'],
                          'xrefs': ['PMID:17084588'],
                          'all_positive_examples': []},
    'config': None,
    'code_statistics': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.6782608695652175 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: [206Po] NAME: polonium-206 atom '
               'REASON: CORRECT Chalcogen: polonium-206\n'
               ' * SMILES: [35S] NAME: sulfur-35 atom REASON: CORRECT '
               'Chalcogen: sulfur-35\n'
               ' * SMILES: [36S] NAME: sulfur-36 atom REASON: CORRECT '
               'Chalcogen: sulfur-36\n'
               ' * SMILES: [S] NAME: sulfur atom REASON: CORRECT Chalcogen: '
               'sulfur\n'
               ' * SMILES: [190Po] NAME: polonium-190 atom REASON: CORRECT '
               'Chalcogen: polonium-190\n'
               ' * SMILES: [198Po] NAME: polonium-198 atom REASON: CORRECT '
               'Chalcogen: polonium-198\n'
               ' * SMILES: [201Po] NAME: polonium-201 atom REASON: CORRECT '
               'Chalcogen: polonium-201\n'
               ' * SMILES: [204Po] NAME: polonium-204 atom REASON: CORRECT '
               'Chalcogen: polonium-204\n'
               ' * SMILES: [209Po] NAME: polonium-209 atom REASON: CORRECT '
               'Chalcogen: polonium-209\n'
               ' * SMILES: [202Po] NAME: polonium-202 atom REASON: CORRECT '
               'Chalcogen: polonium-202\n'
               ' * SMILES: [200Po] NAME: polonium-200 atom REASON: CORRECT '
               'Chalcogen: polonium-200\n'
               ' * SMILES: [191Po] NAME: polonium-191 atom REASON: CORRECT '
               'Chalcogen: polonium-191\n'
               ' * SMILES: [203Po] NAME: polonium-203 atom REASON: CORRECT '
               'Chalcogen: polonium-203\n'
               ' * SMILES: [208Po] NAME: polonium-208 atom REASON: CORRECT '
               'Chalcogen: polonium-208\n'
               ' * SMILES: [17O] NAME: oxygen-17 atom REASON: CORRECT '
               'Chalcogen: oxygen-17\n'
               ' * SMILES: [214Po] NAME: polonium-214 atom REASON: CORRECT '
               'Chalcogen: polonium-214\n'
               ' * SMILES: [216Po] NAME: polonium-216 atom REASON: CORRECT '
               'Chalcogen: polonium-216\n'
               ' * SMILES: [34S] NAME: sulfur-34 atom REASON: CORRECT '
               'Chalcogen: sulfur-34\n'
               ' * SMILES: [Po] NAME: polonium atom REASON: CORRECT Chalcogen: '
               'polonium\n'
               ' * SMILES: [18O] NAME: oxygen-18 atom REASON: CORRECT '
               'Chalcogen: oxygen-18\n'
               ' * SMILES: [217Po] NAME: polonium-217 atom REASON: CORRECT '
               'Chalcogen: polonium-217\n'
               ' * SMILES: [38S] NAME: sulfur-38 atom REASON: CORRECT '
               'Chalcogen: sulfur-38\n'
               ' * SMILES: [210Po] NAME: polonium-210 atom REASON: CORRECT '
               'Chalcogen: polonium-210\n'
               ' * SMILES: [Te] NAME: tellurium atom REASON: CORRECT '
               'Chalcogen: tellurium\n'
               ' * SMILES: [218Po] NAME: polonium-218 atom REASON: CORRECT '
               'Chalcogen: polonium-218\n'
               'False positives: SMILES: [S+6] NAME: sulfur-34(6+) REASON: '
               'WRONGLY CLASSIFIED Chalcogen: sulfur\n'
               ' * SMILES: [Te+][H] NAME: tellanylium REASON: WRONGLY '
               'CLASSIFIED Chalcogen: tellurium\n'
               ' * SMILES: [S-2] NAME: sulfur-34(2-) REASON: WRONGLY '
               'CLASSIFIED Chalcogen: sulfur\n'
               ' * SMILES: [S-] NAME: sulfide(.1-) REASON: WRONGLY CLASSIFIED '
               'Chalcogen: sulfur\n'
               ' * SMILES: [Se][H] NAME: selanyl REASON: WRONGLY CLASSIFIED '
               'Chalcogen: selenium\n'
               ' * SMILES: [H][S+][H] NAME: dihydridosulfur(.1+) REASON: '
               'WRONGLY CLASSIFIED Chalcogen: sulfur\n'
               ' * SMILES: [Se--] NAME: selenide REASON: WRONGLY CLASSIFIED '
               'Chalcogen: selenium\n'
               ' * SMILES: [S+][H] NAME: sulfanylium REASON: WRONGLY '
               'CLASSIFIED Chalcogen: sulfur\n'
               ' * SMILES: [Se-2] NAME: selenium-82(2-) REASON: WRONGLY '
               'CLASSIFIED Chalcogen: selenium\n'
               ' * SMILES: [H][Te]([H])([H])[H] NAME: lambda(4)-tellane '
               'REASON: WRONGLY CLASSIFIED Chalcogen: tellurium\n'
               ' * SMILES: [S--] NAME: sulfide(2-) REASON: WRONGLY CLASSIFIED '
               'Chalcogen: sulfur\n'
               ' * SMILES: [H][S+]([H])[H] NAME: sulfonium REASON: WRONGLY '
               'CLASSIFIED Chalcogen: sulfur\n'
               ' * SMILES: [S+2] NAME: sulfur-34(2+) REASON: WRONGLY '
               'CLASSIFIED Chalcogen: sulfur\n'
               ' * SMILES: [O-][H] NAME: hydroxide REASON: WRONGLY CLASSIFIED '
               'Chalcogen: oxygen\n'
               ' * SMILES: [O][H] NAME: hydroxyl REASON: WRONGLY CLASSIFIED '
               'Chalcogen: oxygen\n'
               ' * SMILES: [H][Te][H] NAME: tellane REASON: WRONGLY CLASSIFIED '
               'Chalcogen: tellurium\n'
               ' * SMILES: [H][Po][H] NAME: polane REASON: WRONGLY CLASSIFIED '
               'Chalcogen: polonium\n'
               ' * SMILES: [H][Se][H] NAME: selane REASON: WRONGLY CLASSIFIED '
               'Chalcogen: selenium\n'
               ' * SMILES: [O-] NAME: oxide(.1-) REASON: WRONGLY CLASSIFIED '
               'Chalcogen: oxygen\n'
               ' * SMILES: [H]O[H] NAME: water REASON: WRONGLY CLASSIFIED '
               'Chalcogen: oxygen\n'
               ' * SMILES: [H][Te]([H])[H] NAME: trihydridotellurium(.) '
               'REASON: WRONGLY CLASSIFIED Chalcogen: tellurium\n'
               ' * SMILES: [Te-][H] NAME: tellanide REASON: WRONGLY CLASSIFIED '
               'Chalcogen: tellurium\n'
               ' * SMILES: [H][Te+]([H])[H] NAME: telluronium REASON: WRONGLY '
               'CLASSIFIED Chalcogen: tellurium\n'
               ' * SMILES: [H][Se]([H])[H] NAME: trihydridoselenium(.) REASON: '
               'WRONGLY CLASSIFIED Chalcogen: selenium\n'
               ' * SMILES: [S-][H] NAME: hydrosulfide REASON: WRONGLY '
               'CLASSIFIED Chalcogen: sulfur\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'Not a single atom'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'Not a single atom'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'Not a single atom'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'Not a single atom'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'Not a single atom'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'Not a single atom'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'Not a single atom'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'Not a single atom'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'Not a single atom'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'Not a single atom'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'Let me analyze the previous attempt and what went wrong:\n'
                 '\n'
                 '1. Main Issue: The program was correctly identifying single '
                 'chalcogen atoms but was also incorrectly classifying '
                 'compounds containing chalcogens (false positives). Looking '
                 'at the false positives, they include:\n'
                 '   - Ions (e.g., [S+6], [S-2], [O-])\n'
                 '   - Hydrides (e.g., [Se][H], [H][Te][H])\n'
                 '   - Compounds (e.g., [H]O[H] water)\n'
                 '\n'
                 '2. Solution Approach:\n'
                 '   - We need to ensure the atom is not bonded to any other '
                 'atoms\n'
                 '   - We need to check that the atom has no charge\n'
                 '   - We need to verify there are no explicit hydrogens '
                 'attached\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 39,
    'num_false_positives': 0,
    'num_true_negatives': 142261,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 1.0,
    'negative_predictive_value': 1.0}