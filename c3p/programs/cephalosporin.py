"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: CHEBI:23066 cephalosporin
"""

from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins are beta-lactam antibiotics characterized by a beta-lactam ring fused to a six-membered dihydrothiazine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Lists to store ring indices
    beta_lactam_rings = []
    dihydrothiazine_rings = []

    # Identify beta-lactam and dihydrothiazine rings
    for ring in atom_rings:
        if len(ring) == 4:
            # Potential beta-lactam ring
            num_N = 0
            has_carbonyl = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    num_N += 1
            if num_N == 1:
                # Check for carbonyl group (C=O) in ring
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 6:
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 8:
                                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    has_carbonyl = True
                if has_carbonyl:
                    beta_lactam_rings.append(set(ring))
        elif len(ring) == 6:
            # Potential dihydrothiazine ring
            num_N = 0
            num_S = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    num_N += 1
                elif atom.GetAtomicNum() == 16:
                    num_S += 1
            if num_N == 1 and num_S == 1:
                dihydrothiazine_rings.append(set(ring))

    # Check for fused rings
    for beta_ring in beta_lactam_rings:
        for dihydro_ring in dihydrothiazine_rings:
            shared_atoms = beta_ring & dihydro_ring
            if len(shared_atoms) >= 2:
                # Rings are fused
                return True, "Contains fused beta-lactam and dihydrothiazine rings characteristic of cephalosporins"

    return False, "Cephalosporin core structure not found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23066',
                          'name': 'cephalosporin',
                          'definition': 'A class of beta-lactam antibiotics '
                                        'differing from the penicillins in '
                                        'having a 6-membered, rather than a '
                                        '5-membered, side ring.  Although '
                                        'cephalosporins are among the most '
                                        'commonly used antibiotics in the '
                                        'treatment of routine infections, and '
                                        'their use is increasing over time, '
                                        'they can cause a range of '
                                        'hypersensitivity reactions, from '
                                        'mild, delayed-onset cutaneous '
                                        'reactions to life-threatening '
                                        'anaphylaxis in patients with '
                                        'immunoglobulin E (IgE)-mediated '
                                        'allergy.',
                          'parents': ['CHEBI:38311'],
                          'xrefs': [   'KEGG:C00875',
                                       'PMID:10069359',
                                       'PMID:11936371',
                                       'PMID:12833570',
                                       'PMID:24269048',
                                       'PMID:3320614',
                                       'PMID:6762896',
                                       'PMID:8426246',
                                       'Wikipedia:Cephalosporin'],
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
               '[C@]12(N(C(=C(CS1)/C=C/C[N+](CC(N)=O)(CC)C)C(=O)[O-])C([C@H]2NC(=O)/C(/C3=NSC(=N3)N)=N\\OCF)=O)[H] '
               'NAME: cefluprenam REASON: MISSED Cephalosporin core structure '
               'not found\n'
               ' * SMILES: '
               '[H][C@]12SCC(COC)=C(N1C(=O)[C@H]2NC(=O)C(=N/OC)\\c1csc(N)n1)C(O)=O '
               'NAME: cefpodoxime REASON: MISSED Cephalosporin core structure '
               'not found\n'
               ' * SMILES: '
               'CC1=NN=C(S1)SCC2=C(N3[C@@H]([C@H](C3=O)NC(=O)CN4C=NN=N4)SC2)C(=O)O '
               'NAME: '
               '(6R,7S)-3-[[(5-methyl-1,3,4-thiadiazol-2-yl)thio]methyl]-8-oxo-7-[[1-oxo-2-(1-tetrazolyl)ethyl]amino]-5-thia-1-azabicyclo[4.2.0]oct-2-ene-2-carboxylic '
               'acid REASON: MISSED Cephalosporin core structure not found\n'
               ' * SMILES: '
               'CO[C@]1(NC(=O)C(C(O)=O)c2ccc(O)cc2)[C@H]2OCC(CSc3nnnn3C)=C(N2C1=O)C(O)=O '
               'NAME: moxalactam REASON: MISSED Cephalosporin core structure '
               'not found\n'
               ' * SMILES: '
               '[H][C@]12SCC(CSc3nc(=O)c(=O)[nH]n3C)=C(N1C(=O)[C@H]2NC(=O)C(=N/OC)\\c1csc(N)n1)C(O)=O '
               'NAME: ceftriaxone REASON: MISSED Cephalosporin core structure '
               'not found\n'
               ' * SMILES: '
               'CC(=O)OCC1=C(N2[C@@H]([C@H](C2=O)NC(=O)CC3=CC=CS3)SC1)C(=O)O '
               'NAME: '
               '(6R,7S)-3-(acetyloxymethyl)-8-oxo-7-[(1-oxo-2-thiophen-2-ylethyl)amino]-5-thia-1-azabicyclo[4.2.0]oct-2-ene-2-carboxylic '
               'acid REASON: MISSED Cephalosporin core structure not found\n'
               ' * SMILES: '
               'CC(C)(C(=O)O)ON=C(C1=CSC(=N1)N)C(=O)N[C@H]2C3N(C2=O)C(=C(CS3)C[N+]4=CC=CC=C4)C(=O)O '
               'NAME: '
               '(7R)-7-[[2-(2-amino-1,3-thiazol-4-yl)-2-(2-carboxypropan-2-yloxyimino)acetyl]amino]-8-oxo-3-(pyridin-1-ium-1-ylmethyl)-5-thia-1-azabicyclo[4.2.0]oct-2-ene-2-carboxylic '
               'acid REASON: MISSED Cephalosporin core structure not found\n'
               ' * SMILES: '
               'CC1=NN=C(S1)SCC2=C(N3C(C(C3=O)NC(=O)CN4C=NN=N4)SC2)C(=O)O '
               'NAME: '
               '3-[[(5-methyl-1,3,4-thiadiazol-2-yl)thio]methyl]-8-oxo-7-[[1-oxo-2-(1-tetrazolyl)ethyl]amino]-5-thia-1-azabicyclo[4.2.0]oct-2-ene-2-carboxylic '
               'acid REASON: MISSED Cephalosporin core structure not found\n'
               ' * SMILES: '
               'C=CC1=C(N2[C@@H]([C@@H](C2=O)NC(=O)C(=NO)C3=CSC(=N3)N)SC1)C(=O)O '
               'NAME: '
               '(6R,7R)-7-[[2-(2-amino-4-thiazolyl)-2-hydroxyimino-1-oxoethyl]amino]-3-ethenyl-8-oxo-5-thia-1-azabicyclo[4.2.0]oct-2-ene-2-carboxylic '
               'acid REASON: MISSED Cephalosporin core structure not found\n'
               ' * SMILES: '
               'C1C(=C(N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CS3)C(=O)O)C[N+]4=CC=C(C=C4)C(=O)N '
               'NAME: '
               '(6R,7R)-3-[(4-carbamoyl-1-pyridin-1-iumyl)methyl]-8-oxo-7-[(1-oxo-2-thiophen-2-ylethyl)amino]-5-thia-1-azabicyclo[4.2.0]oct-2-ene-2-carboxylic '
               'acid REASON: MISSED Cephalosporin core structure not found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=C(C=C3)C(=O)N(C)C)O[C@@H]1CN(C)C(=O)NC4=CC=CC=C4F)[C@H](C)CO',
                                     'name': '4-[(4R,5S)-5-[[[(2-fluoroanilino)-oxomethyl]-methylamino]methyl]-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-8-yl]-N,N-dimethylbenzamide',
                                     'reason': 'Cephalosporin core structure '
                                               'not found'},
                                 {   'smiles': 'C[C@H]1[C@H](CN(C(=O)C2=C(C=CC(=C2)C#N)OC[C@H]3[C@H](CC[C@@H](O3)CCN(C1=O)C)OC)C)OC',
                                     'name': 'LSM-10564',
                                     'reason': 'Cephalosporin core structure '
                                               'not found'},
                                 {   'smiles': 'CCC(C(=O)NC1=NN=C(S1)SCC)(C(F)(F)F)C(F)(F)F',
                                     'name': 'N-[5-(ethylthio)-1,3,4-thiadiazol-2-yl]-2,2-bis(trifluoromethyl)butanamide',
                                     'reason': 'Cephalosporin core structure '
                                               'not found'},
                                 {   'smiles': 'O=C(N[C@@H](CCC(=O)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC1=CC=C(O)C=C1)CC(=O)N',
                                     'name': 'Tyr-Asn-Gln',
                                     'reason': 'Cephalosporin core structure '
                                               'not found'},
                                 {   'smiles': 'OC[C@H]1OC(O)[C@@H](O)[C@@H]1O',
                                     'name': 'D-arabinofuranose',
                                     'reason': 'Cephalosporin core structure '
                                               'not found'},
                                 {   'smiles': 'O[C@@H]([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)[H])C)[H])C)C',
                                     'name': '(5Z,7E)-(3S,22R)- '
                                             '24-nor-9,10-seco-5,7,10(19)-cholatriene-3,22-diol',
                                     'reason': 'Cephalosporin core structure '
                                               'not found'},
                                 {   'smiles': 'O1C(C(CC1=O)C)CCCCCC',
                                     'name': '5-Hexyldihydro-4-methyl-2(3H)-furanone',
                                     'reason': 'Cephalosporin core structure '
                                               'not found'},
                                 {   'smiles': 'O=C1O[C@]23C(=O)N[C@H]([C@@H]2[C@H](C)C(=C[C@@H]3C=CC[C@@H]([C@H]([C@H](C=C1)C)O)C)C)CC4=CC=CC=C4',
                                     'name': '7-deoxy-cytochalasin Z7',
                                     'reason': 'Cephalosporin core structure '
                                               'not found'},
                                 {   'smiles': '[C@]123[C@H](O1)[C@](C(C=C2[C@]4([C@](CC3)(C([C@@H](C(C4)=O)O)(C)C)[H])C)=O)(C=C)C',
                                     'name': 'yucalexin P8',
                                     'reason': 'Cephalosporin core structure '
                                               'not found'},
                                 {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CCC(O)=O)[C@@]4(C)[C@@H](O)C[C@]3([H])[C@@]1(C)CC[C@@H](O)C2',
                                     'name': 'deoxycholic acid',
                                     'reason': 'Cephalosporin core structure '
                                               'not found'}],
    'sample_false_negatives': [   {   'smiles': 'CO[C@]1(NC(=O)C(C(O)=O)c2ccc(O)cc2)[C@H]2OCC(CSc3nnnn3C)=C(N2C1=O)C(O)=O',
                                      'name': 'moxalactam',
                                      'reason': 'Cephalosporin core structure '
                                                'not found'}],
    'attempt': 1,
    'reasoning': 'The previous code failed because the SMARTS pattern used to '
                 'identify the cephalosporin core was incorrect and did not '
                 'accurately represent the core structure common to '
                 'cephalosporins. Cephalosporins are characterized by a '
                 'beta-lactam ring fused to a dihydrothiazine ring (a '
                 'six-membered ring containing sulfur and nitrogen atoms). The '
                 'beta-lactam ring is a four-membered ring containing a '
                 'nitrogen atom and a carbonyl group, while the '
                 'dihydrothiazine ring contains both sulfur and nitrogen '
                 'atoms.\n'
                 '\n'
                 'To improve the program, we need to:\n'
                 '\n'
                 '1. **Identify the beta-lactam ring**: Search for '
                 'four-membered rings containing one nitrogen atom and one '
                 'carbonyl group (C=O).\n'
                 '\n'
                 '2. **Identify the dihydrothiazine ring**: Search for '
                 'six-membered rings containing one sulfur atom and one '
                 'nitrogen atom.\n'
                 '\n'
                 '3. **Check if the two rings are fused**: Determine if the '
                 'beta-lactam ring and the dihydrothiazine ring share two '
                 'adjacent atoms, indicating they are fused together.\n'
                 '\n'
                 'By implementing these steps, we can more accurately '
                 'determine if a given molecule is a cephalosporin based on '
                 'its SMILES string.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 83,
    'num_false_positives': 71,
    'num_true_negatives': 142145,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.538961038961039,
    'recall': 0.9880952380952381,
    'f1': 0.6974789915966387,
    'accuracy': 0.9994940267041462,
    'negative_predictive_value': 0.9999929649796688}