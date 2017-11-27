import sys, os

from rdkit import Chem
# Using last version of standardiser : https://github.com/flatkinson/standardiser
from standardiser import standardise
timeout = -1

_metal_nof = Chem.MolFromSmarts('[Li,Na,K,Rb,Cs,F,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[N,O,F]')
_metal_non = Chem.MolFromSmarts('[Al,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,Hf,Ta,W,Re,Os,Ir,Pt,Au]~[B,C,Si,P,As,Sb,S,Se,Te,Cl,Br,I,At]')
_metals = ['Al','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','Hf','Ta','W','Re','Os','Ir','Pt','Au','Sn']

def disconnect(mol):
    """
    Adapated from molVS standardizer module. Now it returns the list of metals it has disconnected
    """
    for smarts in [_metal_nof, _metal_non]:
        pairs = mol.GetSubstructMatches(smarts)
        rwmol = Chem.RWMol(mol)
        orders = []
        metals = set([])
        for i, j in pairs:
            metalSymbol = mol.GetAtomWithIdx(i).GetSmarts()
            metals.add(metalSymbol)
            orders.append(int(mol.GetBondBetweenAtoms(i, j).GetBondTypeAsDouble()))
            rwmol.RemoveBond(i, j)
        # Adjust neighbouring charges accordingly
        mol = rwmol.GetMol()
        for n, (i, j) in enumerate(pairs):
            chg = orders[n]
            atom1 = mol.GetAtomWithIdx(i)
            atom1.SetFormalCharge(atom1.GetFormalCharge() + chg)
            atom2 = mol.GetAtomWithIdx(j)
            atom2.SetFormalCharge(atom2.GetFormalCharge() - chg)
    
    return mol, metals

def protonate(mol, pH=7.4):
    Chem.MolToMolFile(mol, 'in.sdf')

    stderrf = open (os.devnull, 'w')
    stdoutf = open (os.devnull, 'w')

    outfile = 'out.sdf'

    call = [self.mokaPath+'blabber_sd', 'in.mol',
            '-p',  str(pH),
            '-o',  outfile]

    try:
        retcode = subprocess.call(call,stdout=stdoutf, stderr=stderrf)
    except:
        return (False, 'Blabber execution error')

    stdoutf.close()
    stderrf.close()

    if 'blabber110' in self.mokaPath: # in old blabber versions, error is reported as '0'
        if retcode == 0:
            return (False, 'Blabber 1.0 execution error')
    else:
        if retcode != 0:
            return (False, 'Blabber execution error')

    if not os.path.exists(outfile):
        return (False, 'Blabber output not found')
    if os.stat(outfile).st_size==0:
        return (False, 'Blabber output is empty')

    return (True, Chem.MolFromMolFile(outfile))

def normalize(inF, outF, singleF, failedF, remove_salts= True, keep_nonorganic= False, verbose=False, pH=7.4) :
      
    count = 0        ## count for the whole dataset
    count_inc = 0    ## count for only included molecules
    count_exc = 0    ## count for only excluded molecules
    all_salts = 0    ## count for entries with only salts / solvent
    fail_sanity = 0  ## count for entries that fail sanity check 
    fail_mol = 0     ## count for entries that fail to create mol object 
    fail_prot = 0    ## count for entries that fail protonation

    header = '%s\n' %('\t'.join(['CAS', 'Component', 'Original smiles', 'smiles']))
    fail_header = '%s\n' %('\t'.join(['CAS', 'Original smiles', 'Error']))

    outF.write(header)
    singleF.write(header)
    failedF.write(fail_header)
    
    for line in inF:
        count += 1
        try:
            cas, smi = line.rstrip().split('\t')
        except:
            print ('Failed parsing line:')
            print (line)
            failedF.write(line.rstrip()+'\tFailed parsing line\n')
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            count_exc += 1
            fail_mol += 1
            failedF.write(line.rstrip()+'\tFailed to create molecule object\n')
            continue

        try:
            #mol = standardise.run(mol, keep_nonorganic= keep_nonorganic, remove_salts= remove_salts)
            succ, mol, err = standardise.run(mol, keep_nonorganic= keep_nonorganic)
        except Exception as err:
            err = '{}'.format(err)
            count_exc += 1
            fail_sanity += 1
            failedF.write('{}\t{}\t{}\n'.format(cas, smi, err))
            continue

        i = 1
        if succ:
            count_inc += 1
            nHA = mol.GetNumHeavyAtoms()
            if nHA < 2:
                singleF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, Chem.MolToSmiles(mol)))
            else:
                outF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, Chem.MolToSmiles(mol)))
                #prot, protMol = protonate(mol, pH)
                #if prot:
                #    outF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, Chem.MolToSmiles(protMol)))
                #else:
                #    failedF.write('{}\t{}\t{}\n'.format(cas, smi, protMol))
                #    fail_prot += 1
        else:
            smis = set([Chem.MolToSmiles(moli) for moli in mol])
            if err == 'Multiple non-salt/solvate components':
                for smii in smis:
                    moli = Chem.MolFromSmiles(smii)
                    nHA = moli.GetNumHeavyAtoms()
                    if nHA < 2:
                        singleF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, smii))
                    else:
                        outF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, smii))
                        #prot, protMol = protonate(Chem.MolFromSmiles(smii), pH)
                        #if prot:
                        #    outF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, Chem.MolToSmiles(protMol)))
                        #else:
                        #    failedF.write('{}\t{}\t{}\n'.format(cas, smi, protMol))
                        #    fail_prot += 1
                    i += 1
                count_inc += 1
            elif err == 'No non-salt/solvate components':
                metal = False
                for smii in smis:
                    moli = Chem.MolFromSmiles(smii)
                    nHA = moli.GetNumHeavyAtoms()
                    if nHA == 1 and moli.GetAtomWithIdx(0).GetSymbol() in _metals:
                        singleF.write('{}\t{}\t{}\t{}\n'.format(cas, i, smi, smii))
                        metal = True
                        i += 1
                if metal:
                    count_inc += 1
                else:
                    count_exc += 1
                    all_salts += 1
                    failedF.write('{}\t{}\t{}\n'.format(cas, smi, err))
    
    os.system('rm in.sdf out.sdf')
    print ('the full dataset = {}'.format(count))
    print ('Molecules normalized = {}'.format(count_inc))
    print ('Molecules excluded = {}'.format(count_exc))
    print ('   Fail RDkit mol object = {}'.format(fail_mol))
    print ('   Fail protonation = {}'.format(fail_prot))
    print ('   Fail sanity check = {}'.format(fail_sanity))
    print ('   Only salts / solvent = {}'.format(all_salts))

def std(mol):
    # Standardize and return a dictionary with the smiles as keys
    # and the molecule object and whether it's a metal ion as values
    stdD = {}
    
    # Check single atom compounds, to see if they are metal ions
    if mol.GetNumAtoms() == 1:
        at = mol.GetAtoms()[0].GetSymbol()
        if at in _metals:
            metal = '[%s]' %at
            stdD[metal] = (None, True)
    else:  
        # Extract metal ions from complex compounds
        comp_mol, metals = disconnect(mol)
        for metal in metals:
            metalmol = MolFromSmiles(metal)
            metal = '[%s]' %metalmol.GetAtoms()[0].GetSymbol()
            stdD[metal] = (None, True)

        # For the rest of the molecule, standardize and add
        (passed, std_cmpds, errmessage) = standardise.run(comp_mol)
        if passed:
            stdD[MolToSmiles(std_cmpds)] = (std_cmpds, False)
        elif errmessage == 'Multiple non-salt/solvate components':
            cmpdD = {}
            for cmpd in std_cmpds:
                inchi = MolToInchi(cmpd)
                cmpdD[inchi] = cmpd
            for inchi in cmpdD:
                cmpd = cmpdD[inchi]
                stdD[MolToSmiles(cmpd)] = (cmpd, False)
                
    return stdD
