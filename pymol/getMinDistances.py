# Next time, don't use THR
# -- basicCodeBlock.py
#
from pymol import cmd, stored

AA = {}
residue = {}
one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}


def methylDistances( sel1, sel2 ):
    '''
    DESCRIPTION

    methylDistances finds each methyl in sel1, and returns the shortest distance from that methyl
    to any atom in sel2.


    useage:
    first open getMinDistances.py in pymol,
    then invoke on the PyMOL command line:
    PyMOL> methylDistances sel1, sel2

    '''
   
    print "I will return the minimum distance of each methyl in %s to all atoms in %s as a list."% (sel1, sel2)

    sel1 = sel1 + " AND ((resn ALA and name CB) OR (resn ILE and name CD1) OR (resn LEU and name CD1+CD2) OR (resn MET+MSE and name CE) OR (resn VAL and name CG1+CG2))"
    #sel1 = sel1 + " AND ((resn ALA and name CB) OR (resn ILE and name CD1) OR (resn LEU and name CD1+CD2) OR (resn MET+MSE and name CE) OR (resn THR and name CG2) OR (resn VAL and name CG1+CG2))"

    getMinima(sel1,sel2)
    printOutput()
    return 0
#    return getMinima(sel1, sel2)

def printOutput():
    keys1 = residue.keys()
    keys1.sort()
    for key1 in keys1:
        keys2 = residue[key1].keys()
        keys2.sort()
        for key2 in keys2:
            print "%s%s-%s  %0.4f"%(one_letter[AA[key1]], key1, key2, residue[key1][key2])

def getMinima(moleculeA, moleculeB):
    ids_A = getIds(moleculeA)
    ids_B = getIds(moleculeB)
    for idA in ids_A:
        thisAtom = cmd.get_model("%s and id %s" % (moleculeA, idA)).atom[0]
        resnum = thisAtom.resi
        resname = thisAtom.resn
        atname = thisAtom.name
        if resnum in residue.keys():
            if atname in residue.keys():
                print "Error: %s %s %s already already measured"%(resname,resnum,atname)
            else:
                residue[resnum][atname] = getMinDistance(moleculeA,idA,moleculeB)
        else:
            residue[resnum] = {}
            AA[resnum] = resname
            #print residue
            residue[resnum][atname] = getMinDistance(moleculeA,idA,moleculeB)
        #print resname,resnum,atname,getMinDistance(moleculeA,idA,moleculeB)
    return 0

def getMinDistance(moleculeA, idA, moleculeB):
    ids_B = getIds(moleculeB)
    distances = []
    for idB in ids_B:
        d = distance(moleculeA, idA, moleculeB, idB)
        distances.append(d)
    return min(distances)
 
def distance(a, idA, b, idB):
    atomA = "%s and id %s" % (a, idA)
    atomB = "%s and id %s" % (b, idB)
    #print atomA
    #print atomB
    return cmd.get_distance(atomA, atomB)
 
def getIds(selection):
    my_dict = { 'my_list' : [] }
    cmd.iterate(selection, "my_list.append(ID)", space=my_dict)
    return my_dict['my_list']
 
cmd.extend( "methylDistances", methylDistances );
