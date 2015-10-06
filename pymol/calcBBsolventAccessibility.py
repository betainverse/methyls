# Next time, don't use THR
# -- basicCodeBlock.py
#
from pymol import cmd, stored

def BBsolvent( sel1 ):
    '''
    DESCRIPTION

    calcBBsolventAccessibility calculates the solvent accessibility for each
    backbone N and CO atom in sel1


    useage:
    first open calcBBsolventAccessibility.py in pymol,
    then invoke on the PyMOL command line:
    PyMOL> BBsolvent sel1

    '''
   
    print "Returning the solvent accessibility of each BB N and CO in %s as a table."%sel1

    cmd.set('dot_solvent',1)
    cmd.set('dot_density',4)
    print "resid N C"
    my_dict = { 'my_list' : [] }
    cmd.iterate(sel1, "my_list.append(resi)", space=my_dict)
    residueIDs = sorted(set(map(int,my_dict['my_list'])))
    Nvalues = []
    Cvalues = []
    for rID in residueIDs:
        NASA = cmd.get_area("%s and resid %s and name N"%(sel1,rID))
        CASA = cmd.get_area("%s and resid %s and name C"%(sel1,rID))
        cmd.alter("%s+b and resid %d name N"%(sel1,rID), 'b=%0.4f'%NASA)
        cmd.alter("%s and resid %d name C"%(sel1,rID), 'b=%0.4f'%CASA)
        print rID,NASA,CASA
        Nvalues.append(NASA)
        Cvalues.append(CASA)
    #minN = min(Nvalues)
    #maxN = max(Nvalues)
    #maxC = max(Cvalues)
    #cmd.spectrum("b","blue_white_red","%s or chain B and name N"%sel1, minimum=0,maximum=maxN)
    #cmd.spectrum("b","blue_white_red","%s and name C"%sel1, minimum=0,maximum=maxC)
    return 0

cmd.extend( "BBsolvent", BBsolvent );
