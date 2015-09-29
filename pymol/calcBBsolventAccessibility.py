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
    for rID in residueIDs:
        #print rID
        print rID,cmd.get_area("%s and resid %s and name N"%(sel1,rID)),cmd.get_area("%s and resid %s and name C"%(sel1,rID)) 
    return 0

cmd.extend( "BBsolvent", BBsolvent );
