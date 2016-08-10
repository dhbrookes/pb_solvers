import sys
import numpy as np
import matplotlib.pyplot as plt

'''
Program to rotate and translate a pdb or pqr
'''
fileName = sys,argv[1]
outFile= fileName+ '.rottrans'
rttrnsNam = sys.argv[2]

def getRotTrans(rottranName):
    """Get rotation and translation
       matrices"""
    rttns = open(rottranName).readlines()
    ct, ct3 = 0, 0

    rot = np.zeros([len(rttns)/3, 3, 3])
    trans = np.zeros([len(rttns)/3, 3])

    for line in rttns:
       temp = line.split()
       for i in range(3):
           rot[ct][ct3][i] = float(temp[i+1])
       trans[ct][ct3] = float(temp[4])
       ct3 += 1
       if ct3 == 3 :
         ct3 = 0
         ct += 1

    return( rot, trans)

#-----------------------------------------------------------------------
def RotPDB(fileName, rot, trans, outName):
    """Rotate and translate a pdb"""
    lines = open(fileName).readlines()
    f = open(outName, 'w')

    pos = np.matrix((0.0,0.0,0.0))
    com = np.matrix((0.0,0.0,0.0))
    ct = 0.

    shp = rot.shape
    nmove = shp[0]
    
    for rt in range(nmove):
        rotation = np.matrix((rot[rt][0], rot[rt][1], rot[rt][2]))
        transl = np.matrix((trans[rt]))
        for line in lines:
            if "ATOM" in line:
                tmp = line.split()
                # For PDB
                if ("pdb" in fileName):
                    pos[0,0], pos[0,1], pos[0,2] = float(tmp[6]), \
                                                            float(tmp[7]), \
                                                            float(tmp[8])
                # For PQR
                else:
                    pos[0,0], pos[0,1], pos[0,2] = float(tmp[5]), \
                                                            float(tmp[6]), \
                                                            float(tmp[7])
                # computing center of mass
                com += pos
                ct += 1.
 
                out = rotation*pos.T
                out += transl.T
                #    print pos
                #    print rot*pos.T
                #    print rot*pos.T + tran.T
 
                f.write("{0:s}{1:8.3f}{2:8.3f}{3:8.3f}{4:s}\n".
                            format( line[0:30], out[0,0], out[1,0],
                                        out[2,0], line[54:-1]))
 
    f.close()
    return com/ct

rt, trns = getRotTrans(rttrnsNam)
#print( "This is rot")
#print(rt)
#print( "This is trns")
#print(trns)
RotPDB(fileName, rt, trns, outFile)


