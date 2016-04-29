import numpy as np
import matplotlib.pyplot as plt

'''
Program to rotate and translate a pdb
'''
dirName='/Users/lfelberg/Desktop/'
fileName = dirName + '1BRS_D.pdb'
outFile= dirName + '1brs.pdb'

#-----------------------------------------------------------------------
def RotPDB(fileName, outName):
    """Rotate and translate a pdb"""
    lines = open(fileName).readlines()
    f = open(outName, 'w')

    rot = np.matrix(((1.0, 0.0, 0.0), (0.0, -1.0, 0.0),
                             (0.0, 0.0, -1.0)))
    tran = np.matrix( (0.0, 0.0, 26.1733) )

    pos = np.matrix((0.0,0.0,0.0))
    com = np.matrix((0.0,0.0,0.0))
    ct = 0.

    for line in lines:
        if "ATOM" in line:
            tmp = line.split()
            pos[0,0], pos[0,1], pos[0,2] = float(tmp[6]), \
                                                    float(tmp[7]), \
                                                    float(tmp[8])
            # computing center of mass
            com += pos
            ct += 1.

            out = rot*pos.T
            out += tran.T
            #if int(tmp[1]) < 3:
            #    print pos
            #    print rot*pos.T
            #    print rot*pos.T + tran.T

            f.write("{0:s}{1:8.3f}{2:8.3f}{3:8.3f}{4:s}\n".
                        format( line[0:30], out[0,0], out[1,0],
                                    out[2,0], line[54:-1]))

        else:
            f.write(line)

    f.close()
    return com/ct

print RotPDB(fileName, outFile)
