import numpy as np
import math
import matplotlib.pyplot as plt

'''
Program to write VMD script for
'''
dirName='/Users/lfelberg/PBSAM/pb_solvers/pbam/'\
                    'pbam_test_files/energyforce_test/'
fileName = dirName + 'energyforce.kT_0.05M.out'
outFile= dirName + '2sp_0.05M.kT.'

def cleanString(strg):
    '''Removes , [ and ] from a string'''
    strg = strg.replace('[', '')
    strg = strg.replace(']', '')
    strg = strg.replace(',', '')

    return strg

#-----------------------------------------------------------------------
def FileOpen(fileName):
    """Gets data from output of energy force calc"""
    lines = open(fileName).readlines()

    rad, pos, energy, force, torque = [], [], [], [], []

    for line in lines:
        temp = line.split()
        if 'MOLECULE' in temp[0]:
            rad.append(float(temp[3]))
        elif 'POSITION' in temp[0]:
            temp = [float(cleanString(x)) for x in temp[1:]]
            pos.append(temp)
        elif 'ENERGY' in temp[0]:
            energy.append(float(temp[1]))
        elif 'FORCE' in temp[0]:
            temp = [float(cleanString(x)) for x in temp[2:]]
            force.append(temp)
        elif 'TORQUE' in temp[0]:
            temp = [float(cleanString(x)) for x in temp[2:]]
            torque.append(temp)

    return(rad, pos, energy, force, torque)

rad, ps, nrg, frc, tor = FileOpen(fileName)
scale = 0.25
torscal = .2
rad = np.array(rad)
ps = np.array(ps)
nrg = np.array(nrg)
frc = np.array(frc)/scale
tor = np.array(tor)/torscal

vmd_scr = open(outFile+'force', 'w')
for i in range(len(nrg)):
    pos = [ 0, 0, rad[i]]
    vmd_scr.write('draw color black\n')
    strr = 'draw cylinder {{{0}}} '.format(
        ' '.join(map(str,tuple(ps[i]+pos))))
    tmp = ' '.join(map(str,tuple(ps[i]+pos+frc[i]*0.8)))
    strr += '{{{0:s}}} radius {1}\n'.format(
        tmp, rad[i]*0.2)
    vmd_scr.write(strr)

    strr = 'draw cone {{{0}}} '.format(
        ' '.join(map(str,tuple(ps[i]+pos+frc[i]*0.80))))
    tmp = ' '.join(map(str,tuple(ps[i]+pos+frc[i])))
    strr += '{{{0:s}}} radius {1}\n'.format(
        tmp, rad[i]*0.3)
    vmd_scr.write(strr)

vmd_scr.close()


vmd_scr = open(outFile+'torque', 'w')
for i in range(len(nrg)):
    vmd_scr.write('draw color black\n')
    strr = 'draw cylinder {{{0}}} '.format(
        ' '.join(map(str,tuple(ps[i]))))

    torn = np.linalg.norm(tor[i])/(rad[i]*2.)
    tmp = ' '.join(map(str,tuple(ps[i]+tor[i]/torn*0.8)))
    strr += '{{{0:s}}} radius {1}\n'.format(
        tmp, rad[i]*0.2)
    vmd_scr.write(strr)

    strr = 'draw cone {{{0}}} '.format(
        ' '.join(map(str,tuple(ps[i]+tor[i]/torn*0.8))))
    tmp = ' '.join(map(str,tuple(ps[i]+tor[i]/torn)))
    strr += '{{{0:s}}} radius {1}\n'.format(
        tmp, rad[i]*0.3)
    vmd_scr.write(strr)

    torn = np.linalg.norm(tor[i])
    plan = torn/(rad[i]*1.5)
    sprad = 0.2  ## rad of sphere for curved arr
    crad = 1.1

    npts = 360
    xyVec = np.array([0,0,1])
    xx = np.cross(tor[i], xyVec)
    ang = math.acos(np.dot(xyVec, tor[i]) / (
                                np.linalg.norm(xyVec) *
                                np.linalg.norm(tor[i])))

    print xyVec, xx
    axi = xx/np.linalg.norm(xx)
    print ang, axi

    c = math.cos(ang)
    s = math.sin(ang)
    t = 1.0 - c
    x = axi[0]
    y = axi[1]
    z = axi[2]
    rot = np.array(([t*x*x+c, t*x*y-z*s, t*x*z+y*s],
                            [t*x*y+z*s, t*y*y+c, t*y*z-x*s],
                            [t*x*z-y*s, t*y*z+x*s, t*z*z+c]))

    circXY = np.zeros((npts, 3))
    circROT = np.zeros((npts, 3))
    for pt in range(npts):
        circXY[pt][0] = crad * math.cos(math.radians(pt))
        circXY[pt][1] = crad * math.sin(math.radians(pt))
        circXY[pt][2] = 0.0
        circROT[pt] = np.dot(rot, circXY[pt])
        circXY[pt] += ps[i]+tor[i]/plan
        circROT[pt] += ps[i]+tor[i]/plan

    ct = 0
    for pt in np.arange(0.0, torn-4*sprad, sprad):
        strr = 'draw sphere {{{0}}} '.format(
                ' '.join(map(str,tuple(circROT[ct]))) )
        strr += 'radius {0}\n'.format(sprad*1.1)
        vmd_scr.write(strr)
        ct+=1

    strr = 'draw cone {{{0}}} '.format(
        ' '.join(map(str,tuple(circROT[0]))))
    tmp = ' '.join(map(str,tuple(circROT[330])))
    strr += '{{{0:s}}} radius {1}\n'.format(
        tmp, rad[i]*0.2)
    vmd_scr.write(strr)

vmd_scr.close()
