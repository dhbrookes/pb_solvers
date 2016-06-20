cimport PBAMWrap
cimport PBAMStruct
from cython.view cimport array as cvarray

import logging

_log = logging.getLogger()

cdef class PBAM_Solver:
  #cdef int _ngiter
  cdef double temp_ # system temperature
  cdef double epsilons_ # solvent dielectric
  cdef double epsiloni_  # molecule dielectric

  def __cinit__(self):
    self.temp_ = 0
    self.epsilons_ = 0
    self.epsiloni_ = 0


  def __init__(self, temp, epsilons, epsiloni):
    self.temp_ = temp
    self.epsilons_ = epsilons
    self.epsiloni_ = epsiloni

    import os
    self._pid = os.getpid()

  cdef _run_pbam(self, molecules):
    cdef int i
    cdef int[:] natm
    cdef double[:, :, :] xyzrc
    cdef PBAMStruct.PBAMInput pbamin
    cdef PBAMStruct.PBAMOutput pbamout

    #TODO: change from 1 molecule to many
    nmol = 1
    natm = cvarray(shape=(nmol), itemsize=sizeof(int), format="i")
    natm[0] = <int> len(molecules['atoms'])
    xyzrc = cvarray(shape=(nmol, PBAMStruct.AT_MAX, PBAMStruct.XYZRCWIDTH),
        itemsize=sizeof(double), format="d")

    #TODO: change from 1 molecule to many
    for i, atom in enumerate(molecules['atoms']):
      xyzrc[0, i, 0] = atom['pos'][0]
      xyzrc[0, i, 1] = atom['pos'][1]
      xyzrc[0, i, 2] = atom['pos'][2]
      xyzrc[0, i, 3] = atom['radius']
      xyzrc[0, i, 4] = atom['charge']

    pbamin.temp_ = self.temp_;
    pbamin.sdiel_ = self.epsilons_;
    pbamin.idiel_ = self.epsiloni_;

    pbamout = PBAMWrap.runPBAMSphinxWrap(
      # trying to figure out next line... seems like a pointer to the xyzr vect
        <double (*)[PBAMStruct.AT_MAX][PBAMStruct.XYZRCWIDTH]> &xyzrc[0,0,0],
        nmol, <int (*)> &natm[0], pbamin)

    return {'energy': pbamout.energies_,
        'force': pbamout.forces_}

  def run_solv(self, molecules):
    '''Method called from plugin.py, seems to be wrapping
       to internal stuff'''
    results = []
    _log.info("({}): PBAM flow solver started.".format(self._pid))
    results.append(self._run_pbam(molecules))
    _log.info("({}): PBAM flow solver done.".format(self._pid))
    return results
