/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(coord/gran,ComputeCollisionAtom)
ComputeStyle(collision/atom,ComputeCollisionAtom)

#else

#ifndef LMP_COMPUTE_COLLISION_ATOM_H
#define LMP_COMPUTE_COLLISION_ATOM_H

#include "compute.h"
#include "fix_property_atom.h"

namespace LAMMPS_NS {

class ComputeCollisionAtom : public Compute {
 public:
  ComputeCollisionAtom(class LAMMPS *, int, char **);
  ~ComputeCollisionAtom();
  virtual void init(); //NP modified C.K.
  void init_list(int, class NeighList *);
  virtual void compute_peratom(); //NP modified C.K.
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

  int nvalues;
 protected: //NP modified C.K.
  int nmax;
  class NeighList *list;
   //0 number of collisions
  //1 sum of normal parts of relative velocity with respect to particles contact (center to center vector)
   //2 sum of tangential parts of relative velocity with respect to particles contact (center to center vector)
   //3 sum of radii that collided with particle
  double **array; 
  double skin; //NP modified C.K.
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute collision/atom requires atom style sphere

Self-explanatory.

E: Compute collision/atom requires a pair style be defined

Self-explantory.

W: More than one compute collision/atom

It is not efficient to use compute collision/atom more than once.

*/
