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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "compute_collision_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "vector_liggghts.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCollisionAtom::ComputeCollisionAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
	
  //NP modified C.K. begin
  if (narg < 3) error->all(FLERR,"Illegal compute contact/atom command");

  skin = 0.;

  if(narg > 3)
  {
      if (narg < 5) error->all(FLERR,"Illegal compute contact/atom command");
      if(strcmp("skin",arg[3])) error->all(FLERR,"Illegal compute contact/atom command, expecting keyword 'skin'");
      skin = atof(arg[4]);
  }
  //NP modified C.K. end
  nvalues = 4;
  peratom_flag = 1;
  size_peratom_cols = nvalues;
  comm_reverse = 1;

  nmax = 0;
  
  array = NULL;

  // error checks

  if (!atom->sphere_flag && !atom->superquadric_flag)
    error->all(FLERR,"Compute contact/atom requires atom style sphere or superquadric!");
}

/* ---------------------------------------------------------------------- */

ComputeCollisionAtom::~ComputeCollisionAtom()
{
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputeCollisionAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute collision/atom requires a pair style be defined");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"collision/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute collision/atom");

  // need an occasional neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->gran = 1;
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeCollisionAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
  /*NL*/ //if (screen) fprintf(screen,"list ptr %d\n",ptr);
}

/* ---------------------------------------------------------------------- */

void ComputeCollisionAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,radsumsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double normal[3],vrel[3],vreln[3],vrelt[3];
  double vrelnmag0,vreltmag0;
  
  invoked_peratom = update->ntimestep;

  // grow contact array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(array);
    nmax = atom->nmax;    
    memory->create(array,nmax,nvalues,"collision/atom:array");
    array_atom=array;
  }

  // invoke neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute number of contacts for each atom in group
  // contact if distance <= sum of radii
  // tally for both I and J

  double **x = atom->x;
  double **v = atom->v;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for (i = 0; i < nall; i++) 
	  for (j=0;j<nvalues;j++)
			array[i][j] = 0.0;
  

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      radi = radius[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radsum = radi + radius[j] + skin; //NP modified C.K.
        radsumsq = radsum*radsum;
        if (rsq <= radsumsq) {
			// Count collisions
          array[i] [0]+= 1.0;
          array[j][0] += 1.0;
		  // Contact normal
        vectorConstruct3D(normal,delx,dely,delz);        
        vectorNormalize3D(normal);
        // Relative velocity
        vectorSubtract3D(v[i],v[j],vrel);    
        vrelnmag0 = vectorDot3D(vrel,normal); 
        vectorScalarMult3D(normal,vrelnmag0,vreln);
        vectorSubtract3D(vrel,vreln,vrelt);
        vreltmag0=vectorMag3D(vrelt);
        // save for atom i  
        array[i] [1]+=abs(vrelnmag0);
        array[i] [2]+=abs(vreltmag0);
        array[i] [3]+=radius[j];
		// save for atom j
		array[j] [1]+=abs(vrelnmag0);
        array[j] [2]+=abs(vreltmag0);
        array[j] [3]+=radi;
        }
      }
    }
  }

  // communicate ghost atom counts between neighbor procs if necessary

  if (force->newton_pair) comm->reverse_comm_compute(this);
}

/* ---------------------------------------------------------------------- */

int ComputeCollisionAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = array[i][0];
    buf[m++] = array[i][1];
    buf[m++] = array[i][2];
    buf[m++] = array[i][3];
}
  return 1;
}

/* ---------------------------------------------------------------------- */

void ComputeCollisionAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    array[j][0] += buf[m++];
    array[j][1] += buf[m++];
    array[j][2] += buf[m++];
    array[j][3] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCollisionAtom::memory_usage()
{
  double bytes = nvalues*nmax * sizeof(double);
  return bytes;
}
