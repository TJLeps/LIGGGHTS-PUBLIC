/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    Thomas Leps
    University of Maryland College Park
    tjleps@gmail.com
------------------------------------------------------------------------- */
#include "fix_magnetic.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixMagnetic::FixMagnetic(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix magnetic command");

  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else {
    ex = atof(arg[3]);
    xstyle = CONSTANT;
  }

  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else {
    ey = atof(arg[4]);
    ystyle = CONSTANT;
  }

  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else {
    ez = atof(arg[5]);
    zstyle = CONSTANT;
  }

  maxatom = 0;
  hfield = NULL;
}

/* ---------------------------------------------------------------------- */

FixMagnetic::~FixMagnetic()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  memory->destroy(hfield);
}

/* ---------------------------------------------------------------------- */

int FixMagnetic::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMagnetic::init()
{
  //if (!atom->q_flag) error->all(FLERR,"Fix hfield requires atom attribute q");

  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all(FLERR,
                             "Variable name for fix hfield does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix hfield is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR,
                             "Variable name for fix hfield does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix hfield is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR,
                             "Variable name for fix hfield does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix hfield is invalid style");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  //Neighbor List

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // error checks on coarsegraining
  if(force->cg_active())
    error->cg(FLERR,this->style);
}

/* ---------------------------------------------------------------------- */

void FixMagnetic::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixMagnetic::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ----------------------------------------------------------------------*/


void FixMagnetic::post_force(int vflag)
{
  //std::cout<<"made it this far";
  int i,j,k,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh,*slist;
  double dx, dy, dz, fx, fy, fz, rsq, r, mir, mjr, mumu, A, K, muR;
  double *rad = atom->radius;
  double **x = atom->x;
  double **mu = atom->mu;
  double **f = atom->f;
  double *q = atom->q;
  double p4 = M_PI*4;
  double u = p4*1e-7;
  double C = 3*(20000-1)*.0002*.0002*.0002*p4/3/(20000+2);
  //std::cout<<"C = "<<C<<"\n";
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  // reallocate hfield array if necessary

  if (varflag == ATOM && nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(hfield);
    memory->create(hfield,maxatom,3,"hfield:hfield");
  }
  //std::cout<<"made it this far";
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  // std::cout<<ex<<","<<ey<<","<<ez;
  double muc[inum][3];
  if (varflag == CONSTANT) {
    for (k = 0; k < inum; k++) {
    	  if (mask[k] & groupbit) {
    	    muc[k][0] = mu[k][0];
    	    muc[k][1] = mu[k][1];
    	    muc[k][2] = mu[k][2];
    	    C = 3*(20000-1)*rad[i]*rad[i]*rad[i]*p4/3/(20000+2);
    	    mu[k][0] = C*ex/u;
    	    mu[k][1] = C*ey/u;
    	    mu[k][2] = C*ez/u;
    	  }
    }
    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit) {
        jlist = firstneigh[i];
        jnum = numneigh[i];
        for (jj = 0; jj<jnum; jj++)  {

          j =jlist[jj];
          j &= NEIGHMASK;
          dx = x[i][0] - x[j][0];
          dy = x[i][1] - x[j][1];
          dz = x[i][2] - x[j][2];
          rsq = dx*dx + dy*dy + dz*dz;
          r = sqrt(rsq);
          A = C/p4/r/rsq;
          dx /= r;
          dy /= r;
          dz /= r;
          
          //std::cout<<i<<j<<"rh = "<<dx<<","<<dy<<","<<dz;
          
          mjr = muc[j][0]*dx+muc[j][1]*dy+muc[j][2]*dz;
          mir = muc[i][0]*dx+muc[i][1]*dy+muc[i][2]*dz;
          //std::cout<<"mjdotr = "<<mjr;
          //std::cout<<"mu = "<<mu[i][0]<<","<<mu[i][1]<<","<<mu[i][2];
          
          mu[i][0] += A*(3*mjr*dx-muc[j][0]);
          mu[i][1] += A*(3*mjr*dy-muc[j][1]);
          mu[i][2] += A*(3*mjr*dz-muc[j][2]);
          mu[j][0] += A*(3*mir*dx-muc[i][0]);
          mu[j][1] += A*(3*mir*dy-muc[i][1]);
          mu[j][2] += A*(3*mir*dz-muc[i][2]);
          
          //std::cout<<i<<j<<"mu = "<<mu[i][0]<<","<<mu[i][1]<<","<<mu[i][ 2]<<"\n";
        }

        mumu = mu[i][0]*mu[i][0]+mu[i][1]*mu[i][1]+mu[i][2]*mu[i][2];
        if (mumu > C*C*4/u/u){
        	muR=sqrt(C*C*4/u/u/mumu);
        	mu[i][0]=mu[i][0]*muR;
        	mu[i][1]=mu[i][1]*muR;
        	mu[i][2]=mu[i][2]*muR;
        }
        //std::cout<<"mut = "<<mu[i][0]<<","<<mu[i][1]<<","<<mu[i][2]<<"\n";
      }
    }

      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        if (mask[i] & groupbit) {
          jlist = firstneigh[i];
          jnum = numneigh[i];
          for (jj = 0; jj<jnum; jj++)  {

            j =jlist[jj];
            j &= NEIGHMASK;
            dx = x[i][0] - x[j][0];
            dy = x[i][1] - x[j][1];
            dz = x[i][2] - x[j][2];
            rsq = dx*dx + dy*dy + dz*dz;
            r = sqrt(rsq);
            //std::cout<<"r = "<<r<<" ";

            K = 3e-7/rsq/rsq;

            dx /= r;
            dy /= r;
            dz /= r;

            mir = mu[i][0]*dx+mu[i][1]*dy+mu[i][2]*dz;
            mjr = mu[j][0]*dx+mu[j][1]*dy+mu[j][2]*dz;
            mumu = mu[i][0]*mu[j][0]+mu[i][1]*mu[j][1]+mu[i][2]*mu[j][2];

            //std::cout<<"fm = "<<fx<<","<<fy<<","<<fz<<" ";

            f[i][0] += K*(mir*mu[j][0]+mjr*mu[i][0]+(mumu-5*mjr*mir)*dx);
            f[i][1] += K*(mir*mu[j][1]+mjr*mu[i][1]+(mumu-5*mjr*mir)*dy);
            f[i][2] += K*(mir*mu[j][2]+mjr*mu[i][2]+(mumu-5*mjr*mir)*dz);
            
          }
          //std::cout<<"fmf = "<<fx<<","<<fy<<","<<fz<<"\n";
          //f[i][0] += fx;
          //f[i][1] += fy;
          //f[i][2] += fz;
        }
      //std::cout<<"f = "<<f[i][0]<<","<<f[i][1]<<","<<f[i][ 2]<<" ";
      //std::cout<<"fm = "<<fx<<","<<fy<<","<<fz<<" ";
      //std::cout<<"mu = "<<mu[i][0]<<","<<mu[i][1]<<","<<mu[i][ 2]<<" ";
      //std::cout<<"K = "<<K<<"A = "<<A<<"C = "<<C<<"\n";
    }
  }
  //std::cout<<"\n \n \n \n \n";
  // variable hfield, wrap with clear/add

  /*} else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) ex = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM && hfield)
      input->variable->compute_atom(xvar,igroup,&hfield[0][0],3,0);
    if (ystyle == EQUAL) ey = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM && hfield)
      input->variable->compute_atom(yvar,igroup,&hfield[0][1],3,0);
    if (zstyle == EQUAL) ez = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM && hfield)
      input->variable->compute_atom(zvar,igroup,&hfield[0][2],3,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (xstyle == ATOM) f[i][0] += q[i]*hfield[i][0];
        else f[i][0] += q[i]*ex;
        if (ystyle == ATOM) f[i][1] += q[i]*hfield[i][1];
        else f[i][1] += q[i]*ey;
        if (zstyle == ATOM) f[i][2] += q[i]*hfield[i][2];
        else f[i][2] += q[i]*ez;
      }
  }*/
}

/* ---------------------------------------------------------------------- */

void FixMagnetic::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixMagnetic::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}
