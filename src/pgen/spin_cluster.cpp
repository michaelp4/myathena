//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file spin_cluster.cpp
//  \brief Galaxy cluster with added spin.
//
// REFERENCE: Hernquist 1990.

// C++ headers
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

bool _log_on = false;
void Grav(MeshBlock *pmb, const Real time, const Real dt,
          const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
          AthenaArray<Real> &cons);

void MeshBlock::log_info(std::string msg)
{
  if (log_on)
    std::cout << std::endl
              << "*** " + msg + " ***" << std::endl;
}
void log_info_pmb(MeshBlock *pmb, std::string msg)
{
  // if (pmb->log_on)
  if (_log_on)
    std::cout << std::endl
              << "*** " + msg + " ***" << std::endl;
}

void TemperatureCondition(MeshBlock *pmb, const Real time, const Real dt,
                          const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
                          AthenaArray<Real> &cons)
{
  Real numerator = 0.0;
  Real denominator = 0.0;
  for (int k = pmb->ks; k <= pmb->ke; k++)
  {
    for (int j = pmb->js; j <= pmb->je; j++)
    {
      for (int i = pmb->is; i <= pmb->ie; i++)
      {
        Real x = pmb->pcoord->x1v(i);
        Real y = pmb->pcoord->x2v(j);
        Real z = pmb->pcoord->x3v(k);
        Real den = prim(IDN, k, j, i);
        Real temperature = 2.0 / 3.0 * prim(IEN, k, j, i) / prim(IDN, k, j, i);
        numerator += temperature * den;
        denominator += den;
      }
    }
  }
  try
  {
    //log_info_pmb(pmb, "***temprature average:" + std::to_string(numerator / denominator) + "***\n");

    // if the average of the temprature is under 100 eV end the simulation
    if (numerator / denominator < 100)
    {
      Globals::is_running = false;
    }
  }
  catch (...)
  {
    std::cout << "***error in temperature check***";
  }
}
void Grav(MeshBlock *pmb, const Real time, const Real dt,
          const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
          AthenaArray<Real> &cons)
{
  // Setting the Gravitational constant
  Real G = 0.00430091 * pow(10.0, 7.0); // Units: pc (parsec) / solar mass * (km/s)^2
  Real scale_length = 676;

  // Real x0   = pin->GetOrAddReal("problem","x1_0",0.0);
  // Real y0   = pin->GetOrAddReal("problem","x2_0",0.0);
  // Real z0   = pin->GetOrAddReal("problem","x3_0",0.0);
  Real x0 = 0.0;
  Real y0 = 0.0;
  Real z0 = 0.0;

  for (int k = pmb->ks; k <= pmb->ke; k++)
  {
    for (int j = pmb->js; j <= pmb->je; j++)
    {
      for (int i = pmb->is; i <= pmb->ie; i++)
      {
        Real x = pmb->pcoord->x1v(i);
        Real y = pmb->pcoord->x2v(j);
        Real z = pmb->pcoord->x3v(k);
        Real den = prim(IDN, k, j, i);
        Real rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        Real tot_mass = pow(10.0, 5.0);

        Real dPhi = (G * tot_mass) / pow(rad + scale_length, 2.0);
        Real force = -dPhi * den;
        Real dMomentum = force * dt;
        Real dIM1 = dMomentum * x / rad;
        Real dIM2 = dMomentum * y / rad;
        Real dIM3 = dMomentum * z / rad;
        cons(IM1, k, j, i) += dIM1;
        cons(IM2, k, j, i) += dIM2;
        cons(IM3, k, j, i) += dIM3;
        Real velocity_x = cons(IM1, k, j, i) / den;
        Real velocity_y = cons(IM2, k, j, i) / den;
        Real velocity_z = cons(IM3, k, j, i) / den;
        cons(IEN, k, j, i) += dIM1 * velocity_x + dIM2 * velocity_y + dIM3 * velocity_z;

        // logs:
        // std::cout<<std::endl<<"momentum1:"<<cons(IM1,k,j,i)<<std::endl;
        // std::cout<<std::endl<<"momentum2:"<<cons(IM2,k,j,i)<<std::endl;
        // std::cout<<std::endl<<"momentum3:"<<cons(IM3,k,j,i)<<std::endl;
      }
    }
  }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  if (pin->GetOrAddReal("problem", "add_grav", false))
  {
    log_info_pmb(pblock, "*** Grav function was added ***");
    EnrollUserExplicitSourceFunction(Grav);
  }
  if (pin->GetOrAddReal("problem", "add_temperature_condition", false))
  {
    log_info_pmb(pblock, "*** TemperatureCondition function was added ***");
    EnrollUserExplicitSourceFunction(TemperatureCondition);
  }
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Setting the Gravitational constant
  Real G = 0.00430091 * pow(10.0, 7.0); // Units: kpc (kilo parsec) / 10^10 solar mass * (km/s)^2

  log_on = pin->GetOrAddBoolean("problem", "log_on", true);
  _log_on = log_on;
  Real tot_mass = pin->GetOrAddReal("problem", "tot_mass", pow(10.0, 5.0));
  Real scale_length = pin->GetOrAddReal("problem", "scale_length", 676);
  Real angular_velocity = pin->GetOrAddReal("problem", "angular_velocity", 0.0);

  Real rout = pin->GetReal("problem", "radius");
  Real rin = rout - pin->GetOrAddReal("problem", "ramp", 0.0);
  Real pa = pin->GetOrAddReal("problem", "pamb", 1.0);
  Real da = pin->GetOrAddReal("problem", "damb", 1.0);
  Real prat = pin->GetReal("problem", "prat");
  Real drat = pin->GetOrAddReal("problem", "drat", 1.0);
  Real b0, angle;
  if (MAGNETIC_FIELDS_ENABLED)
  {
    b0 = pin->GetReal("problem", "b0");
    angle = (PI / 180.0) * pin->GetReal("problem", "angle");
  }
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  // get coordinates of center of blast, and convert to Cartesian if necessary
  Real x1_0 = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0 = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0 = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  if (COORDINATE_SYSTEM == "cartesian")
  {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  }
  else
  {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  // log stuff intialization
  // std::string density_log = "";
  // std::string kinetic_energy_log = "";

  // setup uniform ambient medium with spherical over-pressured region
  for (int k = ks; k <= ke; k++)
  {
    for (int j = js; j <= je; j++)
    {
      for (int i = is; i <= ie; i++)
      {
        Real rad = 0;
        // if (COORDINATE_SYSTEM == "cartesian") {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        // }

        // Using the function from Hernquist 1990 for mass density and eq of state
        Real den = (tot_mass / (2 * PI)) * (scale_length / rad) * (1 / pow(rad + scale_length, 3.0));
        phydro->u(IDN, k, j, i) = den;
        Real rad_to_scale_ratio = rad / scale_length;
        Real radial_velocity_avg_squared = ((G * tot_mass) / (12 * scale_length)) * ((12 * rad * pow(rad + scale_length, 3.0) / pow(scale_length, 4.0)) * log((rad + scale_length) / rad) - (rad / (rad + scale_length)) * (25 + 52 * rad_to_scale_ratio + 42 * pow(rad_to_scale_ratio, 2.0) + 12 * pow(rad_to_scale_ratio, 3.0)));
        Real pressure = den * radial_velocity_avg_squared;
        std::cout<<"pressure:"+pressure<< std::endl;
        std::cout<<"radius:"+rad<< std::endl;
        Real velocity_squared = pow(angular_velocity * x, 2.0)+ pow(angular_velocity * y, 2.0);
        Real kinetic_energy = pressure / gm1 + 0.5 * den * velocity_squared;
        phydro->u(IM1, k, j, i) = - angular_velocity * y * den;
        phydro->u(IM2, k, j, i) = angular_velocity * x * den;
        phydro->u(IM3, k, j, i) = 0.0;
        phydro->u(IEN, k, j, i) = kinetic_energy;
        // logs
        // if (k == ks)
        // {
        //   density_log += std::to_string(den) + ",";
        //   kinetic_energy_log += std::to_string(kinetic_energy) + ",";
        // }
      }
      // logs
      // if (k == ks)
      // {
      //   density_log += "\n";
      //   kinetic_energy_log += "\n";
      // }
    }
  }
  // log_info("density matrix:");
  // log_info(density_log);
  // log_info("kinetic enegry matrix:");
  // log_info(kinetic_energy_log);
  // log_info("finished initializing spin_cluster");
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Check radius of sphere to make sure it is round
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  return;
  if (!pin->GetOrAddBoolean("problem", "compute_error", false))
    return;

  // analysis - check shape of the spherical blast wave
  int is = pblock->is, ie = pblock->ie;
  int js = pblock->js, je = pblock->je;
  int ks = pblock->ks, ke = pblock->ke;
  AthenaArray<Real> pr;
  pr.InitWithShallowSlice(pblock->phydro->w, 4, IPR, 1);

  // get coordinate location of the center, convert to Cartesian
  Real x1_0 = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0 = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0 = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  if (COORDINATE_SYSTEM == "cartesian")
  {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  }
  else if (COORDINATE_SYSTEM == "cylindrical")
  {
    x0 = x1_0 * std::cos(x2_0);
    y0 = x1_0 * std::sin(x2_0);
    z0 = x3_0;
  }
  else if (COORDINATE_SYSTEM == "spherical_polar")
  {
    x0 = x1_0 * std::sin(x2_0) * std::cos(x3_0);
    y0 = x1_0 * std::sin(x2_0) * std::sin(x3_0);
    z0 = x1_0 * std::cos(x2_0);
  }
  else
  {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ParameterInput" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // find indices of the center
  int ic, jc, kc;
  for (ic = is; ic <= ie; ic++)
    if (pblock->pcoord->x1f(ic) > x1_0)
      break;
  ic--;
  for (jc = pblock->js; jc <= pblock->je; jc++)
    if (pblock->pcoord->x2f(jc) > x2_0)
      break;
  jc--;
  for (kc = pblock->ks; kc <= pblock->ke; kc++)
    if (pblock->pcoord->x3f(kc) > x3_0)
      break;
  kc--;

  // search pressure maximum in each direction
  Real rmax = 0.0, rmin = 100.0, rave = 0.0;
  int nr = 0;
  for (int o = 0; o <= 6; o++)
  {
    int ios = 0, jos = 0, kos = 0;
    if (o == 1)
      ios = -10;
    else if (o == 2)
      ios = 10;
    else if (o == 3)
      jos = -10;
    else if (o == 4)
      jos = 10;
    else if (o == 5)
      kos = -10;
    else if (o == 6)
      kos = 10;
    for (int d = 0; d < 6; d++)
    {
      Real pmax = 0.0;
      int imax, jmax, kmax;
      if (d == 0)
      {
        if (ios != 0)
          continue;
        jmax = jc + jos, kmax = kc + kos;
        for (int i = ic; i >= is; i--)
        {
          if (pr(kmax, jmax, i) > pmax)
          {
            pmax = pr(kmax, jmax, i);
            imax = i;
          }
        }
      }
      else if (d == 1)
      {
        if (ios != 0)
          continue;
        jmax = jc + jos, kmax = kc + kos;
        for (int i = ic; i <= ie; i++)
        {
          if (pr(kmax, jmax, i) > pmax)
          {
            pmax = pr(kmax, jmax, i);
            imax = i;
          }
        }
      }
      else if (d == 2)
      {
        if (jos != 0)
          continue;
        imax = ic + ios, kmax = kc + kos;
        for (int j = jc; j >= js; j--)
        {
          if (pr(kmax, j, imax) > pmax)
          {
            pmax = pr(kmax, j, imax);
            jmax = j;
          }
        }
      }
      else if (d == 3)
      {
        if (jos != 0)
          continue;
        imax = ic + ios, kmax = kc + kos;
        for (int j = jc; j <= je; j++)
        {
          if (pr(kmax, j, imax) > pmax)
          {
            pmax = pr(kmax, j, imax);
            jmax = j;
          }
        }
      }
      else if (d == 4)
      {
        if (kos != 0)
          continue;
        imax = ic + ios, jmax = jc + jos;
        for (int k = kc; k >= ks; k--)
        {
          if (pr(k, jmax, imax) > pmax)
          {
            pmax = pr(k, jmax, imax);
            kmax = k;
          }
        }
      }
      else
      { // if (d==5) {
        if (kos != 0)
          continue;
        imax = ic + ios, jmax = jc + jos;
        for (int k = kc; k <= ke; k++)
        {
          if (pr(k, jmax, imax) > pmax)
          {
            pmax = pr(k, jmax, imax);
            kmax = k;
          }
        }
      }

      Real xm, ym, zm;
      Real x1m = pblock->pcoord->x1v(imax);
      Real x2m = pblock->pcoord->x2v(jmax);
      Real x3m = pblock->pcoord->x3v(kmax);
      if (COORDINATE_SYSTEM == "cartesian")
      {
        xm = x1m;
        ym = x2m;
        zm = x3m;
      }
      else if (COORDINATE_SYSTEM == "cylindrical")
      {
        xm = x1m * std::cos(x2m);
        ym = x1m * std::sin(x2m);
        zm = x3m;
      }
      else
      { // if (COORDINATE_SYSTEM == "spherical_polar") {
        xm = x1m * std::sin(x2m) * std::cos(x3m);
        ym = x1m * std::sin(x2m) * std::sin(x3m);
        zm = x1m * std::cos(x2m);
      }
      Real rad = std::sqrt(SQR(xm - x0) + SQR(ym - y0) + SQR(zm - z0));
      if (rad > rmax)
        rmax = rad;
      if (rad < rmin)
        rmin = rad;
      rave += rad;
      nr++;
    }
  }
  rave /= static_cast<Real>(nr);

  // use physical grid spacing at center of blast
  Real dr_max;
  Real x1c = pblock->pcoord->x1v(ic);
  Real dx1c = pblock->pcoord->dx1f(ic);
  Real x2c = pblock->pcoord->x2v(jc);
  Real dx2c = pblock->pcoord->dx2f(jc);
  Real dx3c = pblock->pcoord->dx3f(kc);
  if (COORDINATE_SYSTEM == "cartesian")
  {
    dr_max = std::max(std::max(dx1c, dx2c), dx3c);
  }
  else if (COORDINATE_SYSTEM == "cylindrical")
  {
    dr_max = std::max(std::max(dx1c, x1c * dx2c), dx3c);
  }
  else
  { // if (COORDINATE_SYSTEM == "spherical_polar") {
    dr_max = std::max(std::max(dx1c, x1c * dx2c), x1c * std::sin(x2c) * dx3c);
  }
  Real deform = (rmax - rmin) / dr_max;

  // only the root process outputs the data
  if (Globals::my_rank == 0)
  {
    std::string fname;
    fname.assign("blastwave-shape.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = fopen(fname.c_str(), "r")) != NULL)
    {
      if ((pfile = freopen(fname.c_str(), "a", pfile)) == NULL)
      {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl
            << "Blast shape output file could not be opened" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }

      // The file does not exist -- open the file in write mode and add headers
    }
    else
    {
      if ((pfile = fopen(fname.c_str(), "w")) == NULL)
      {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl
            << "Blast shape output file could not be opened" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
    fprintf(pfile, "# Offset blast wave test in %s coordinates:\n", COORDINATE_SYSTEM);
    fprintf(pfile, "# Rmax       Rmin       Rave        Deformation\n");
    fprintf(pfile, "%e  %e  %e  %e \n", rmax, rmin, rave, deform);
    fclose(pfile);
  }

  pr.DeleteAthenaArray();
  return;
}
