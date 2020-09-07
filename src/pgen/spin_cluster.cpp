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

void log_info(ParameterInput *pin, std::string msg)
{
  if (pin->GetOrAddBoolean("problem", "log_on", true))
    std::cout << std::endl
              << "*** " + msg + " ***" << std::endl;
}

void Grav(MeshBlock *pmb, const Real dt, const AthenaArray<Real> &prim,
          AthenaArray<Real> &cons, Real G,Real tot_mass,Real scale_length,
          Real rad,Real den,Real x,Real y,Real z,Real k,Real j,Real i){
  Real dPhi = (G * tot_mass) / pow(rad + scale_length, 2.0);
  Real force = -dPhi * den;
  Real dMomentum = force * dt;
  Real dIM1 = dMomentum * x / rad;
  Real dIM2 = dMomentum * y / rad;
  Real dIM3 = dMomentum * z / rad;
  cons(IM1, k, j, i) += dIM1;
  cons(IM2, k, j, i) += dIM2;
  cons(IM3, k, j, i) += dIM3;
  Real velocity_x = prim(IVX, k, j, i);
  Real velocity_y = prim(IVY, k, j, i);
  Real velocity_z = prim(IVZ, k, j, i);
  cons(IEN, k, j, i) += dIM1 * velocity_x + dIM2 * velocity_y + dIM3 * velocity_z;
}


void Cooling(AthenaArray<Real> &cons, const Real dt, Real k,Real j,Real i,
             Real den, Real pressure, Real rad, Real cooling_param, Real no_cooling_radius,
             Real E_floor, Real time, Real log_on){
  Real cooled_energy = 2.52 * pow(10.0, 7.0) * pow(den, 1.5) * pow(pressure, 0.5) * dt * cooling_param;
  if(rad <= no_cooling_radius || rad >= 1000 || cons(IEN, k, j, i) < cooled_energy) {
    return;
  }
  Real temperature = 72.8 * pressure / den;
  if(log_on > 0.0 && E_floor > cons(IEN, k, j, i) - cooled_energy){
     std::cout<< "*** Energy:" << cons(IEN, k, j, i)<< std::endl 
              << " Cooled energy:" << cooled_energy << std::endl
              << " Radoius: " << rad << " kpc" <<std::endl
              << " Time: " << time <<  std::endl
              // << " Temperature: " << temperature << " ***" << std::endl
              ;
  }
  cons(IEN, k, j, i) = std::fmax(E_floor, cons(IEN, k, j, i) - cooled_energy);
}
void TempCondition(Mesh* mesh){
  if (mesh->dt < pow(10,-7)){
      Globals::is_running = false;
  }
}

void LogTemp(Real* numerator, Real* denominator, Real den, 
                   Real energy){
  Real temperature = 2.0 / 3.0 * energy / den;
  *numerator += temperature * den;
  *denominator += den;
}
void SpinSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
          const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
          AthenaArray<Real> &cons){
  // Setting the Gravitational constant
  Real G = 0.00430091 * pow(10.0, 7.0); // Units: pc (parsec) / solar mass * (km/s)^2
  ParameterInput *pin =pmb->phydro->pin;
  Real tot_mass = pin->GetOrAddReal("problem", "tot_mass", pow(10.0, 5.0));
  Real scale_length = pin->GetOrAddReal("problem", "scale_length", 676);
  Real angular_velocity = pin->GetOrAddReal("problem", "angular_velocity", 0.0);
  Real x0   = pin->GetOrAddReal("problem","x1_0",0.0);
  Real y0   = pin->GetOrAddReal("problem","x2_0",0.0);
  Real z0   = pin->GetOrAddReal("problem","x3_0",0.0);

  Real add_grav = pin->GetOrAddReal("problem", "add_grav", false);
  Real add_temerature_condition = pin->GetOrAddReal("problem", "add_temperature_condition", false);
  Real cooling_param = pin->GetOrAddReal("problem", "cooling_param", false);
  Real no_cooling_radius = pin->GetOrAddReal("problem", "no_cooling_radius", 0.0);
  Real E_floor = pin->GetOrAddReal("problem", "e_floor", pow(10.0,-5.0));
  Real log_on = pin->GetOrAddBoolean("problem", "log_on", 0.0);

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
        Real pressure = prim(IEN, k, j, i);
        Real rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));

        if (add_grav){
          Grav(pmb, dt, prim, cons, G, tot_mass, scale_length,
          rad, den, x, y, z, k, j, i);
        }
        if (add_temerature_condition){
          TempCondition(pmb->pmy_mesh);
        }
        if (cooling_param){
          Cooling(cons, dt, k, j, i, den, pressure, rad, cooling_param, no_cooling_radius,E_floor, time, log_on);
        }
      }
    }
  }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  log_info(pin, "before adding source function");
  EnrollUserExplicitSourceFunction(SpinSourceFunction);
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Setting the Gravitational constant
  Real G = 0.00430091 * pow(10.0, 7.0); // Units: kpc (kilo parsec) / 10^10 solar mass * (km/s)^2

  Real tot_mass = pin->GetOrAddReal("problem", "tot_mass", pow(10.0, 5.0));
  Real scale_length = pin->GetOrAddReal("problem", "scale_length", 676);
  Real angular_velocity = pin->GetOrAddReal("problem", "angular_velocity", 0.0);
  
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
        Real barion_fraction = 0.17;
        Real den = (tot_mass / (2 * PI)) * (scale_length / rad) * (1 / pow(rad + scale_length, 3.0)) * barion_fraction;
        phydro->u(IDN, k, j, i) = den;
        Real rad_to_scale_ratio = rad / scale_length;
        Real radial_velocity_avg_squared = ((G * tot_mass) / (12 * scale_length)) * 
        ((12 * rad * pow(rad + scale_length, 3.0) / pow(scale_length, 4.0)) * log((rad + scale_length) / rad) - 
        (rad / (rad + scale_length)) * (25 + 52 * rad_to_scale_ratio + 42 * pow(rad_to_scale_ratio, 2.0) + 12 * 
        pow(rad_to_scale_ratio, 3.0)));
        Real pressure = den * radial_velocity_avg_squared;

        Real modi_angular_vel = angular_velocity/(1+exp((rad-6500)/300));

        // Adding ceiling velocity (0.2 of the escape velocity)
        Real escape_velocity = std::sqrt(2*G*tot_mass / (rad + scale_length));
        modi_angular_vel = std::fmin(0.2 * escape_velocity / rad, modi_angular_vel);

        Real velocity_squared = pow(modi_angular_vel * x, 2.0)+ pow(modi_angular_vel * y, 2.0);
        Real total_energy = pressure / gm1 + 0.5 * den * velocity_squared; // Theraml energy + Kinetic energy
        phydro->u(IM1, k, j, i) = - modi_angular_vel * y * den;
        phydro->u(IM2, k, j, i) = modi_angular_vel * x * den;
        phydro->u(IM3, k, j, i) = 0.0;
        phydro->u(IEN, k, j, i) = total_energy;
      }
    }
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Check radius of sphere to make sure it is round
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  return;
}
