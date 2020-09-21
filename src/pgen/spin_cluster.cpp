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

void Grav(MeshBlock *pmb, const Real dt, const AthenaArray<Real> &prim,
          AthenaArray<Real> &cons, Real G,
          Real rad,Real x,Real y,Real z,Real k,Real j,Real i){
  Real dPhi = (G * Globals::tot_mass) / pow(rad + Globals::scale_length, 2.0);
  Real den = prim(IDN, k, j, i);
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


void Cooling(AthenaArray<Real> &cons, const AthenaArray<Real> &prim, const Real dt, Real k,Real j,Real i,
             Real rad, Real time,MeshBlock *pmb){
  // Momentum and energy are desities.
  Real prim_rho = prim(IDN, k, j, i); // Primitive mass density
  Real pressure = prim(IPR, k, j, i); // Primitive pressure density
  Real primitive_cooled_energy = 2.52 * pow(10.0, 7.0) * pow(prim_rho, 1.5) * 
                                        pow(pressure, 0.5) * dt * Globals::cooling_param;
  if(rad <= Globals::no_cooling_radius_under || rad >= Globals::no_cooling_radius_above) {
    return;
  }
  Real gamma = pmb->peos->GetGamma();
  Real gm1 = gamma - 1.0;
  Real primative_velocity_squared = (SQR(prim(IVX,k,j,i)) + SQR(prim(IVY,k,j,i)) + SQR(prim(IVZ,k,j,i)));
  Real primative_kinetic_energy = 0.5*prim_rho*primative_velocity_squared;

  Real cons_rho = cons(IDN, k, j, i);
  Real conservative_momnetum_squared = SQR(cons(IM1,k,j,i)) + SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i));
  Real conservative_kinetic_energy = 0.5*conservative_momnetum_squared/cons_rho;
  bool log = true;//(Globals::cooling_param == 0.1 && time > 0.055 && time < 0.06) || (Globals::cooling_param == 1 && time > 0.015 && time < 0.02);
  // if (Globals::log_on > 0 && log &&
  // Globals::E_floor + conservative_kinetic_energy > pressure/gm1 + primative_kinetic_energy - primitive_cooled_energy) {
  //   std::cout << "pressure: "<< pressure <<
  //   " prim_rho: "<< prim_rho <<  
  //   " cons_rho: "<< cons_rho <<  
  //   " rad: "<< rad <<  
  //   " prim_cooled_energy: "<< primitive_cooled_energy <<  
  //   " cons_IEN: "<< cons(IEN, k, j, i) <<  
  //   " cons_kin_energy: "<< conservative_kinetic_energy <<  
  //   " prim_kin_energy: "<< primative_kinetic_energy <<  
  //   " prim_temperature: "<< 72.8 * pressure/prim_rho <<  
  //   " time: "<< time << std::endl;
  // }
  
  // if(Globals::log_on > 0) {
  //   Globals::counter=Globals::counter+1.0;
  //   std::cout << "in cooling function, counter: " <<Globals::counter<< std::endl;
  // }

  cons(IEN, k, j, i) = std::fmax(Globals::E_floor + conservative_kinetic_energy, 
  (pressure/gm1 + primative_kinetic_energy) / (1 + primitive_cooled_energy));
}
void TempCondition(Mesh* mesh){
  if (mesh->dt < pow(10,-8)){
      Globals::is_running = false;
  }
}

void SpinSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
          const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
          AthenaArray<Real> &cons){
  // Setting the Gravitational constant
  Real G = 0.00430091 * pow(10.0, 7.0); // Units: pc (parsec) / solar mass * (km/s)^2
  ParameterInput *pin =pmb->phydro->pin;

  for (int k = pmb->ks; k <= pmb->ke; k++)
  {
    for (int j = pmb->js; j <= pmb->je; j++)
    {
      for (int i = pmb->is; i <= pmb->ie; i++)
      {
        Real x = pmb->pcoord->x1v(i);
        Real y = pmb->pcoord->x2v(j);
        Real z = pmb->pcoord->x3v(k);
        Real rad = std::sqrt(SQR(x - Globals::x0) + SQR(y - Globals::y0) + SQR(z - Globals::z0));
        if (Globals::add_grav){
          Grav(pmb, dt, prim, cons, G, rad, x, y, z, k, j, i);
        }
        if (Globals::add_temerature_condition){
          TempCondition(pmb->pmy_mesh);
        }
        if (Globals::cooling_param){
          Cooling(cons, prim, dt, k, j, i, rad, time, pmb);
        }
      }
    }
  }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
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

  Globals::tot_mass = pin->GetOrAddReal("problem", "tot_mass", pow(10.0, 5.0));
  Globals::scale_length = pin->GetOrAddReal("problem", "scale_length", 676);
  Globals::angular_velocity = pin->GetOrAddReal("problem", "angular_velocity", 0.0);
  Globals::x0   = pin->GetOrAddReal("problem","x1_0",0.0);
  Globals::y0   = pin->GetOrAddReal("problem","x2_0",0.0);
  Globals::z0   = pin->GetOrAddReal("problem","x3_0",0.0);
  Globals::no_cooling_radius_under = pin->GetOrAddReal("problem", "no_cooling_radius_under", 0.0);
  Globals::no_cooling_radius_above = pin->GetOrAddReal("problem", "no_cooling_radius_above", 2000.0);
  Globals::E_floor = pin->GetOrAddReal("problem", "e_floor", pow(10.0,-5.0));
  Globals::log_on = pin->GetOrAddReal("problem", "log_on", 0.0);
  Globals::add_grav = pin->GetOrAddReal("problem", "add_grav", false);
  Globals::add_temerature_condition = pin->GetOrAddReal("problem", "add_temperature_condition", false);
  Globals::cooling_param = pin->GetOrAddReal("problem", "cooling_param", false);
  Globals::counter = 0;  

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
        Real den = (Globals::tot_mass / (2.0 * PI)) * (Globals::scale_length / rad) * (1.0 / pow(rad + Globals::scale_length, 3.0)) * barion_fraction;
        phydro->u(IDN, k, j, i) = den;
        Real rad_to_scale_ratio = rad / Globals::scale_length;
        Real radial_velocity_avg_squared = ((G * Globals::tot_mass) / (12.0 * Globals::scale_length)) * 
        ((12.0 * rad * pow(rad + Globals::scale_length, 3.0) / pow(Globals::scale_length, 4.0)) * log((rad + Globals::scale_length) / rad) - 
        (rad / (rad + Globals::scale_length)) * (25.0 + 52.0 * rad_to_scale_ratio + 42.0 * pow(rad_to_scale_ratio, 2.0) + 12.0 * 
        pow(rad_to_scale_ratio, 3.0)));
        Real pressure = den * radial_velocity_avg_squared;

        Real modi_angular_vel = Globals::angular_velocity/(1.0+exp((rad-6500.0)/300.0));

        // Adding ceiling velocity (0.2 of the escape velocity)
        Real escape_velocity = std::sqrt(2.0 * G * Globals::tot_mass / (rad + Globals::scale_length));
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
