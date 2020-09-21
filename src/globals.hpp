#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file globals.hpp
//  \brief namespace containing external global variables
#include "athena.hpp"

namespace Globals {
  extern int my_rank, nranks;
  // Main loop runnint
  extern bool is_running;
  extern Real tot_mass, 
  scale_length,
  angular_velocity,
  x0,
  y0,
  z0,
  no_cooling_radius_under,
  no_cooling_radius_above,
  E_floor,
  log_on,
  log_up_to_redius,
  add_grav,
  add_temerature_condition,
  cooling_param,
  counter;
}

#endif // GLOBALS_HPP_
