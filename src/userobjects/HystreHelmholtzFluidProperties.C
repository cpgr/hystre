//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HystreHelmholtzFluidProperties.h"
#include "BrentsMethod.h"
#include "libmesh/utility.h"

InputParameters
HystreHelmholtzFluidProperties::validParams()
{
  InputParameters params = HystreSinglePhaseFluidProperties::validParams();
  params.addClassDescription("Base class for Helmholtz free energy fluid EOS");
  return params;
}

HystreHelmholtzFluidProperties::HystreHelmholtzFluidProperties(const InputParameters & parameters)
  : HystreSinglePhaseFluidProperties(parameters)
{
}

Real
HystreHelmholtzFluidProperties::rho_from_p_T(Real p, Real T) const
{
  Real rho;
  // Initial estimate of a bracketing interval for the density
  Real rho_lower = 1.0e-2;
  Real rho_upper = 100.0;

  // The density is found by finding the zero of the pressure
  auto pressure_diff = [&p, &T, this](Real x) { return this->p_from_rho_T(x, T) - p; };

  BrentsMethod::bracket(pressure_diff, rho_lower, rho_upper);
  rho = BrentsMethod::root(pressure_diff, rho_lower, rho_upper);

  return rho;
}

void
HystreHelmholtzFluidProperties::rho_from_p_T(
    Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const
{
  rho = this->rho_from_p_T(p, T);

  // Scale the density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;
  const Real da_dd = dalpha_ddelta(delta, tau);
  const Real d2a_dd2 = d2alpha_ddelta2(delta, tau);

  drho_dp = molarMass() / (_R * T * delta * (2.0 * da_dd + delta * d2a_dd2));
  drho_dT =
      rho * (tau * d2alpha_ddeltatau(delta, tau) - da_dd) / T / (2.0 * da_dd + delta * d2a_dd2);
}

Real
HystreHelmholtzFluidProperties::e_from_p_T(Real p, Real T) const
{
  // Require density first
  const Real rho = rho_from_p_T(p, T);
  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  return _R * T * tau * dalpha_dtau(delta, tau) / molarMass();
}

void
HystreHelmholtzFluidProperties::e_from_p_T(
    Real p, Real T, Real & e, Real & de_dp, Real & de_dT) const
{
  e = this->e_from_p_T(p, T);

  // Require density first
  const Real rho = rho_from_p_T(p, T);
  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  const Real da_dd = dalpha_ddelta(delta, tau);
  const Real d2a_dd2 = d2alpha_ddelta2(delta, tau);
  const Real d2a_ddt = d2alpha_ddeltatau(delta, tau);

  de_dp = tau * d2a_ddt / (rho * (2.0 * da_dd + delta * d2a_dd2));
  de_dT = -_R *
          (delta * tau * d2a_ddt * (da_dd - tau * d2a_ddt) / (2.0 * da_dd + delta * d2a_dd2) +
           tau * tau * d2alpha_dtau2(delta, tau)) /
          molarMass();
}

Real
HystreHelmholtzFluidProperties::c_from_p_T(Real p, Real T) const
{
  // Require density first
  const Real rho = rho_from_p_T(p, T);
  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  const Real da_dd = dalpha_ddelta(delta, tau);

  Real w = 2.0 * delta * da_dd + delta * delta * d2alpha_ddelta2(delta, tau);
  w -= Utility::pow<2>(delta * da_dd - delta * tau * d2alpha_ddeltatau(delta, tau)) /
       (tau * tau * d2alpha_dtau2(delta, tau));

  return std::sqrt(_R * T * w / molarMass());
}

Real
HystreHelmholtzFluidProperties::cp_from_p_T(Real p, Real T) const
{
  // Require density first
  const Real rho = rho_from_p_T(p, T);
  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  const Real da_dd = dalpha_ddelta(delta, tau);

  const Real cp = _R *
                  (-tau * tau * d2alpha_dtau2(delta, tau) +
                   Utility::pow<2>(delta * da_dd - delta * tau * d2alpha_ddeltatau(delta, tau)) /
                       (2.0 * delta * da_dd + delta * delta * d2alpha_ddelta2(delta, tau))) /
                  molarMass();

  return cp;
}

Real
HystreHelmholtzFluidProperties::cv_from_p_T(Real p, Real T) const
{
  // Require density first
  const Real rho = rho_from_p_T(p, T);
  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  return -_R * tau * tau * d2alpha_dtau2(delta, tau) / molarMass();
}

Real
HystreHelmholtzFluidProperties::s_from_p_T(Real p, Real T) const
{
  // Require density first
  const Real rho = rho_from_p_T(p, T);
  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  return _R * (tau * dalpha_dtau(delta, tau) - alpha(delta, tau)) / molarMass();
}

void
HystreHelmholtzFluidProperties::s_from_p_T(
    Real p, Real T, Real & s, Real & ds_dp, Real & ds_dT) const
{
  s = this->s_from_p_T(p, T);

  // Require density first
  const Real rho = rho_from_p_T(p, T);
  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  const Real da_dd = dalpha_ddelta(delta, tau);
  const Real da_dt = dalpha_dtau(delta, tau);
  const Real d2a_dd2 = d2alpha_ddelta2(delta, tau);
  const Real d2a_dt2 = d2alpha_dtau2(delta, tau);
  const Real d2a_ddt = d2alpha_ddeltatau(delta, tau);

  ds_dp = tau * (d2a_ddt - da_dd) / (rho * T * (2.0 * da_dd + delta * d2a_dd2));
  ds_dT = -_R * tau * (da_dt - alpha(delta, tau) + tau * (d2a_dt2 - da_dt)) / (molarMass() * T);
}

Real
HystreHelmholtzFluidProperties::h_from_p_T(Real p, Real T) const
{
  // Require density first
  const Real rho = rho_from_p_T(p, T);
  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  return _R * T * (tau * dalpha_dtau(delta, tau) + delta * dalpha_ddelta(delta, tau)) / molarMass();
}

void
HystreHelmholtzFluidProperties::h_from_p_T(
    Real p, Real T, Real & h, Real & dh_dp, Real & dh_dT) const
{
  h = this->h_from_p_T(p, T);

  // Require density first
  const Real rho = rho_from_p_T(p, T);
  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  const Real da_dd = dalpha_ddelta(delta, tau);
  const Real d2a_dd2 = d2alpha_ddelta2(delta, tau);
  const Real d2a_ddt = d2alpha_ddeltatau(delta, tau);

  dh_dp = (da_dd + delta * d2a_dd2 + tau * d2a_ddt) / (rho * (2.0 * da_dd + delta * d2a_dd2));
  dh_dT = _R *
          (delta * da_dd * (1.0 - tau * d2a_ddt / da_dd) * (1.0 - tau * d2a_ddt / da_dd) /
               (2.0 + delta * d2a_dd2 / da_dd) -
           tau * tau * d2alpha_dtau2(delta, tau)) /
          molarMass();
}

Real
HystreHelmholtzFluidProperties::p_from_rho_T(Real rho, Real T) const
{
  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  return _R * rho * T * delta * dalpha_ddelta(delta, tau) / molarMass();
}

Real
HystreHelmholtzFluidProperties::z_from_p_T(Real p, Real T) const
{
  // Require density first
  const Real rho = rho_from_p_T(p, T);

  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  return 1.0 + delta * dalphar_ddelta(delta, tau);
}

void
HystreHelmholtzFluidProperties::z_from_p_T(
    Real p, Real T, Real & z, Real & dz_dp, Real & dz_dT) const
{
  z = this->z_from_p_T(p, T);

  // Require density first
  const Real rho = rho_from_p_T(p, T);

  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  const Real da_dd = dalpha_ddelta(delta, tau);
  const Real d2a_dd2 = d2alpha_ddelta2(delta, tau);
  const Real d2a_ddt = d2alpha_ddeltatau(delta, tau);

  dz_dp =
      (da_dd + delta * d2a_dd2) / (rho * _R * T * (2.0 * da_dd + delta * d2a_dd2) / molarMass());

  dz_dT = (delta * tau / T * (da_dd + delta * d2a_dd2) * (d2a_ddt - da_dd / tau) -
           tau * molarMass() * d2a_ddt / (_R * T * T)) /
              (2.0 * da_dd + delta * d2a_dd2) -
          tau * delta * d2a_ddt / T;
}

Real
HystreHelmholtzFluidProperties::f_from_p_T(Real p, Real T) const
{
  // Require density first
  const Real rho = rho_from_p_T(p, T);

  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  const Real z = z_from_p_T(p, T);

  return p * std::exp(z - 1.0 - std::log(z) + alphar(delta, tau));
}

void
HystreHelmholtzFluidProperties::f_from_p_T(
    Real p, Real T, Real & f, Real & df_dp, Real & df_dT) const
{
  f = this->f_from_p_T(p, T);

  // Require density first
  const Real rho = rho_from_p_T(p, T);

  // Scale the input density and temperature
  const Real delta = rho / criticalDensity();
  const Real tau = criticalTemperature() / T;

  const Real z = z_from_p_T(p, T);
  const Real dar_dd = dalphar_ddelta(delta, tau);
  const Real dar_dt = dalphar_dtau(delta, tau);
  const Real d2ar_dd2 = d2alphar_ddelta2(delta, tau);
  const Real d2ar_ddt = d2alphar_ddeltatau(delta, tau);
  const Real da_dd = dalpha_ddelta(delta, tau);
  const Real d2a_dd2 = d2alpha_ddelta2(delta, tau);
  const Real d2a_ddt = d2alpha_ddeltatau(delta, tau);
  const Real dz_dd = dar_dd + delta * d2ar_dd2;
  const Real dz_dt = delta * d2ar_ddt;
  const Real dp_dd = rho * _R * T / molarMass() * (2.0 * da_dd + delta * d2a_dd2);
  const Real dp_dt =
      _R * rho * criticalTemperature() * delta / (molarMass() * tau) * (d2a_ddt - da_dd / tau);
  const Real df_dd = (dp_dd / p + (dz_dd - dz_dd / z + dar_dd)) * f;
  const Real df_dt = (dp_dt / p + (dz_dt - dz_dt / z + dar_dt)) * f;

  df_dp = df_dd / dp_dd;
  df_dT = tau / T * (df_dd * dp_dt - df_dt * dp_dd) / dp_dd;
}

Real
HystreHelmholtzFluidProperties::psi_from_p_T(Real p, Real T) const
{
  mooseAssert(p >= 0.0,
              "HystreHelmholtzFluidProperties::psi_from_p_T: p must be greater than zero");
  return this->f_from_p_T(p, T) / p;
}

void
HystreHelmholtzFluidProperties::psi_from_p_T(
    Real p, Real T, Real & psi, Real & dpsi_dp, Real & dpsi_dT) const
{
  mooseAssert(p >= 0.0,
              "HystreHelmholtzFluidProperties::psi_from_p_T: p must be greater than zero");

  Real f, df_dp, df_dT;
  this->f_from_p_T(p, T, f, df_dp, df_dT);

  psi = f / p;
  dpsi_dp = df_dp / p - f / p / p;
  dpsi_dT = df_dT / p;
}

Real
HystreHelmholtzFluidProperties::alpha(Real delta, Real tau) const
{
  return alpha0(delta, tau) + alphar(delta, tau);
}

Real
HystreHelmholtzFluidProperties::dalpha_ddelta(Real delta, Real tau) const
{
  return dalpha0_ddelta(delta, tau) + dalphar_ddelta(delta, tau);
}

Real
HystreHelmholtzFluidProperties::dalpha_dtau(Real delta, Real tau) const
{
  return dalpha0_dtau(delta, tau) + dalphar_dtau(delta, tau);
}

Real
HystreHelmholtzFluidProperties::d2alpha_ddelta2(Real delta, Real tau) const
{
  return d2alpha0_ddelta2(delta, tau) + d2alphar_ddelta2(delta, tau);
}

Real
HystreHelmholtzFluidProperties::d2alpha_dtau2(Real delta, Real tau) const
{
  return d2alpha0_dtau2(delta, tau) + d2alphar_dtau2(delta, tau);
}

Real
HystreHelmholtzFluidProperties::d2alpha_ddeltatau(Real delta, Real tau) const
{
  return d2alpha0_ddeltatau(delta, tau) + d2alphar_ddeltatau(delta, tau);
}
