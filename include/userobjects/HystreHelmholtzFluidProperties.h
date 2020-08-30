//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HystreSinglePhaseFluidProperties.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

/**
 * Base class equation of state for fluids that use a Helmholtz free energy
 * alpha(delta, tau), where delta is a scaled density and tau is a scaled temperature.
 *
 * To implement a fluid using the Helmholtz EOS, inherit from this base class and
 * override the method alpha(delta, tau) and its first and second derivatives wrt delta
 * and tau. Thermophysical properties such as enthalpy, internal energy etc, will be
 * calculated from the Helmholtz free energy.
 *
 * Transport properties such as viscosity and thermal conductivity will need to be implemented
 * in any derived class as required.
 */
class HystreHelmholtzFluidProperties : public HystreSinglePhaseFluidProperties
{
public:
  static InputParameters validParams();

  HystreHelmholtzFluidProperties(const InputParameters & parameters);

  virtual Real rho_from_p_T(Real p, Real T) const override;

  virtual void
  rho_from_p_T(Real p, Real T, Real & rho, Real & drho_dp, Real & drho_dT) const override;

  virtual Real e_from_p_T(Real p, Real T) const override;
  virtual void e_from_p_T(Real p, Real T, Real & e, Real & de_dp, Real & de_dT) const override;

  virtual Real c_from_p_T(Real p, Real T) const override;

  virtual Real cp_from_p_T(Real p, Real T) const override;

  using HystreSinglePhaseFluidProperties::cp_from_p_T;

  virtual Real cv_from_p_T(Real p, Real T) const override;

  virtual Real s_from_p_T(Real p, Real T) const override;
  virtual void s_from_p_T(Real p, Real T, Real & s, Real & ds_dp, Real & ds_dT) const override;

  virtual Real h_from_p_T(Real p, Real T) const override;
  virtual void h_from_p_T(Real p, Real T, Real & h, Real & dh_dp, Real & dh_dT) const override;

  /**
   * Pressure as a function of density and temperature
   *
   * @param rho density (kg/m^3)
   * @param T temperature (K)
   * @return pressure (Pa)
   */
  virtual Real p_from_rho_T(Real rho, Real T) const;

  virtual Real z_from_p_T(Real p, Real T) const override;
  virtual void z_from_p_T(Real p, Real T, Real & z, Real & dz_dp, Real & dz_dT) const override;

  virtual Real f_from_p_T(Real p, Real T) const override;
  virtual void f_from_p_T(Real p, Real T, Real & f, Real & df_dp, Real & df_dT) const override;

  virtual Real psi_from_p_T(Real p, Real T) const override;
  virtual void
  psi_from_p_T(Real p, Real T, Real & psi, Real & dpsi_dp, Real & dpsi_dT) const override;

protected:
  /**
   * Ideal gas part of Helmholtz free energy
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return alpha0 ideal gas part of Helmholtz free energy
   */
  virtual Real alpha0(Real delta, Real tau) const = 0;

  /**
   * Derivative of ideal gas part of Helmholtz free energy wrt delta
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return derivative of ideal gas part of Helmholtz free energy wrt delta
   */
  virtual Real dalpha0_ddelta(Real delta, Real tau) const = 0;

  /**
   * Derivative of ideal gas part of Helmholtz free energy wrt tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return derivative of ideal gas part of Helmholtz free energy wrt tau
   */
  virtual Real dalpha0_dtau(Real delta, Real tau) const = 0;

  /**
   * Second derivative of ideal gas part of Helmholtz free energy wrt delta
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of ideal gas part of Helmholtz free energy wrt delta
   */
  virtual Real d2alpha0_ddelta2(Real delta, Real tau) const = 0;

  /**
   * Second derivative of ideal gas part of Helmholtz free energy wrt tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of ideal gas part of Helmholtz free energy wrt tau
   */
  virtual Real d2alpha0_dtau2(Real delta, Real tau) const = 0;

  /**
   * Second derivative of ideal gas part of Helmholtz free energy wrt delta and tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of ideal gas part of Helmholtz free energy wrt delta and tau
   */
  virtual Real d2alpha0_ddeltatau(Real delta, Real tau) const = 0;

  /**
   * Residual part of Helmholtz free energy
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return alpha0 residual part of Helmholtz free energy
   */
  virtual Real alphar(Real delta, Real tau) const = 0;

  /**
   * Derivative of residual gas part of Helmholtz free energy wrt delta
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return derivative of residual part of Helmholtz free energy wrt delta
   */
  virtual Real dalphar_ddelta(Real delta, Real tau) const = 0;

  /**
   * Derivative of residual part of Helmholtz free energy wrt tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return derivative of residual part of Helmholtz free energy wrt tau
   */
  virtual Real dalphar_dtau(Real delta, Real tau) const = 0;

  /**
   * Second derivative of residual part of Helmholtz free energy wrt delta
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of residual part of Helmholtz free energy wrt delta
   */
  virtual Real d2alphar_ddelta2(Real delta, Real tau) const = 0;

  /**
   * Second derivative of residual part of Helmholtz free energy wrt tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of residual part of Helmholtz free energy wrt tau
   */
  virtual Real d2alphar_dtau2(Real delta, Real tau) const = 0;

  /**
   * Second derivative of residual part of Helmholtz free energy wrt delta and tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of residual part of Helmholtz free energy wrt delta and tau
   */
  virtual Real d2alphar_ddeltatau(Real delta, Real tau) const = 0;

  /**
   * Helmholtz free energy
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return alpha Helmholtz free energy
   */
  virtual Real alpha(Real delta, Real tau) const;

  /**
   * Derivative of Helmholtz free energy wrt delta
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return derivative of Helmholtz free energy wrt delta
   */
  virtual Real dalpha_ddelta(Real delta, Real tau) const;

  /**
   * Derivative of Helmholtz free energy wrt tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return derivative of Helmholtz free energy wrt tau
   */
  virtual Real dalpha_dtau(Real delta, Real tau) const;

  /**
   * Second derivative of Helmholtz free energy wrt delta
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of Helmholtz free energy wrt delta
   */
  virtual Real d2alpha_ddelta2(Real delta, Real tau) const;

  /**
   * Second derivative of Helmholtz free energy wrt tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of Helmholtz free energy wrt tau
   */
  virtual Real d2alpha_dtau2(Real delta, Real tau) const;

  /**
   * Second derivative of Helmholtz free energy wrt delta and tau
   *
   * @param delta scaled density (-)
   * @param tau scaled temperature (-)
   * @return second derivative of Helmholtz free energy wrt delta and tau
   */
  virtual Real d2alpha_ddeltatau(Real delta, Real tau) const;
};

#pragma GCC diagnostic pop
