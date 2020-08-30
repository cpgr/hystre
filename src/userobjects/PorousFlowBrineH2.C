//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PorousFlowBrineH2.h"
#include "BrineFluidProperties.h"
#include "HystreSinglePhaseFluidProperties.h"
#include "SinglePhaseFluidProperties.h"
#include "MathUtils.h"
#include "Conversion.h"
#include "libmesh/utility.h"

registerMooseObject("hystreApp", PorousFlowBrineH2);

// defineLegacyParams(PorousFlowBrineH2);

InputParameters
PorousFlowBrineH2::validParams()
{
  InputParameters params = PorousFlowFluidStateMultiComponentBase::validParams();
  params.addRequiredParam<UserObjectName>("brine_fp", "The name of the user object for brine");
  params.addRequiredParam<UserObjectName>("h2_fp", "The name of the user object for H2");
  params.addParam<unsigned int>("salt_component", 2, "The component number of salt");
  params.addClassDescription("Fluid state class for brine and H2");
  return params;
}

PorousFlowBrineH2::PorousFlowBrineH2(const InputParameters & parameters)
  : PorousFlowFluidStateMultiComponentBase(parameters),
    _salt_component(getParam<unsigned int>("salt_component")),
    _brine_fp(getUserObject<BrineFluidProperties>("brine_fp")),
    _h2_fp(getUserObject<HystreSinglePhaseFluidProperties>("h2_fp")),
    _H2O_fp(_brine_fp.getComponent(BrineFluidProperties::WATER)),
    _MH2O(_brine_fp.molarMassH2O()),
    _invMH2O(1.0 / _MH2O),
    _ln_invMH2O(std::log(_invMH2O)),
    _MH2(_h2_fp.molarMass()),
    _Mnacl(_brine_fp.molarMassNaCl()),
    _Zmin(1.0e-4)
{
  // Check that the correct FluidProperties UserObjects have been provided
  if (_h2_fp.fluidName() != "hydrogen")
    paramError("h2_fp", "A valid H2 FluidProperties UserObject must be provided");

  if (_brine_fp.fluidName() != "brine")
    paramError("brine_fp", "A valid Brine FluidProperties UserObject must be provided");

  // Set the number of phases and components, and their indexes
  _num_phases = 2;
  _num_components = 3;
  _gas_phase_number = 1 - _aqueous_phase_number;
  _gas_fluid_component = 3 - _aqueous_fluid_component - _salt_component;

  // Check that _aqueous_phase_number is <= total number of phases
  if (_aqueous_phase_number >= _num_phases)
    paramError("liquid_phase_number",
               "This value is larger than the possible number of phases ",
               _num_phases);

  // Check that _aqueous_fluid_component is <= total number of fluid components
  if (_aqueous_fluid_component >= _num_components)
    paramError("liquid_fluid_component",
               "This value is larger than the possible number of fluid components",
               _num_components);

  // Check that the salt component index is not identical to the liquid fluid component
  if (_salt_component == _aqueous_fluid_component)
    paramError(
        "salt_component",
        "The value provided must be different from the value entered in liquid_fluid_component");

  // Check that _salt_component is <= total number of fluid components
  if (_salt_component >= _num_components)
    paramError("salt_component",
               "The value provided is larger than the possible number of fluid components",
               _num_components);

  _empty_fsp = FluidStateProperties(_num_components);
}

std::string
PorousFlowBrineH2::fluidStateName() const
{
  return "brine-h2";
}

void
PorousFlowBrineH2::thermophysicalProperties(Real pressure,
                                            Real temperature,
                                            Real Xnacl,
                                            Real Z,
                                            unsigned int qp,
                                            std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Check whether the input temperature is within the region of validity
  checkVariables(pressure, temperature);

  // AD versions of primary variables
  DualReal p = pressure;
  Moose::derivInsert(p.derivatives(), _pidx, 1.0);
  DualReal T = temperature;
  Moose::derivInsert(T.derivatives(), _Tidx, 1.0);
  DualReal ZH2 = Z;
  Moose::derivInsert(ZH2.derivatives(), _Zidx, 1.0);
  DualReal X = Xnacl;
  Moose::derivInsert(X.derivatives(), _Xidx, 1.0);

  // Clear all of the FluidStateProperties data
  clearFluidStateProperties(fsp);

  FluidStatePhaseEnum phase_state;
  massFractions(p, T, X, ZH2, phase_state, fsp);

  switch (phase_state)
  {
    case FluidStatePhaseEnum::GAS:
    {
      // Set the gas saturations
      gas.saturation = 1.0;

      // Calculate gas properties
      gasProperties(p, T, fsp);

      break;
    }

    case FluidStatePhaseEnum::LIQUID:
    {
      // Calculate the liquid properties
      const DualReal liquid_pressure = p - _pc.capillaryPressure(1.0, qp);
      liquidProperties(liquid_pressure, T, X, fsp);

      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Calculate the gas and liquid properties in the two phase region
      twoPhaseProperties(p, T, X, ZH2, qp, fsp);

      break;
    }
  }

  // Liquid saturations can now be set
  liquid.saturation = 1.0 - gas.saturation;

  // Save pressures to FluidStateProperties object
  gas.pressure = p;
  liquid.pressure = p - _pc.capillaryPressure(liquid.saturation, qp);
}

void
PorousFlowBrineH2::massFractions(const DualReal & pressure,
                                 const DualReal & temperature,
                                 const DualReal & Xnacl,
                                 const DualReal & Z,
                                 FluidStatePhaseEnum & phase_state,
                                 std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  DualReal XH2 = 0.0;
  DualReal YH2O = 0.0;
  DualReal YH2 = 0.0;

  // If the amount of H2 is less than the smallest solubility, then all H2 will
  // be dissolved, and the equilibrium mass fractions do not need to be computed
  if (Z < _Zmin)
    phase_state = FluidStatePhaseEnum::LIQUID;

  else
  {
    // Equilibrium mass fraction of H2 in liquid and H2O in gas phases
    equilibriumMassFractions(pressure, temperature, Xnacl, XH2, YH2O);

    YH2 = 1.0 - YH2O;

    // Determine which phases are present based on the value of z
    phaseState(Z.value(), XH2.value(), YH2.value(), phase_state);
  }

  // The equilibrium mass fractions calculated above are only correct in the two phase
  // state. If only liquid or gas phases are present, the mass fractions are given by
  // the total mass fraction z
  DualReal XH2O = 0.0;

  switch (phase_state)
  {
    case FluidStatePhaseEnum::LIQUID:
    {
      XH2 = Z;
      YH2 = 0.0;
      XH2O = 1.0 - Z;
      YH2O = 0.0;
      Moose::derivInsert(XH2.derivatives(), _pidx, 0.0);
      Moose::derivInsert(XH2.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(XH2.derivatives(), _Xidx, 0.0);
      Moose::derivInsert(XH2.derivatives(), _Zidx, 1.0);
      Moose::derivInsert(YH2.derivatives(), _pidx, 0.0);
      Moose::derivInsert(YH2.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(YH2.derivatives(), _Xidx, 0.0);
      break;
    }

    case FluidStatePhaseEnum::GAS:
    {
      XH2 = 0.0;
      YH2 = Z;
      YH2O = 1.0 - Z;
      Moose::derivInsert(XH2.derivatives(), _pidx, 0.0);
      Moose::derivInsert(XH2.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(XH2.derivatives(), _Xidx, 0.0);
      Moose::derivInsert(YH2.derivatives(), _pidx, 0.0);
      Moose::derivInsert(YH2.derivatives(), _Tidx, 0.0);
      Moose::derivInsert(YH2.derivatives(), _Xidx, 0.0);
      Moose::derivInsert(YH2.derivatives(), _Zidx, 1.0);
      break;
    }

    case FluidStatePhaseEnum::TWOPHASE:
    {
      // Keep equilibrium mass fractions
      XH2O = 1.0 - XH2;
      break;
    }
  }

  // Save the mass fractions in the FluidStateProperties object
  liquid.mass_fraction[_aqueous_fluid_component] = XH2O;
  liquid.mass_fraction[_gas_fluid_component] = XH2;
  liquid.mass_fraction[_salt_component] = Xnacl;
  gas.mass_fraction[_aqueous_fluid_component] = YH2O;
  gas.mass_fraction[_gas_fluid_component] = YH2;
}

void
PorousFlowBrineH2::gasProperties(const DualReal & pressure,
                                 const DualReal & temperature,
                                 std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Gas density, viscosity and enthalpy are approximated with pure H2 - no correction due
  // to the small amount of water vapor is made
  DualReal H2_density, H2_viscosity;
  _h2_fp.rho_mu_from_p_T(pressure, temperature, H2_density, H2_viscosity);

  const DualReal H2_enthalpy = _h2_fp.h_from_p_T(pressure, temperature);

  // Save the values to the FluidStateProperties object. Note that derivatives wrt z are 0
  gas.density = H2_density;
  gas.viscosity = H2_viscosity;
  gas.enthalpy = H2_enthalpy;

  mooseAssert(gas.density.value() > 0.0, "Gas density must be greater than zero");
  gas.internal_energy = gas.enthalpy - pressure / gas.density;
}

void
PorousFlowBrineH2::liquidProperties(const DualReal & pressure,
                                    const DualReal & temperature,
                                    const DualReal & Xnacl,
                                    std::vector<FluidStateProperties> & fsp) const
{
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];

  // The liquid density includes the density increase due to dissolved H2
  const DualReal brine_density = _brine_fp.rho_from_p_T_X(pressure, temperature, Xnacl);

  // Mass fraction of H2 in liquid phase
  const DualReal XH2 = liquid.mass_fraction[_gas_fluid_component];

  // The partial density of H2
  const DualReal H2_partial_density = _MH2 / molarVolume(pressure, temperature);

  // The liquid density
  const DualReal liquid_density = 1.0 / (XH2 / H2_partial_density + (1.0 - XH2) / brine_density);

  // Assume that liquid viscosity is just the brine viscosity
  const DualReal liquid_viscosity = _brine_fp.mu_from_p_T_X(pressure, temperature, Xnacl);

  // Liquid enthalpy (note: not yet including contribution due to the enthalpy of dissolution)
  const DualReal brine_enthalpy = _brine_fp.h_from_p_T_X(pressure, temperature, Xnacl);

  // Enthalpy of H2
  const DualReal H2_enthalpy = _h2_fp.h_from_p_T(pressure, temperature);

  // Enthalpy of liquid
  const DualReal liquid_enthalpy = (1.0 - XH2) * brine_enthalpy + XH2 * (H2_enthalpy);

  // Save the values to the FluidStateProperties object
  liquid.density = liquid_density;
  liquid.viscosity = liquid_viscosity;
  liquid.enthalpy = liquid_enthalpy;

  mooseAssert(liquid.density.value() > 0.0, "Liquid density must be greater than zero");
  liquid.internal_energy = liquid.enthalpy - pressure / liquid.density;
}

DualReal
PorousFlowBrineH2::saturation(const DualReal & pressure,
                              const DualReal & temperature,
                              const DualReal & Xnacl,
                              const DualReal & Z,
                              std::vector<FluidStateProperties> & fsp) const
{
  auto & gas = fsp[_gas_phase_number];
  auto & liquid = fsp[_aqueous_fluid_component];

  // Approximate liquid density as saturation isn't known yet, by using the gas
  // pressure rather than the liquid pressure. This does result in a small error
  // in the calculated saturation, but this is below the error associated with
  // the correlations. A more accurate saturation could be found iteraviely,
  // at the cost of increased computational expense

  // Gas density
  const DualReal gas_density = _h2_fp.rho_from_p_T(pressure, temperature);

  // Approximate liquid density as saturation isn't known yet
  const DualReal brine_density = _brine_fp.rho_from_p_T_X(pressure, temperature, Xnacl);

  // Mass fraction of H2 in liquid phase
  const DualReal XH2 = liquid.mass_fraction[_gas_fluid_component];

  // The liquid density
  const DualReal liquid_density = brine_density;

  const DualReal YH2 = gas.mass_fraction[_gas_fluid_component];

  // Set mass equilibrium constants used in the calculation of vapor mass fraction
  const DualReal K0 = YH2 / XH2;
  const DualReal K1 = (1.0 - YH2) / (1.0 - XH2);
  const DualReal vapor_mass_fraction = vaporMassFraction(Z, K0, K1);

  // The gas saturation in the two phase case
  const DualReal saturation = vapor_mass_fraction * liquid_density /
                              (gas_density + vapor_mass_fraction * (liquid_density - gas_density));

  return saturation;
}

void
PorousFlowBrineH2::twoPhaseProperties(const DualReal & pressure,
                                      const DualReal & temperature,
                                      const DualReal & Xnacl,
                                      const DualReal & Z,
                                      unsigned int qp,
                                      std::vector<FluidStateProperties> & fsp) const
{
  auto & gas = fsp[_gas_phase_number];

  // Calculate all of the gas phase properties, as these don't depend on saturation
  gasProperties(pressure, temperature, fsp);

  // The gas saturation in the two phase case
  gas.saturation = saturation(pressure, temperature, Xnacl, Z, fsp);

  // The liquid pressure and properties can now be calculated
  const DualReal liquid_pressure = pressure - _pc.capillaryPressure(1.0 - gas.saturation, qp);
  liquidProperties(liquid_pressure, temperature, Xnacl, fsp);
}

void
PorousFlowBrineH2::equilibriumMassFractions(const DualReal & pressure,
                                            const DualReal & temperature,
                                            const DualReal & Xnacl,
                                            DualReal & XH2,
                                            DualReal & YH2O) const
{
  // Mole fractions at equilibrium
  DualReal xH2, yH2O;
  equilibriumMoleFractions(pressure, temperature, Xnacl, xH2, yH2O);

  // The mass fraction of H2O in gas (assume no salt in gas phase) and derivatives
  // wrt p, T, and X
  YH2O = yH2O * _MH2O / (yH2O * _MH2O + (1.0 - yH2O) * _MH2);

  // NaCl molality (mol/kg)
  const DualReal mnacl = Xnacl / (1.0 - Xnacl) / _Mnacl;

  // The molality of H2 in 1kg of H2O
  const DualReal mH2 = xH2 * (2.0 * mnacl + _invMH2O) / (1.0 - xH2);
  // The mass fraction of H2 in brine is then
  const DualReal denominator = (1.0 + mnacl * _Mnacl + mH2 * _MH2);
  XH2 = mH2 * _MH2 / denominator;
}

void
PorousFlowBrineH2::equilibriumMoleFractions(const DualReal & pressure,
                                            const DualReal & temperature,
                                            const DualReal & Xnacl,
                                            DualReal & xH2,
                                            DualReal & yH2O) const
{
  // To calculate equilibrium mole fractions, need NaCl molality
  const DualReal mnacl = Xnacl / (1.0 - Xnacl) / _Mnacl;

  // Molality of H2 in liquid phase
  const DualReal mH2 = equilibriumMolality(pressure, temperature, mnacl);

  // Mole fraction of H2 in liquid phase
  xH2 = mH2 / (_invMH2O + mH2);

  // Mole fraction of H2O in gas phase
  yH2O = this->yH2O(pressure, temperature);
}

Real
PorousFlowBrineH2::totalMassFraction(
    Real pressure, Real temperature, Real Xnacl, Real saturation, unsigned int qp) const
{
  // Check whether the input pressure and temperature are within the region of validity
  checkVariables(pressure, temperature);

  // As we do not require derivatives, we can simply ignore their initialisation
  const DualReal p = pressure;
  const DualReal T = temperature;
  const DualReal X = Xnacl;

  // FluidStateProperties data structure
  std::vector<FluidStateProperties> fsp(_num_phases, FluidStateProperties(_num_components));
  FluidStateProperties & liquid = fsp[_aqueous_phase_number];
  FluidStateProperties & gas = fsp[_gas_phase_number];

  // Calculate equilibrium mass fractions in the two-phase state
  DualReal XH2, YH2O;
  equilibriumMassFractions(p, T, X, XH2, YH2O);

  // Save the mass fractions in the FluidStateMassFractions object
  const DualReal YH2 = 1.0 - YH2O;
  liquid.mass_fraction[_aqueous_fluid_component] = 1.0 - XH2;
  liquid.mass_fraction[_gas_fluid_component] = XH2;
  gas.mass_fraction[_aqueous_fluid_component] = YH2O;
  gas.mass_fraction[_gas_fluid_component] = YH2;

  // Gas properties
  gasProperties(pressure, temperature, fsp);

  // Liquid properties
  const DualReal liquid_saturation = 1.0 - saturation;
  const DualReal liquid_pressure = p - _pc.capillaryPressure(liquid_saturation, qp);
  liquidProperties(liquid_pressure, T, X, fsp);

  // The total mass fraction of ncg (z) can now be calculated
  const DualReal Z = (saturation * gas.density * YH2 + liquid_saturation * liquid.density * XH2) /
                     (saturation * gas.density + liquid_saturation * liquid.density);

  return Z.value();
}

DualReal
PorousFlowBrineH2::henryConstant(const DualReal & temperature) const
{
  // Henry's constant for dissolution in water
  const DualReal lnKh_H2O = _a[0] * Utility::pow<2>(temperature) + _a[1] * temperature + _a[2] +
                            _a[3] / temperature + _a[4] / Utility::pow<2>(temperature);

  // Note: ln(Kh) provided in Li et al
  return std::exp(lnKh_H2O);
}

DualReal
PorousFlowBrineH2::molarVolume(const DualReal & pressure, const DualReal & temperature) const
{
  // Pressure in MPa
  const DualReal p = pressure * 1.0e-6;
  const DualReal V =
      51.1904 - 0.208062 * temperature + 3.4427e-4 * Utility::pow<2>(temperature) - 0.022 * p;

  // V is in cm^3/mol, return m^3/mol
  return V * 1.0e-6;
}

DualReal
PorousFlowBrineH2::molarVolumeCH4(const DualReal & /*pressure*/, const DualReal & temperature) const
{
  const DualReal V = 36.75 + 0.1276 * temperature - 8.012e-4 * Utility::pow<2>(temperature) +
                     1.284e-6 * Utility::pow<3>(temperature);

  // V is in cm^3/mol, return m^3/mol
  return V * 1.0e-6;
}

DualReal
PorousFlowBrineH2::poyntingFactor(const DualReal & pressure, const DualReal & temperature) const
{
  // Pressure in MPa
  const DualReal p = pressure * 1.0e-6;

  const DualReal PF = _b[0] * p / temperature + _b[1] * p + _b[2] * p * temperature +
                      _b[3] * Utility::pow<2>(p) / temperature;

  return PF;
}

DualReal
PorousFlowBrineH2::activityCoefficient(const DualReal & temperature, const DualReal & mnacl) const
{
  // Note: ln(gamma) provided in Li et al
  return std::exp((0.64485 - 0.00142 * temperature) * mnacl);
}

DualReal
PorousFlowBrineH2::yH2O(const DualReal & pressure, const DualReal & temperature) const
{
  mooseAssert(pressure > 0.0, "PorousFlowBrineH2::yH2O(): pressure must be greater than zero");

  return _H2O_fp.vaporPressure(temperature) / pressure;
}

DualReal
PorousFlowBrineH2::fugacityCoefficient(const DualReal & pressure,
                                       const DualReal & temperature) const
{
  return _h2_fp.psi_from_p_T(pressure, temperature);
}

DualReal
PorousFlowBrineH2::fugacityCoefficientH2CH4(const DualReal & pressure,
                                            const DualReal & temperature,
                                            const DualReal & yH2,
                                            const DualReal & yCH4) const
{
  const DualReal Bmixt = yH2 * yH2 * binaryInteractionCoeffH2(temperature) +
                         2.0 * yH2 * yCH4 * binaryInteractionCoeffH2CH4(temperature) +
                         yCH4 * yCH4 * binaryInteractionCoeffCH4(temperature);

  const DualReal lnpsi = (2.0 * yH2 * binaryInteractionCoeffH2(temperature) +
                          2.0 * yCH4 * binaryInteractionCoeffH2CH4(temperature) - Bmixt) *
                         pressure / (_R * temperature);

  return std::exp(lnpsi);
}

DualReal
PorousFlowBrineH2::fugacityCoefficientH2N2(const DualReal & pressure,
                                           const DualReal & temperature,
                                           const DualReal & yH2,
                                           const DualReal & yN2) const
{
  const DualReal Bmixt = yH2 * yH2 * binaryInteractionCoeffH2(temperature) +
                         2.0 * yH2 * yN2 * binaryInteractionCoeffH2N2(temperature) +
                         yN2 * yN2 * binaryInteractionCoeffN2(temperature);

  const DualReal lnpsi = (2.0 * yH2 * binaryInteractionCoeffH2(temperature) +
                          2.0 * yN2 * binaryInteractionCoeffH2N2(temperature) - Bmixt) *
                         pressure / (_R * temperature);

  return std::exp(lnpsi);
}

DualReal
PorousFlowBrineH2::equilibriumMolality(const DualReal & pressure,
                                       const DualReal & temperature,
                                       const DualReal & mnacl) const
{
  const DualReal lnyH2 = std::log(1.0 - yH2O(pressure, temperature));
  const DualReal lnphi = std::log(fugacityCoefficient(pressure, temperature));
  const DualReal lnKh = std::log(henryConstant(temperature));
  const DualReal PF = poyntingFactor(pressure, temperature);
  const DualReal lngamma = std::log(activityCoefficient(temperature, mnacl));

  const DualReal lnmH2 =
      lnyH2 + std::log(pressure * 1.0e-6) + lnphi - lnKh - PF - lngamma + _ln_invMH2O;

  return std::exp(lnmH2);
}

DualReal
PorousFlowBrineH2::equilibriumMolalityH2N2(const DualReal & pressure,
                                           const DualReal & temperature,
                                           const DualReal & mnacl,
                                           const DualReal & yH2,
                                           const DualReal & yN2) const
{
  const DualReal lnphi = std::log(fugacityCoefficientH2N2(pressure, temperature, yH2, yN2));
  const DualReal lnKh = std::log(henryConstant(temperature));
  const DualReal PF = poyntingFactor(pressure, temperature);
  const DualReal lngamma = std::log(activityCoefficient(temperature, mnacl));

  const DualReal lnmH2 =
      std::log(yH2) + std::log(pressure * 1.0e-6) + lnphi - lnKh - PF - lngamma + _ln_invMH2O;

  return std::exp(lnmH2);
}

DualReal
PorousFlowBrineH2::equilibriumMolalityH2CH4(const DualReal & pressure,
                                            const DualReal & temperature,
                                            const DualReal & mnacl,
                                            const DualReal & yH2,
                                            const DualReal & yCH4) const
{
  const DualReal lnphi = std::log(fugacityCoefficientH2CH4(pressure, temperature, yH2, yCH4));
  const DualReal lnKh = std::log(henryConstant(temperature));
  const DualReal PF = poyntingFactor(pressure, temperature);
  const DualReal lngamma = std::log(activityCoefficient(temperature, mnacl));

  const DualReal lnmH2 =
      std::log(yH2) + std::log(pressure * 1.0e-6) + lnphi - lnKh - PF - lngamma + _ln_invMH2O;

  return std::exp(lnmH2);
}

DualReal
PorousFlowBrineH2::binaryInteractionCoeff(const DualReal & temperature,
                                          const std::array<Real, 4> & a) const
{
  const DualReal B = a[0] + a[1] * temperature + a[2] * Utility::pow<2>(temperature) +
                     a[3] * Utility::pow<3>(temperature);

  // Return m^3/mol rather than cm^3/mol
  return B * 1.0e-6;
}

DualReal
PorousFlowBrineH2::binaryInteractionCoeffH2(const DualReal & temperature) const
{
  return binaryInteractionCoeff(temperature, _a_H2);
}

DualReal
PorousFlowBrineH2::binaryInteractionCoeffN2(const DualReal & temperature) const
{
  return binaryInteractionCoeff(temperature, _a_N2);
}

DualReal
PorousFlowBrineH2::binaryInteractionCoeffCH4(const DualReal & temperature) const
{
  return binaryInteractionCoeff(temperature, _a_CH4);
}

DualReal
PorousFlowBrineH2::binaryInteractionCoeffH2N2(const DualReal & temperature) const
{
  return binaryInteractionCoeff(temperature, _a_H2_N2);
}

DualReal
PorousFlowBrineH2::binaryInteractionCoeffH2CH4(const DualReal & temperature) const
{
  return binaryInteractionCoeff(temperature, _a_H2_CH4);
}

void
PorousFlowBrineH2::checkVariables(Real pressure, Real temperature) const
{
  // The calculation of mass fractions is valid from 0C <= T <= 100C, and
  // pressure less than 60 MPa
  if (temperature < 273.15 || temperature > 373.15)
    mooseException(name() + ": temperature " + Moose::stringify(temperature) +
                   " is outside range 273.15 K <= T <= 373.15 K");

  if (pressure > 60.0e7)
    mooseException(name() + ": pressure " + Moose::stringify(pressure) +
                   " must be less than 60 MPa");
}
