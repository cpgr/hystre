//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVBasicH2CH4Brine.h"
#include "PorousFlowDictator.h"
#include "PorousFlowCapillaryPressure.h"
#include "BrineFluidProperties.h"
#include "SinglePhaseFluidProperties.h"

registerMooseObject("hystreApp", FVBasicH2CH4Brine);

InputParameters
FVBasicH2CH4Brine::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredCoupledVar("liquid_porepressure", "Liquid pressure");
  params.addRequiredCoupledVar("gas_saturation", "Gas saturation");
  params.addRequiredCoupledVar("YCH4", "Mass fraction of CH4 in gas phase");
  params.addCoupledVar(
      "temperature", 293, "The fluid temperature (C or K, depending on temperature_unit)");
  params.addCoupledVar("xnacl", 0, "The salt mass fraction in the brine (kg/kg)");
  MooseEnum unit_choice("Kelvin=0 Celsius=1", "Kelvin");
  params.addParam<MooseEnum>(
      "temperature_unit", unit_choice, "The unit of the temperature variable");
  params.addRequiredParam<UserObjectName>("capillary_pressure",
                                          "Name of the Capillary Pressure UserObject");
  params.addRequiredParam<UserObjectName>("brine_fp", "The name of the user object for brine");
  params.addRequiredParam<UserObjectName>("h2_fp", "The name of the user object for H2");
  params.addRequiredParam<UserObjectName>("ch4_fp", "The name of the user object for CH4");
  params.addClassDescription("Class for mixed H2-CH4-brine with no dissolution");
  return params;
}

FVBasicH2CH4Brine::FVBasicH2CH4Brine(const InputParameters & parameters)
  : Material(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _liquid_porepressure(adCoupledValue("liquid_porepressure")),
    _gas_saturation(adCoupledValue("gas_saturation")),
    _YCH4(adCoupledValue("YCH4")),
    _temperature(adCoupledValue("temperature")),
    _Xnacl(adCoupledValue("xnacl")),
    _pressure(declareADProperty<std::vector<Real>>("pressure")),
    _saturation(declareADProperty<std::vector<Real>>("saturation")),
    _mass_frac(declareADProperty<std::vector<std::vector<Real>>>("mass_fractions")),
    _fluid_density(declareADProperty<std::vector<Real>>("density")),
    _fluid_viscosity(declareADProperty<std::vector<Real>>("viscosity")),
    _num_phases(_dictator.numPhases()),
    _num_components(_dictator.numComponents()),
    _T_c2k(getParam<MooseEnum>("temperature_unit") == 0 ? 0.0 : 273.15),
    _pc(getUserObject<PorousFlowCapillaryPressure>("capillary_pressure")),
    _brine_fp(getUserObject<BrineFluidProperties>("brine_fp")),
    _h2_fp(getUserObject<SinglePhaseFluidProperties>("h2_fp")),
    _ch4_fp(getUserObject<SinglePhaseFluidProperties>("ch4_fp"))
{
}

void
FVBasicH2CH4Brine::thermophysicalProperties()
{
  // The FluidProperty objects use temperature in K
  const ADReal Tk = _temperature[_qp] + _T_c2k;

  // Set the size of all vectors
  setMaterialVectorSize();

  // Liquid phase (phase 0)
  _saturation[_qp][0] = 1.0 - _gas_saturation[_qp];
  _pressure[_qp][0] = _liquid_porepressure[_qp];
  _fluid_density[_qp][0] = _brine_fp.rho_from_p_T_X(_liquid_porepressure[_qp], Tk, _Xnacl[_qp]);
  _fluid_viscosity[_qp][0] = _brine_fp.mu_from_p_T_X(_liquid_porepressure[_qp], Tk, _Xnacl[_qp]);
  _mass_frac[_qp][0][0] = 1.0; // All liquid is brine

  // Gas phase (phase 1)
  _saturation[_qp][1] = _gas_saturation[_qp];
  _pressure[_qp][1] = _liquid_porepressure[_qp] + _pc.capillaryPressure(1.0 - _gas_saturation[_qp]);

  _fluid_density[_qp][1] = (1.0 - _YCH4[_qp]) * _h2_fp.rho_from_p_T(_pressure[_qp][1], Tk) +
                           _YCH4[_qp] * _ch4_fp.rho_from_p_T(_pressure[_qp][1], Tk);

  _fluid_viscosity[_qp][1] = (1.0 - _YCH4[_qp]) * _h2_fp.mu_from_p_T(_pressure[_qp][1], Tk) +
                             _YCH4[_qp] * _ch4_fp.mu_from_p_T(_pressure[_qp][1], Tk);

  _mass_frac[_qp][1][1] = 1.0 - _YCH4[_qp]; // H2 in gas phase
  _mass_frac[_qp][1][2] = _YCH4[_qp];       // CH4 in gas phase
}

void
FVBasicH2CH4Brine::initQpStatefulProperties()
{
  thermophysicalProperties();
}

void
FVBasicH2CH4Brine::computeQpProperties()
{
  thermophysicalProperties();
}

void
FVBasicH2CH4Brine::setMaterialVectorSize() const
{
  _pressure[_qp].assign(_num_phases, 0.0);
  _saturation[_qp].assign(_num_phases, 0.0);
  _fluid_density[_qp].assign(_num_phases, 0.0);
  _fluid_viscosity[_qp].assign(_num_phases, 0.0);
  _mass_frac[_qp].resize(_num_phases);

  for (unsigned int ph = 0; ph < _num_phases; ++ph)
    _mass_frac[_qp][ph].resize(_num_components);
}
