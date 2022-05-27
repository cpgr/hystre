//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADPorousFlowPropertyAux.h"

registerMooseObject("hystreApp", ADPorousFlowPropertyAux);

InputParameters
ADPorousFlowPropertyAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  MooseEnum property_enum("pressure saturation temperature density viscosity mass_fraction relperm "
                          "capillary_pressure enthalpy internal_energy porosity permeability ");
  params.addRequiredParam<MooseEnum>(
      "property", property_enum, "The fluid property that this auxillary kernel is to calculate");
  params.addParam<unsigned int>("phase", 0, "The index of the phase this auxillary kernel acts on");
  params.addParam<unsigned int>(
      "liquid_phase", 0, "The index of the liquid phase (used for capillary pressure)");
  params.addParam<unsigned int>(
      "gas_phase", 1, "The index of the gas phase (used for capillary pressure)");
  params.addParam<unsigned int>(
      "fluid_component", 0, "The index of the fluid component this auxillary kernel acts on");
  params.addRangeCheckedParam<unsigned int>(
      "row", 0, "row>=0 & row<=2", "Row of permeability tensor to output");
  params.addRangeCheckedParam<unsigned int>(
      "column", 0, "column>=0 & column<=2", "Column of permeability tensor to output");
  params.addClassDescription("AuxKernel to provide access to properties evaluated at quadpoints. "
                             "Note that elemental AuxVariables must be used, so that these "
                             "properties are integrated over each element.");
  return params;
}

ADPorousFlowPropertyAux::ADPorousFlowPropertyAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _property_enum(getParam<MooseEnum>("property").getEnum<PropertyEnum>()),
    _phase(getParam<unsigned int>("phase")),
    _liquid_phase(getParam<unsigned int>("liquid_phase")),
    _gas_phase(getParam<unsigned int>("gas_phase")),
    _fluid_component(getParam<unsigned int>("fluid_component")),
    _k_row(getParam<unsigned int>("row")),
    _k_col(getParam<unsigned int>("column"))
{
  // Check that the phase and fluid_component are valid
  if (_phase >= _dictator.numPhases())
    paramError("phase",
               "Phase number entered is greater than the number of phases specified in the "
               "Dictator. Remember that indexing starts at 0");

  if (_fluid_component >= _dictator.numComponents())
    paramError("fluid_component",
               "Fluid component number entered is greater than the number of fluid components "
               "specified in the Dictator. Remember that indexing starts at 0");

  // Check the parameters used to calculate capillary pressure
  if (_property_enum == PropertyEnum::CAPILLARY_PRESSURE)
  {
    if (_liquid_phase >= _dictator.numPhases())
      paramError(
          "liquid_phase",
          "Liquid phase number entered is greater than the number of phases specified in the "
          "Dictator. Remember that indexing starts at 0");

    if (_gas_phase >= _dictator.numPhases())
      paramError("gas_phase",
                 "Gas phase number entered is greater than the number of phases specified in the "
                 "Dictator. Remember that indexing starts at 0");

    if (_liquid_phase == _gas_phase)
      paramError("liquid_phase", "Liquid phase number entered cannot be equal to gas_phase");
  }

  // Only get material properties required by this instance of the AuxKernel
  switch (_property_enum)
  {
    case PropertyEnum::PRESSURE:
      _pressure = &getADMaterialProperty<std::vector<Real>>("pressure");
      break;

    case PropertyEnum::SATURATION:
      _saturation = &getADMaterialProperty<std::vector<Real>>("saturation");
      break;

    case PropertyEnum::TEMPERATURE:
      _temperature = &getADMaterialProperty<Real>("temperature");
      break;

    case PropertyEnum::DENSITY:
      _fluid_density = &getADMaterialProperty<std::vector<Real>>("density");
      break;

    case PropertyEnum::VISCOSITY:
      _fluid_viscosity = &getADMaterialProperty<std::vector<Real>>("viscosity");
      break;

    case PropertyEnum::MASS_FRACTION:
      _mass_fractions = &getADMaterialProperty<std::vector<std::vector<Real>>>("mass_fractions");
      break;

    case PropertyEnum::RELPERM:
      _relative_permeability = &getADMaterialProperty<std::vector<Real>>("relperm");
      break;

    case PropertyEnum::CAPILLARY_PRESSURE:
      _pressure = &getADMaterialProperty<std::vector<Real>>("pressure");
      break;

    case PropertyEnum::ENTHALPY:
      _enthalpy = &getADMaterialProperty<std::vector<Real>>("enthalpy");
      break;

    case PropertyEnum::INTERNAL_ENERGY:
      _internal_energy = &getADMaterialProperty<std::vector<Real>>("internal_energy");
      break;

    case PropertyEnum::POROSITY:
      _porosity = &getMaterialProperty<Real>("porosity");
      break;

    case PropertyEnum::PERMEABILITY:
      _permeability = &getMaterialProperty<RealTensorValue>("permeability");
      break;
  }
}

Real
ADPorousFlowPropertyAux::computeValue()
{
  Real property = 0.0;

  switch (_property_enum)
  {
    case PropertyEnum::PRESSURE:
      property = (*_pressure)[_qp][_phase].value();
      break;

    case PropertyEnum::SATURATION:
      property = (*_saturation)[_qp][_phase].value();
      break;

    case PropertyEnum::TEMPERATURE:
      property = (*_temperature)[_qp].value();
      break;

    case PropertyEnum::DENSITY:
      property = (*_fluid_density)[_qp][_phase].value();
      break;

    case PropertyEnum::VISCOSITY:
      property = (*_fluid_viscosity)[_qp][_phase].value();
      break;

    case PropertyEnum::MASS_FRACTION:
      property = (*_mass_fractions)[_qp][_phase][_fluid_component].value();
      break;

    case PropertyEnum::RELPERM:
      property = (*_relative_permeability)[_qp][_phase].value();
      break;

    case PropertyEnum::CAPILLARY_PRESSURE:
      property = ((*_pressure)[_qp][_gas_phase] - (*_pressure)[_qp][_liquid_phase]).value();
      break;

    case PropertyEnum::ENTHALPY:
      property = (*_enthalpy)[_qp][_phase].value();
      break;

    case PropertyEnum::INTERNAL_ENERGY:
      property = (*_internal_energy)[_qp][_phase].value();
      break;

    case PropertyEnum::POROSITY:
      property = (*_porosity)[_qp];
      break;

    case PropertyEnum::PERMEABILITY:
      property = (*_permeability)[_qp](_k_row, _k_col);
      break;
  }

  return property;
}
