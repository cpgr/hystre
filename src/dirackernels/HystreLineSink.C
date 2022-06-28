//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HystreLineSink.h"
#include "libmesh/utility.h"

registerADMooseObject("hystreApp", HystreLineSink);

InputParameters
HystreLineSink::validParams()
{
  InputParameters params = PorousFlowLineGeometry::validParams();
  params.addRequiredParam<UserObjectName>(
      "SumQuantityUO",
      "User Object of type=PorousFlowSumQuantity in which to place the total "
      "outflow from the line sink for each time step.");
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addParam<Real>("flux", 0, "Mass flux");
  params.addParam<unsigned int>("fluid_phase", 0, "The fluid phase for the line sink");
  params.addParam<unsigned int>(
      "mass_fraction_component", 0, "The index corresponding to a fluid component");
  params.addParam<bool>("multiply_by_mass_frac", false, "Multiply flux by mass fraction");
  params.addParam<bool>("multiply_by_relperm", false, "Multiply flux by relative permeability");
  params.addClassDescription("Approximates a line sink in the mesh by a sequence of weighted Dirac "
                             "points whose positions are read from a file");
  return params;
}

HystreLineSink::HystreLineSink(const InputParameters & parameters)
  : PorousFlowLineGeometry(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _total_outflow_mass(
        const_cast<PorousFlowSumQuantity &>(getUserObject<PorousFlowSumQuantity>("SumQuantityUO"))),
    _flux(getParam<Real>("flux")),
    _phase(getParam<unsigned int>("fluid_phase")),
    _fluid_component(getParam<unsigned int>("mass_fraction_component")),
    _multiply_by_mass_frac(getParam<bool>("multiply_by_mass_frac")),
    _mass_fractions(getADMaterialProperty<std::vector<std::vector<Real>>>("mass_fractions")),
    _multiply_by_relperm(getParam<bool>("multiply_by_relperm")),
    _relperm(getADMaterialProperty<std::vector<Real>>("relperm"))
{
  // zero the outflow mass
  _total_outflow_mass.zero();
}

void
HystreLineSink::addPoints()
{
  // This function gets called just before the DiracKernel is evaluated
  // so this is a handy place to zero this out.
  _total_outflow_mass.zero();

  PorousFlowLineGeometry::addPoints();
}

Real
HystreLineSink::computeQpResidual()
{
  // Get the ID we initially assigned to this point
  const unsigned current_dirac_ptid = currentPointCachedID();

  Real outflow = 0.0;

  if (current_dirac_ptid > 0)
    // contribution from half-segment "behind" this point (must have >1 point for
    // current_dirac_ptid>0)
    outflow += _half_seg_len[current_dirac_ptid - 1];

  if (current_dirac_ptid + 1 < _zs.size() || _zs.size() == 1)
    // contribution from half-segment "ahead of" this point, or we only have one point
    outflow += _half_seg_len[current_dirac_ptid];

  outflow *= _flux;

  if (outflow == 0.0)
    return 0.0;

  if (_multiply_by_mass_frac)
    outflow *= _mass_fractions[_qp][_phase][_fluid_component].value();

  if (_multiply_by_relperm)
    outflow *= _relperm[_qp][_phase].value();

  _total_outflow_mass.add(outflow * _dt);

  return outflow;
}
