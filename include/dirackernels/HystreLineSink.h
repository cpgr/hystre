//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowLineGeometry.h"
#include "PorousFlowSumQuantity.h"
#include "PorousFlowDictator.h"

/**
 * Approximates a line sink by a sequence of Dirac Points
 */
class HystreLineSink : public PorousFlowLineGeometry
{
public:
  static InputParameters validParams();

  HystreLineSink(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  /// Add Dirac Points to the borehole
  virtual void addPoints() override;

  const PorousFlowDictator & _dictator;
  PorousFlowSumQuantity & _total_outflow_mass;
  const Real _flux;
  const unsigned int _phase;
  const unsigned int _fluid_component;
  const bool _multiply_by_mass_frac;
  const ADMaterialProperty<std::vector<std::vector<Real>>> & _mass_fractions;
};
