/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#pragma once

#include "ElementIntegralPostprocessor.h"

/**
 * Postprocessor produces the mass of a given fluid component in a region
 */
class PorousFlowDictator;

class FVPorousFlowFluidMass : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();

  FVPorousFlowFluidMass(const InputParameters & parameters);

  virtual Real computeQpIntegral() override;

  const PorousFlowDictator & _dictator;
  const unsigned int _fluid_component;
  std::vector<unsigned int> _phases;
  const Real _saturation_threshold;
  const MaterialProperty<Real> & _porosity;
  const ADMaterialProperty<std::vector<Real>> & _density;
  const ADMaterialProperty<std::vector<Real>> & _saturation;
  const ADMaterialProperty<std::vector<std::vector<Real>>> & _mass_fractions;
};
