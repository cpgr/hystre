/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#pragma once

#include "FVTimeKernel.h"

class PorousFlowDictator;
class FVPorousFlowMassTimeDerivative : public FVTimeKernel
{
public:
  static InputParameters validParams();
  FVPorousFlowMassTimeDerivative(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  const PorousFlowDictator & _dictator;
  const unsigned int _num_phases;
  const unsigned int _fluid_component;

  const MaterialProperty<Real> & _porosity;
  const ADMaterialProperty<std::vector<Real>> & _density;
  const ADMaterialProperty<std::vector<Real>> & _saturation;
  const ADMaterialProperty<std::vector<std::vector<Real>>> & _mass_fractions;
  const MaterialProperty<std::vector<Real>> & _density_old;
  const MaterialProperty<std::vector<Real>> & _saturation_old;
  const MaterialProperty<std::vector<std::vector<Real>>> & _mass_fractions_old;
};
