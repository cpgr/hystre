/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#pragma once

#include "Material.h"

class PorousFlowDictator;

class FVDiffusivityConstant : public Material
{
public:
  static InputParameters validParams();
  FVDiffusivityConstant(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

private:
  const PorousFlowDictator & _dictator;

  ADMaterialProperty<std::vector<Real>> & _tortuosity;
  ADMaterialProperty<std::vector<std::vector<Real>>> & _diffusion_coeff;
  const std::vector<Real> _input_tortuosity;
  const std::vector<Real> _input_diffusion_coeff;
  const unsigned int _num_phases;
  const unsigned int _num_components;
};
