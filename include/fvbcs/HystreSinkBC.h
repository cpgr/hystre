/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#pragma once

#include "FVFluxBC.h"

class PorousFlowDictator;

class HystreSinkBC : public FVFluxBC
{
public:
  static InputParameters validParams();
  HystreSinkBC(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty<std::vector<Real>> & _density;
  const ADMaterialProperty<std::vector<Real>> & _viscosity;
  const ADMaterialProperty<std::vector<Real>> & _relperm;
  const ADMaterialProperty<std::vector<std::vector<Real>>> & _mass_fractions;
  const MaterialProperty<RealTensorValue> & _permeability;
  const ADMaterialProperty<std::vector<Real>> & _pressure;
  const RealVectorValue & _gravity;
  const PorousFlowDictator & _dictator;
  const unsigned int _num_phases;
  const unsigned int _phase;
  const unsigned int _fluid_component;
  const Real _paquifer;
  const Real _dist;
};
