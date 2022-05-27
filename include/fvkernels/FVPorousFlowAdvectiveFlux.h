/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#pragma once

#include "FVFluxKernel.h"

class PorousFlowDictator;

class FVPorousFlowAdvectiveFlux : public FVFluxKernel
{
public:
  static InputParameters validParams();
  FVPorousFlowAdvectiveFlux(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;

  const PorousFlowDictator & _dictator;
  const unsigned int _num_phases;
  const unsigned int _fluid_component;

  const ADMaterialProperty<std::vector<Real>> & _density;
  const ADMaterialProperty<std::vector<Real>> & _density_neighbor;

  const ADMaterialProperty<std::vector<Real>> & _viscosity;
  const ADMaterialProperty<std::vector<Real>> & _viscosity_neighbor;

  const ADMaterialProperty<std::vector<Real>> & _relperm;
  const ADMaterialProperty<std::vector<Real>> & _relperm_neighbor;

  const ADMaterialProperty<std::vector<std::vector<Real>>> & _mass_fractions;
  const ADMaterialProperty<std::vector<std::vector<Real>>> & _mass_fractions_neighbor;

  const MaterialProperty<RealTensorValue> & _permeability;
  const MaterialProperty<RealTensorValue> & _permeability_neighbor;

  const ADMaterialProperty<std::vector<Real>> & _pressure;
  const ADMaterialProperty<std::vector<Real>> & _pressure_neighbor;

  const RealVectorValue & _gravity;
};
