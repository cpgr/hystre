/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#include "FVPorousFlowOutflowBC.h"
#include "PorousFlowDictator.h"
#include "MooseUtils.h"

registerADMooseObject("hystreApp", FVPorousFlowOutflowBC);

InputParameters
FVPorousFlowOutflowBC::validParams()
{
  InputParameters params = FVFluxBC::validParams();
  params.addRequiredParam<UserObjectName>("PorousFlowDictator",
                                          "The PorousFlowDictator UserObject");
  RealVectorValue g(0, 0, -9.81);
  params.addParam<RealVectorValue>("gravity", g, "Gravity vector. Defaults to (0, 0, -9.81)");
  params.addParam<unsigned int>("fluid_component", 0, "The fluid component");
  return params;
}

FVPorousFlowOutflowBC::FVPorousFlowOutflowBC(const InputParameters & params)
  : FVFluxBC(params),
    _density(getADMaterialProperty<std::vector<Real>>("density")),
    _density_neighbor(getNeighborADMaterialProperty<std::vector<Real>>("density")),
    _viscosity(getADMaterialProperty<std::vector<Real>>("viscosity")),
    _viscosity_neighbor(getNeighborADMaterialProperty<std::vector<Real>>("viscosity")),
    _relperm(getADMaterialProperty<std::vector<Real>>("relperm")),
    _relperm_neighbor(getNeighborADMaterialProperty<std::vector<Real>>("relperm")),
    _mass_fractions(getADMaterialProperty<std::vector<std::vector<Real>>>("mass_fractions")),
    _mass_fractions_neighbor(
        getNeighborADMaterialProperty<std::vector<std::vector<Real>>>("mass_fractions")),
    _permeability(getMaterialProperty<RealTensorValue>("permeability")),
    _permeability_neighbor(getNeighborMaterialProperty<RealTensorValue>("permeability")),
    _pressure(getADMaterialProperty<std::vector<Real>>("pressure")),
    _pressure_neighbor(getNeighborADMaterialProperty<std::vector<Real>>("pressure")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _fluid_component(getParam<unsigned int>("fluid_component"))
{
}

ADReal
FVPorousFlowOutflowBC::computeQpResidual()
{
  ADReal flux = 0.0;

  for (unsigned int p = 0; p < _num_phases; ++p)
  {
    ADReal p_elem = _pressure[_qp][p];
    ADReal p_neighbor = _pressure_neighbor[_qp][p];

    const ADRealVectorValue gradp =
        (p_elem - p_neighbor) * _face_info->eCF() / _face_info->dCFMag();

    const ADRealTensorValue mobility_element = _mass_fractions[_qp][p][_fluid_component] *
                                               _relperm[_qp][p] * _permeability[_qp] *
                                               _density[_qp][p] / _viscosity[_qp][p];

    const ADRealTensorValue mobility_neighbor =
        _mass_fractions_neighbor[_qp][p][_fluid_component] * _relperm_neighbor[_qp][p] *
        _permeability_neighbor[_qp] * _density_neighbor[_qp][p] / _viscosity_neighbor[_qp][p];

    const auto pressure_grad = gradp + _density[_qp][p] * _gravity;

    ADRealTensorValue mobility_upwind;
    interpolate(Moose::FV::InterpMethod::Upwind,
                mobility_upwind,
                mobility_element,
                mobility_neighbor,
                pressure_grad,
                *_face_info,
                true);

    flux += mobility_upwind * pressure_grad * _normal;
  }

  // Restrict flux to be positive (no inflow allowed)
  return (flux.value() > 0.0 ? flux : 0.0);
}
