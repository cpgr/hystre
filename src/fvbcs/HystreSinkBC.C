/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#include "HystreSinkBC.h"
#include "PorousFlowDictator.h"
#include "MooseUtils.h"

registerADMooseObject("hystreApp", HystreSinkBC);

InputParameters
HystreSinkBC::validParams()
{
  InputParameters params = FVFluxBC::validParams();
  params.addRequiredParam<UserObjectName>("PorousFlowDictator",
                                          "The PorousFlowDictator UserObject");
  RealVectorValue g(0, 0, -9.81);
  params.addParam<RealVectorValue>("gravity", g, "Gravity vector. Defaults to (0, 0, -9.81)");
  params.addParam<unsigned int>("phase", 0, "The fluid phase");
  params.addParam<unsigned int>("fluid_component", 0, "The fluid component");
  params.addParam<Real>("aquifer_pressure", 0, "Pressure in the boundary aquifer");
  params.addRangeCheckedParam<Real>(
      "distance", 1, "distance > 0", "Distance to the infinite aquifer");
  return params;
}

HystreSinkBC::HystreSinkBC(const InputParameters & params)
  : FVFluxBC(params),
    _density(getADMaterialProperty<std::vector<Real>>("density")),
    _viscosity(getADMaterialProperty<std::vector<Real>>("viscosity")),
    _relperm(getADMaterialProperty<std::vector<Real>>("relperm")),
    _mass_fractions(getADMaterialProperty<std::vector<std::vector<Real>>>("mass_fractions")),
    _permeability(getMaterialProperty<RealTensorValue>("permeability")),
    _pressure(getADMaterialProperty<std::vector<Real>>("pressure")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _phase(getParam<unsigned int>("phase")),
    _fluid_component(getParam<unsigned int>("fluid_component")),
    _paquifer(getParam<Real>("aquifer_pressure")),
    _dist(getParam<Real>("distance"))
{
}

ADReal
HystreSinkBC::computeQpResidual()
{
  ADReal flux = 0.0;

  auto p = _pressure[_qp][_phase];

  const auto gradp = (p - _paquifer) / _dist;

  const auto mobility = _mass_fractions[_qp][_phase][_fluid_component] * _relperm[_qp][_phase] *
                        _permeability[_qp] * _density[_qp][_phase] / _viscosity[_qp][_phase] *
                        _normal * _normal;

  const auto pressure_grad = gradp + _density[_qp][_phase] * _gravity * _normal;

  flux = mobility * pressure_grad;

  return flux;
}
