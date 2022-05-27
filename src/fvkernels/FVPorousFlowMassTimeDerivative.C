/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#include "FVPorousFlowMassTimeDerivative.h"
#include "PorousFlowDictator.h"

registerADMooseObject("hystreApp", FVPorousFlowMassTimeDerivative);

InputParameters
FVPorousFlowMassTimeDerivative::validParams()
{
  InputParameters params = FVTimeKernel::validParams();
  params.addRequiredParam<UserObjectName>("PorousFlowDictator",
                                          "The PorousFlowDictator UserObject");
  params.addParam<unsigned int>("fluid_component", 0, "The fluid component");
  return params;
}

FVPorousFlowMassTimeDerivative::FVPorousFlowMassTimeDerivative(const InputParameters & parameters)
  : FVTimeKernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _fluid_component(getParam<unsigned int>("fluid_component")),
    _porosity(getMaterialProperty<Real>("porosity")),
    _density(getADMaterialProperty<std::vector<Real>>("density")),
    _saturation(getADMaterialProperty<std::vector<Real>>("saturation")),
    _mass_fractions(getADMaterialProperty<std::vector<std::vector<Real>>>("mass_fractions")),
    _density_old(getMaterialPropertyOld<std::vector<Real>>("density")),
    _saturation_old(getMaterialPropertyOld<std::vector<Real>>("saturation")),
    _mass_fractions_old(getMaterialPropertyOld<std::vector<std::vector<Real>>>("mass_fractions"))
{
}

ADReal
FVPorousFlowMassTimeDerivative::computeQpResidual()
{
  ADReal mass = 0.0;
  Real mass_old = 0.0;

  for (unsigned p = 0; p < _num_phases; ++p)
  {
    mass += _density[_qp][p] * _saturation[_qp][p] * _mass_fractions[_qp][p][_fluid_component];
    mass_old += _density_old[_qp][p] * _saturation_old[_qp][p] *
                _mass_fractions_old[_qp][p][_fluid_component];
  }

  return _porosity[_qp] * (mass - mass_old) / _dt;
}
