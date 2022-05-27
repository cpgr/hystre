/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#include "FVPorousFlowDispersiveFlux.h"
#include "PorousFlowDictator.h"

registerADMooseObject("hystreApp", FVPorousFlowDispersiveFlux);

InputParameters
FVPorousFlowDispersiveFlux::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  RealVectorValue g(0, 0, -9.81);
  params.addParam<RealVectorValue>("gravity", g, "Gravity vector. Defaults to (0, 0, -9.81)");
  params.addRequiredParam<UserObjectName>("PorousFlowDictator",
                                          "The PorousFlowDictator UserObject");
  params.addParam<unsigned int>("fluid_component", 0, "The fluid component");
  return params;
}

FVPorousFlowDispersiveFlux::FVPorousFlowDispersiveFlux(const InputParameters & params)
  : FVFluxKernel(params),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _num_phases(_dictator.numPhases()),
    _fluid_component(getParam<unsigned int>("fluid_component")),
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
    _porosity(getMaterialProperty<Real>("porosity")),
    _porosity_neighbor(getNeighborMaterialProperty<Real>("porosity")),
    _tortuosity(getADMaterialProperty<std::vector<Real>>("tortuosity")),
    _tortuosity_neighbor(getNeighborADMaterialProperty<std::vector<Real>>("tortuosity")),
    _diffusion_coeff(getADMaterialProperty<std::vector<std::vector<Real>>>("diffusion_coeff")),
    _diffusion_coeff_neighbor(
        getNeighborADMaterialProperty<std::vector<std::vector<Real>>>("diffusion_coeff")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _identity_tensor(RankTwoTensor::initIdentity)
{
}

ADReal
FVPorousFlowDispersiveFlux::computeQpResidual()
{
  ADReal flux = 0.0;

  for (unsigned int p = 0; p < _num_phases; ++p)
  {
    ADReal X_elem = _mass_fractions[_qp][p][_fluid_component];
    ADReal X_neighbor = _mass_fractions_neighbor[_qp][p][_fluid_component];

    const ADRealVectorValue gradX =
        (X_elem - X_neighbor) * _face_info->eCF() / _face_info->dCFMag();

    ADReal coeff_ave;
    const auto coeff = _porosity[_qp] * _tortuosity[_qp][p] * _density[_qp][p] *
                       _diffusion_coeff[_qp][p][_fluid_component];

    const auto coeff_neighbor = _porosity_neighbor[_qp] * _tortuosity_neighbor[_qp][p] *
                                _density_neighbor[_qp][p] *
                                _diffusion_coeff_neighbor[_qp][p][_fluid_component];

    interpolate(
        Moose::FV::InterpMethod::Average, coeff_ave, coeff, coeff_neighbor, *_face_info, true);

    flux += coeff_ave * _identity_tensor * gradX * _normal;
  }

  return flux;
}
