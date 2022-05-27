/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#include "FVDiffusivityConstant.h"
#include "PorousFlowDictator.h"

registerMooseObject("hystreApp", FVDiffusivityConstant);

InputParameters
FVDiffusivityConstant::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredParam<std::vector<Real>>(
      "diffusion_coeff",
      "List of diffusion coefficients.  Order is i) component 0 in phase 0; ii) "
      "component 1 in phase 0 ...; component 0 in phase 1; ... component k in "
      "phase n (m^2/s");
  params.addRequiredParam<std::vector<Real>>(
      "tortuosity", "List of tortuosities. Order is i) phase 0; ii) phase 1; etc");
  return params;
}

FVDiffusivityConstant::FVDiffusivityConstant(const InputParameters & parameters)
  : Material(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _tortuosity(declareADProperty<std::vector<Real>>("tortuosity")),
    _diffusion_coeff(declareADProperty<std::vector<std::vector<Real>>>("diffusion_coeff")),
    _input_tortuosity(getParam<std::vector<Real>>("tortuosity")),
    _input_diffusion_coeff(getParam<std::vector<Real>>("diffusion_coeff")),
    _num_phases(_dictator.numPhases()),
    _num_components(_dictator.numComponents())
{
  // Check that the correct number of input parameters have been entered
  if (_input_diffusion_coeff.size() != _num_phases * _num_components)
    paramError("diffusion_coeff",
               "The number of diffusion coefficients entered is not equal to the number of phases "
               "multiplied by the number of fluid components");

  if (_input_tortuosity.size() != _num_phases)
    paramError("tortuosity",
               "The number of tortuosity values entered is not equal to the number of phases "
               "specified in the Dictator");
}

void
FVDiffusivityConstant::computeQpProperties()
{
  _diffusion_coeff[_qp].resize(_num_phases);
  _tortuosity[_qp].resize(_num_phases);

  for (unsigned int p = 0; p < _num_phases; ++p)
  {
    _diffusion_coeff[_qp][p].resize(_num_components);
    _tortuosity[_qp][p] = _input_tortuosity[p];

    for (unsigned int c = 0; c < _num_components; ++c)
      _diffusion_coeff[_qp][p][c] = _input_diffusion_coeff[p + c];
  }
}
