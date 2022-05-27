/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#include "FVPorousFlowFluidMass.h"
#include "PorousFlowDictator.h"

registerADMooseObject("hystreApp", FVPorousFlowFluidMass);

InputParameters
FVPorousFlowFluidMass::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names.");
  params.addParam<std::vector<unsigned int>>("phases",
                                             "The index of the fluid phase that this "
                                             "Postprocessor is restricted to.  Multiple "
                                             "indices can be entered");
  params.addRangeCheckedParam<Real>("saturation_threshold",
                                    1.0,
                                    "saturation_threshold >= 0 & saturation_threshold <= 1",
                                    "The saturation threshold below which the mass is calculated "
                                    "for a specific phase. Default is 1.0. Note: only one "
                                    "phase_index can be entered");
  params.addParam<unsigned int>("fluid_component", 0, "The fluid component whose mass to compute");
  return params;
}

FVPorousFlowFluidMass::FVPorousFlowFluidMass(const InputParameters & params)
  : ElementIntegralPostprocessor(params),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _fluid_component(getParam<unsigned int>("fluid_component")),
    _phases(getParam<std::vector<unsigned int>>("phases")),
    _saturation_threshold(getParam<Real>("saturation_threshold")),
    _porosity(getMaterialProperty<Real>("porosity")),
    _density(getADMaterialProperty<std::vector<Real>>("density")),
    _saturation(getADMaterialProperty<std::vector<Real>>("saturation")),
    _mass_fractions(getADMaterialProperty<std::vector<std::vector<Real>>>("mass_fractions"))
{
  const unsigned int num_phases = _dictator.numPhases();
  const unsigned int num_components = _dictator.numComponents();

  // Check that the number of components entered is not greater than the total number of components
  if (_fluid_component >= num_components)
    paramError(
        "fluid_component",
        "The Dictator proclaims that the number of components in this simulation is ",
        num_components,
        " whereas you have used a component index of ",
        _fluid_component,
        ". Remember that indexing starts at 0. The Dictator does not take such mistakes lightly.");

  // Check that the number of phases entered is not more than the total possible phases
  if (_phases.size() > num_phases)
    paramError("phase",
               "The Dictator decrees that the number of phases in this simulation is ",
               num_phases,
               " but you have entered ",
               _phases.size(),
               " phases.");

  // Using saturation_threshold only makes sense for a specific phase_index
  if (_saturation_threshold < 1.0 && _phases.size() != 1)
    paramError("saturation_threshold",
               "A single phase_index must be entered when prescribing a saturation_threshold");
}

Real
FVPorousFlowFluidMass::computeQpIntegral()
{
  Real mass = 0.0;

  for (auto p : _phases)
  {
    if (_saturation[_qp][p] <= _saturation_threshold)
      mass += _porosity[_qp] * _mass_fractions[_qp][p][_fluid_component].value() *
              _saturation[_qp][p].value() * _density[_qp][p].value();
  }

  // Clip to zero to avoid small negative masses due to precision
  return (mass > 0.0 ? mass : 0.0);
}
