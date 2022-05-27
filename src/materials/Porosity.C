/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#include "Porosity.h"

registerMooseObject("hystreApp", Porosity);

InputParameters
Porosity::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("porosity", "Porosity");
  params.addClassDescription("Porosity");
  return params;
}

Porosity::Porosity(const InputParameters & parameters)
  : Material(parameters),
    _porosity(declareProperty<Real>("porosity")),
    _porosity_value(coupledValue("porosity"))
{
}

void
Porosity::computeQpProperties()
{
  _porosity[_qp] = _porosity_value[_qp];
}
