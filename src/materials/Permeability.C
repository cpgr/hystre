/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#include "Permeability.h"

registerMooseObject("hystreApp", Permeability);

InputParameters
Permeability::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("perm_xx", "The xx component of the permeability tensor");
  params.addCoupledVar("perm_yy", 0.0, "The yy component of the permeability tensor");
  params.addCoupledVar("perm_zz", 0.0, "The zz component of the permeability tensor");
  params.addClassDescription("Permeability");
  return params;
}

Permeability::Permeability(const InputParameters & parameters)
  : Material(parameters),
    _permeability(declareProperty<RealTensorValue>("permeability")),
    _perm_xx(coupledValue("perm_xx")),
    _perm_yy(parameters.isParamSetByUser("perm_yy") ? coupledValue("perm_yy") : _perm_xx),
    _perm_zz(parameters.isParamSetByUser("perm_zz") ? coupledValue("perm_zz") : _perm_xx)
{
}

void
Permeability::computeQpProperties()
{
  RealTensorValue permeability(
      _perm_xx[_qp], 0.0, 0.0, 0.0, _perm_yy[_qp], 0.0, 0.0, 0.0, _perm_zz[_qp]);

  _permeability[_qp] = permeability;
}
