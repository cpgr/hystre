/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#include "FVRelPermBC.h"

registerMooseObject("hystreApp", FVRelPermBC);

InputParameters
FVRelPermBC::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<Real>("w_coeff", "The Brooks-Corey exponent of the wetting phase");
  params.addRequiredParam<Real>("nw_coeff", "The Brooks-Corey exponent of the non-wetting phase");
  params.addRangeCheckedParam<Real>(
      "swirr", 0, "swirr >= 0 & swirr < 1", "The irreducible saturation of the wetting phase");
  params.addRangeCheckedParam<Real>("snwirr",
                                    0,
                                    "snwirr >= 0 & snwirr < 1",
                                    "The irreducible saturation of the non-wetting phase");
  params.addParam<Real>(
      "krw_end", 1, "The endpoint relative permeability the wetting phase (default is 1)");
  params.addParam<Real>(
      "krnw_end", 1, "The endpoint relative permeability the non-wetting phase (default is 1)");
  return params;
}

FVRelPermBC::FVRelPermBC(const InputParameters & parameters)
  : Material(parameters),
    _relperm(declareADProperty<std::vector<Real>>("relperm")),
    _saturation(getADMaterialProperty<std::vector<Real>>("saturation")),
    _w_coeff(getParam<Real>("w_coeff")),
    _nw_coeff(getParam<Real>("nw_coeff")),
    _krw_end(getParam<Real>("krw_end")),
    _krnw_end(getParam<Real>("krnw_end")),
    _swirr(getParam<Real>("swirr")),
    _snwirr(getParam<Real>("snwirr"))
{
}

void
FVRelPermBC::computeQpProperties()
{
  // Hardcode two phases for now
  _relperm[_qp].resize(2);

  const ADReal s = (_saturation[_qp][0] - _swirr) / (1.0 - _swirr - _snwirr);

  if (MooseUtils::absoluteFuzzyGreaterEqual(s.value(), 1.0))
  {
    _relperm[_qp][0] = _krw_end;
    _relperm[_qp][1] = 0.0;
  }

  else if (MooseUtils::absoluteFuzzyLessEqual(s.value(), 0.0)) // We are at sw == swirr
  {
    _relperm[_qp][0] = 0.0;
    _relperm[_qp][1] = _krnw_end;
  }
  else
  {
    _relperm[_qp][0] = _krw_end * std::pow(s, _w_coeff);
    _relperm[_qp][1] = _krnw_end * std::pow(1.0 - s, _nw_coeff);
  }
}
