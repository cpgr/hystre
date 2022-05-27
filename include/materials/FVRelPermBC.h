/*****************************************************************/
/*    FINCH - FINite volume Capillary Heterogeneity modelling    */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#pragma once

#include "Material.h"

class FVRelPermBC : public Material
{
public:
  static InputParameters validParams();
  FVRelPermBC(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  ADMaterialProperty<std::vector<Real>> & _relperm;
  const ADMaterialProperty<std::vector<Real>> & _saturation;

private:
  const Real _w_coeff;
  const Real _nw_coeff;
  const Real _krw_end;
  const Real _krnw_end;
  const Real _swirr;
  const Real _snwirr;
};
