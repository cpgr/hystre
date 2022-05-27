/*****************************************************************/
/*           HYSTRE - HYdrogen STorage in REservoirs             */
/*                                                               */
/*           All contents are licensed under MIT/BSD             */
/*              See LICENSE for full restrictions                */
/*****************************************************************/

#pragma once

#include "Material.h"

class PorousFlowDictator;
class PorousFlowCapillaryPressure;
class BrineFluidProperties;
class SinglePhaseFluidProperties;

class FVBasicH2CH4Brine : public Material
{
public:
  static InputParameters validParams();
  FVBasicH2CH4Brine(const InputParameters & parameters);

protected:
  void thermophysicalProperties();
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Size material property vectors and initialise with zeros
  void setMaterialVectorSize() const;

  const PorousFlowDictator & _dictator;

  const ADVariableValue & _liquid_porepressure;
  const ADVariableValue & _gas_saturation;
  const ADVariableValue & _YCH4;
  const ADVariableValue & _temperature;
  const ADVariableValue & _Xnacl;

  /// Phase pressures
  ADMaterialProperty<std::vector<Real>> & _pressure;
  /// Phase saturations
  ADMaterialProperty<std::vector<Real>> & _saturation;
  /// Mass fractions
  ADMaterialProperty<std::vector<std::vector<Real>>> & _mass_frac;
  /// Density of each phase
  ADMaterialProperty<std::vector<Real>> & _fluid_density;
  /// Viscosity of each phase
  ADMaterialProperty<std::vector<Real>> & _fluid_viscosity;

private:
  const unsigned int _num_phases;
  const unsigned int _num_components;
  /// Conversion from degrees Celsius to degrees Kelvin
  const Real _T_c2k;
  /// Capillary pressure UserObject
  const PorousFlowCapillaryPressure & _pc;
  /// Fluid properties UserObject for water
  const BrineFluidProperties & _brine_fp;
  /// Fluid properties UserObject for the H2
  const SinglePhaseFluidProperties & _h2_fp;
  /// Fluid properties UserObject for the CH4
  const SinglePhaseFluidProperties & _ch4_fp;
};
