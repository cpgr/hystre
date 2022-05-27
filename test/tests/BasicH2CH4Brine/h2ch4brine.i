# Tests correct calculation of properties in FVBasicH2CH4Brine

[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  temperature = 30
  gravity = '0 -9.81 0'
[]

[Variables]
  [pliq]
    initial_condition = 10e6
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
  [sgas]
    initial_condition = 0.5
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
  [ych4]
    initial_condition = 0.8
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
[]

[AuxVariables]
  [xnacl]
    initial_condition = 0.1
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
  [pressure_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [pressure_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [saturation_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [saturation_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [density_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [density_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [viscosity_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [viscosity_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [x0_water]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [x0_gas]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [x1_water]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [x1_gas]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [x2_water]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [x2_gas]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[AuxKernels]
  [pressure_water]
    type = ADPorousFlowPropertyAux
    variable = pressure_water
    property = pressure
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [pressure_gas]
    type = ADPorousFlowPropertyAux
    variable = pressure_gas
    property = pressure
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [saturation_water]
    type = ADPorousFlowPropertyAux
    variable = saturation_water
    property = saturation
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [saturation_gas]
    type = ADPorousFlowPropertyAux
    variable = saturation_gas
    property = saturation
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [density_water]
    type = ADPorousFlowPropertyAux
    variable = density_water
    property = density
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [density_gas]
    type = ADPorousFlowPropertyAux
    variable = density_gas
    property = density
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [viscosity_water]
    type = ADPorousFlowPropertyAux
    variable = viscosity_water
    property = viscosity
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [viscosity_gas]
    type = ADPorousFlowPropertyAux
    variable = viscosity_gas
    property = viscosity
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [x0_water]
    type = ADPorousFlowPropertyAux
    variable = x0_water
    property = mass_fraction
    phase = 0
    fluid_component = 0
    execute_on = 'initial timestep_end'
  []
  [x0_gas]
    type = ADPorousFlowPropertyAux
    variable = x0_gas
    property = mass_fraction
    phase = 1
    fluid_component = 0
    execute_on = 'initial timestep_end'
  []
  [x1_water]
    type = ADPorousFlowPropertyAux
    variable = x1_water
    property = mass_fraction
    phase = 0
    fluid_component = 1
    execute_on = 'initial timestep_end'
  []
  [x1_gas]
    type = ADPorousFlowPropertyAux
    variable = x1_gas
    property = mass_fraction
    phase = 1
    fluid_component = 1
    execute_on = 'initial timestep_end'
  []
  [x2_water]
    type = ADPorousFlowPropertyAux
    variable = x2_water
    property = mass_fraction
    phase = 0
    fluid_component = 2
    execute_on = 'initial timestep_end'
  []
  [x2_gas]
    type = ADPorousFlowPropertyAux
    variable = x2_gas
    property = mass_fraction
    phase = 1
    fluid_component = 2
    execute_on = 'initial timestep_end'
  []
[]

[FVKernels]
  [mass0]
    type = FVPorousFlowMassTimeDerivative
    variable = pliq
    fluid_component = 0
  []
  [flux0]
    type = FVPorousFlowAdvectiveFlux
    variable = pliq
    fluid_component = 0
  []
  [mass1]
    type = FVPorousFlowMassTimeDerivative
    variable = sgas
    fluid_component = 1
  []
  [flux1]
    type = FVPorousFlowAdvectiveFlux
    variable = pliq
    fluid_component = 1
  []
  [mass2]
    type = FVPorousFlowMassTimeDerivative
    variable = ych4
    fluid_component = 2
  []
  [flux2]
    type = FVPorousFlowAdvectiveFlux
    variable = pliq
    fluid_component = 2
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pliq sgas ych4'
    number_fluid_phases = 2
    number_fluid_components = 3
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  []
[]

[Modules]
  [FluidProperties]
    [h2]
      type = HydrogenFluidProperties
    []
    [ch4]
      type = MethaneFluidProperties
    []
    [brine]
      type = BrineFluidProperties
    []
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [brineh2ch4]
    type = FVBasicH2CH4Brine
    liquid_porepressure = pliq
    gas_saturation = sgas
    YCH4 = ych4
    temperature_unit = Celsius
    xnacl = xnacl
    capillary_pressure = pc
    brine_fp = brine
    h2_fp = h2
    ch4_fp = ch4
    pc = pc
  []
  [porosity]
    type = Porosity
    porosity = 0.2
  []
  [permeability]
    type = Permeability
    perm_xx = 1e-12
  []
  [relperm]
    type = FVRelPermBC
    nw_coeff = 2
    w_coeff = 2
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 1
  end_time = 1
  nl_abs_tol = 1e-12
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [density_water]
    type = ElementIntegralVariablePostprocessor
    variable = density_water
    execute_on = 'initial timestep_end'
  []
  [density_gas]
    type = ElementIntegralVariablePostprocessor
    variable = density_gas
    execute_on = 'initial timestep_end'
  []
  [viscosity_water]
    type = ElementIntegralVariablePostprocessor
    variable = viscosity_water
    execute_on = 'initial timestep_end'
  []
  [viscosity_gas]
    type = ElementIntegralVariablePostprocessor
    variable = viscosity_gas
    execute_on = 'initial timestep_end'
  []
  [x0_water]
    type = ElementIntegralVariablePostprocessor
    variable = x0_water
    execute_on = 'initial timestep_end'
  []
  [x1_water]
    type = ElementIntegralVariablePostprocessor
    variable = x1_water
    execute_on = 'initial timestep_end'
  []
  [x2_water]
    type = ElementIntegralVariablePostprocessor
    variable = x2_water
    execute_on = 'initial timestep_end'
  []
  [x0_gas]
    type = ElementIntegralVariablePostprocessor
    variable = x0_gas
    execute_on = 'initial timestep_end'
  []
  [x1_gas]
    type = ElementIntegralVariablePostprocessor
    variable = x1_gas
    execute_on = 'initial timestep_end'
  []
  [x2_gas]
    type = ElementIntegralVariablePostprocessor
    variable = x2_gas
    execute_on = 'initial timestep_end'
  []
  [sg]
    type = ElementIntegralVariablePostprocessor
    variable = saturation_gas
    execute_on = 'initial timestep_end'
  []
  [sw]
    type = ElementIntegralVariablePostprocessor
    variable = saturation_water
    execute_on = 'initial timestep_end'
  []
  [pwater]
    type = ElementIntegralVariablePostprocessor
    variable = pressure_water
    execute_on = 'initial timestep_end'
  []
  [pgas]
    type = ElementIntegralVariablePostprocessor
    variable = pressure_gas
    execute_on = 'initial timestep_end'
  []
  [x0mass]
    type = FVPorousFlowFluidMass
    fluid_component = 0
    phases = '0 1'
    execute_on = 'initial timestep_end'
  []
  [x1mass]
    type = FVPorousFlowFluidMass
    fluid_component = 1
    phases = '0 1'
    execute_on = 'initial timestep_end'
  []
  [x2mass]
    type = FVPorousFlowFluidMass
    fluid_component = 2
    phases = '0 1'
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  exodus = true
  execute_on = 'initial timestep_end'
  perf_graph = true
[]
