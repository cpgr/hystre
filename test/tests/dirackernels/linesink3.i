[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmin = 0
  xmax = 4
  ymin = 0
  ymax = 2
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [pp]
    family = MONOMIAL
    order = CONSTANT
    fv = true
    initial_condition = 1E7
  []
  [sgas]
    initial_condition = 0.5
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
  [Ych4]
    initial_condition = 0.9
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
[]

[AuxVariables]
  [xnacl]
    initial_condition = 0
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
[]

[FVKernels]
  [mass0]
    type = FVPorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  []
  [mass1]
    type = FVPorousFlowMassTimeDerivative
    fluid_component = 1
    variable = sgas
  []
  [mass2]
    type = FVPorousFlowMassTimeDerivative
    fluid_component = 2
    variable = Ych4
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp sgas Ych4'
    number_fluid_phases = 2
    number_fluid_components = 3
  []
  [total_outflow_mass_h2]
    type = PorousFlowSumQuantity
  []
  [total_outflow_mass_ch4]
    type = PorousFlowSumQuantity
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    m = 0.5
    alpha = 1e-7
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
    liquid_porepressure = pp
    gas_saturation = sgas
    YCH4 = Ych4
    temperature_unit = Celsius
    xnacl = xnacl
    capillary_pressure = pc
    brine_fp = brine
    h2_fp = h2
    ch4_fp = ch4
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

[DiracKernels]
  [pls]
    type = HystreLineSink
    flux = 1
    fluid_phase = 1
    mass_fraction_component = 1
    point_file = well.bh
    SumQuantityUO = total_outflow_mass_h2
    variable = sgas
    line_length = 1
    multiply_by_mass_frac = true
  []
  [pls2]
    type = HystreLineSink
    flux = 1
    fluid_phase = 1
    mass_fraction_component = 2
    point_file = well.bh
    SumQuantityUO = total_outflow_mass_ch4
    variable = Ych4
    line_length = 1
    multiply_by_mass_frac = true
  []
[]

[Postprocessors]
  [pls_report_h2]
    type = PorousFlowPlotQuantity
    uo = total_outflow_mass_h2
  []
  [pls_report_ch4]
    type = PorousFlowPlotQuantity
    uo = total_outflow_mass_ch4
  []
  [fluid_mass0]
    type = FVPorousFlowFluidMass
    execute_on = 'initial TIMESTEP_END'
    fluid_component = 0
    phases = '0 1'
  []
  [fluid_mass1]
    type = FVPorousFlowFluidMass
    phases = '0 1'
    fluid_component = 1
    execute_on = 'initial TIMESTEP_END'
  []
  [fluid_mass2]
    type = FVPorousFlowFluidMass
    phases = '0 1'
    fluid_component = 2
    execute_on = 'initial TIMESTEP_END'
  []
[]

[Preconditioning]
  [usual]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  end_time = 5
  dt = 1
  solve_type = NEWTON
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = false
  csv = true
  execute_on = timestep_end
[]
