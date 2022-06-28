[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 1
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 1
  []
  [aquifer]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    bottom_left = '0 0 0'
    top_right = '10 1 0'
    block_id = 1
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [pp]
    family = MONOMIAL
    order = CONSTANT
    fv = true
    initial_condition = 1E7
  []
[]

[FVBCs]
  [left]
    type = HystreSinkBC
    boundary = left
    variable = pp
    fluid_component = 0
    phase = 0
    aquifer_pressure = 1e7
    distance = 1e2
  []
[]

[DiracKernels]
  [produce]
    type = ConstantPointSource
    point = '90 0.5 0'
    variable = pp
    value = -1e-2
  []
[]

[AuxVariables]
  [Ych4]
    initial_condition = 1
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
  [sgas]
    initial_condition = 0
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
  [xnacl]
    initial_condition = 0
    family = MONOMIAL
    order = CONSTANT
    fv = true
  []
[]

[FVKernels]
  [mass]
    type = FVPorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  []
  [flux]
    type = FVPorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pp
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp'
    number_fluid_phases = 2
    number_fluid_components = 3
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
    # gas_saturation = sgas
    gas_saturation = 0
    # YCH4 = Ych4
    YCH4 = 0
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
    block = 0
  []
  [permeability1]
    type = Permeability
    perm_xx = 1e-12
    block = 1
  []
  [relperm]
    type = FVRelPermBC
    nw_coeff = 2
    w_coeff = 2
  []
[]

[Postprocessors]
  [fluid_mass0]
    type = FVPorousFlowFluidMass
    execute_on = 'initial timestep_end'
    fluid_component = 0
    phases = '0 1'
  []
  [p00]
    type = PointValue
    variable = pp
    point = '15 0 0'
    execute_on = 'initial timestep_end'
  []
  [p01]
    type = PointValue
    variable = pp
    point = '90 0 0'
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [usual]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason'
    # petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -ksp_max_it'
    # petsc_options_value = 'bcgs bjacobi 1E-10 1E-10 10000 30'
  []
[]

[Executioner]
  type = Transient
  end_time = 2e3
  dt = 100
  solve_type = NEWTON
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = false
  csv = true
  execute_on = timestep_end
[]
