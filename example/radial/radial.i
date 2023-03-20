# H2 injection in a brine aquifer
# 2D radial model

################################################################################
# Injection schedule
# Inject 10,000 kgs H2 over 5 days, then store for 30 days, then produce for 10 days
inject_mass = 1.E4
inject_time = 5
store_time = 30
produce_mass = 1.E4
produce_time = 10
end_time = '${fparse inject_time + store_time + produce_time}'

start_injection = 0
end_injection = '${fparse start_injection + inject_time}'

start_production = '${fparse end_injection + store_time}'
end_production = '${fparse start_production + produce_time}'

# Time in seconds (used in inputs below)
inject_time_s = '${fparse inject_time * 86400}'
produce_time_s = '${fparse produce_time * 86400}'
start_injection_s = '${fparse start_injection * 86400}'
end_injection_s = '${fparse end_injection * 86400}'
start_production_s = '${fparse start_production * 86400}'
end_production_s = '${fparse end_production * 86400}'
end_time_s = '${fparse end_time * 86400}'

################################################################################
# Mesh details (geometry in RZ coordinates)
# borehole radius (m)
bh_r = 0.1
# model radius (m)
max_r = 1000
# aquifer base
min_y = 0
# aquifer top (m)
max_y = 30
# number of elements in radial direction
num_r = 100
# mesh bias in radial direction
bias_r = 1.04
# number of elements in vertical direction
num_y = 30
# mesh bias in vertical direction in aquifer top
bias_y = 0.97
# depth (m) - for hydrostatic pressure initial condition
depth = 1000

################################################################################
# Aquifer properties
# Insitu temperature at centre of aquifer (degC)
temp = 50
# Vertical geothermal gradient (K/m).  A positive number means temperature increases towards the centre of the earth
geothermal_gradient = 0
# Gravity
gravity = -9.81
# Porosity
porosity = 0.25
# Permeability (m^2)
hor_perm = 1E-13
ver_perm = 1E-14
# Salinity (mass fraction)
xnacl = 0.1

################################################################################
# Input file
[Mesh]
  [aquifer]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${num_r}
    xmin = ${bh_r}
    xmax = ${max_r}
    bias_x = ${bias_r}
    bias_y = ${bias_y}
    ny = ${num_y}
    ymin = ${min_y}
    ymax = ${max_y}
  []
[]

[Problem]
  coord_type = RZ
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 ${gravity} 0'
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pgas zi'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    alpha = 1e-3
    m = 0.4
    sat_lr = 0.1
  []
  [fs]
    type = PorousFlowBrineH2
    brine_fp = brine
    h2_fp = h2tab
    capillary_pressure = pc
  []
[]

[Variables]
  [pgas]
  []
  [zi]
    initial_condition = 0
    scaling = 1e4
  []
[]

[ICs]
  [pgas]
    type = FunctionIC
    variable = pgas
    function = insitu_pressure
  []
  [temperature]
    type = FunctionIC
    variable = temperature
    function = insitu_temperature
  []
[]

[BCs]
  [outer_boundary_porepressure]
    type = FunctionDirichletBC
    preset = true
    variable = pgas
    function = insitu_pressure
    boundary = right
  []
  [inject_h2]
    type = PorousFlowSink
    variable = zi
    boundary = left
    flux_function = injection_rate
    fluid_phase = 1
    enable = false
  []
  [produce_h2]
    type = PorousFlowSink
    variable = zi
    boundary = left
    flux_function = production_rate
    fluid_phase = 1
    enable = false
    use_relperm = true
    mass_fraction_component = 1
  []
  [produce_brine]
    type = PorousFlowSink
    variable = pgas
    boundary = left
    flux_function = production_rate
    fluid_phase = 0
    enable = false
    use_relperm = true
    mass_fraction_component = 0
  []
[]

[Controls]
  [inject]
    type = TimePeriod
    enable_objects = 'BCs::inject_h2'
    start_time = ${start_injection_s}
    end_time = ${end_injection_s}
    set_sync_times = true
    execute_on = 'initial timestep_begin'
    implicit = false
  []
  [produce]
    type = TimePeriod
    enable_objects = 'BCs::produce_h2 BCs::produce_brine'
    start_time = ${start_production_s}
    end_time = ${end_production_s}
    set_sync_times = true
    execute_on = 'initial timestep_begin'
    implicit = false
  []
[]

[Functions]
  [insitu_pressure]
    type = ParsedFunction
    value = '(y-${depth}) * 1100 * ${gravity} + 1E5'
  []
  [insitu_temperature]
    type = ParsedFunction
    value = '${temp} + (${depth} - y) * ${geothermal_gradient}'
  []
  [injection_rate]
    type = ParsedFunction
    value = '-${inject_mass}/(2 * pi * ${bh_r} * (${max_y} - ${min_y}) * ${inject_time_s})'
  []
  [production_rate]
    type = ParsedFunction
    value = '${produce_mass}/(2 * pi * ${bh_r} * (${max_y} - ${min_y}) * ${produce_time_s})'
  []
[]

[Kernels]
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pgas
  []
  [flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pgas
  []
  [mass1]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = zi
  []
  [flux1]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = zi
  []
[]

[AuxVariables]
  [temperature]
  []
  [liquid_density]
    order = CONSTANT
    family = MONOMIAL
  []
  [gas_density]
    order = CONSTANT
    family = MONOMIAL
  []
  [pressure_liquid]
    order = CONSTANT
    family = MONOMIAL
  []
  [gas_saturation]
    order = CONSTANT
    family = MONOMIAL
  []
  [xH2]
    order = CONSTANT
    family = MONOMIAL
  []
  [yH20]
    order = CONSTANT
    family = MONOMIAL
  []
  [xnacl]
    initial_condition = ${xnacl}
  []
[]

[AuxKernels]
  [liquid_density]
    type = PorousFlowPropertyAux
    variable = liquid_density
    property = density
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [gas_density]
    type = PorousFlowPropertyAux
    variable = gas_density
    property = density
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [pressure_liquid]
    type = PorousFlowPropertyAux
    variable = pressure_liquid
    property = pressure
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [gas_saturation]
    type = PorousFlowPropertyAux
    variable = gas_saturation
    property = saturation
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [xH2]
    type = PorousFlowPropertyAux
    variable = xH2
    property = mass_fraction
    phase = 0
    fluid_component = 1
    execute_on = 'initial timestep_end'
  []
  [yH20]
    type = PorousFlowPropertyAux
    variable = yH20
    property = mass_fraction
    phase = 1
    fluid_component = 0
    execute_on = 'initial timestep_end'
  []
[]

[FluidProperties]
  [h2]
    type = HystreHydrogenFluidProperties
  []
  [h2tab]
    type = HystreTabulatedFluidProperties
    fp = h2
    save_file = false
  []
  [water]
    type = Water97FluidProperties
  []
  [watertab]
    type = TabulatedFluidProperties
    fp = water
    save_file = false
  []
  [brine]
    type = BrineFluidProperties
    water_fp = watertab
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temperature
  []
  [brineh2]
    type = PorousFlowFluidState
    gas_porepressure = pgas
    z = zi
    temperature_unit = Celsius
    xnacl = xnacl
    capillary_pressure = pc
    fluid_state = fs
  []
  [porosity]
    type = PorousFlowPorosityConst
    porosity = ${porosity}
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '${hor_perm} 0 0   0 ${ver_perm} 0   0 0 ${ver_perm}'
  []
  [relperm_liquid]
    type = PorousFlowRelativePermeabilityVG
    m = 0.5
    phase = 0
    s_res = 0.1
    sum_s_res = 0.15
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 1
    s_res = 0.05
    sum_s_res = 0.15
  []
[]

[Postprocessors]
  [massH2]
    type = PorousFlowFluidMass
    fluid_component = 1
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
    petsc_options_value = 'asm lu NONZERO'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  end_time = ${end_time_s}
  l_max_its = 200
  nl_rel_tol = 1e-06
  nl_abs_tol = 1e-7
  nl_max_its = 10
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e2
    growth_factor = 2
    optimal_iterations = 8
    iteration_window = 3
  []
  dtmax = 8.64e4
[]

[Outputs]
  exodus = true
  csv = true
  perf_graph = true
  print_linear_residuals = false
[]
