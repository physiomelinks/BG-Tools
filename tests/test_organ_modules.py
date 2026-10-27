#===============================================================================
#
#  Organ modules: the Na/K-ATPase, arithmetic signal blocks, ports of several variables, a reference version's
#  form, causality, a model's unit system and values in other units
#
#===============================================================================

import json
from pathlib import Path

import lxml.etree as etree
import pytest

from .cellml_eval import cellml_from_bgrdf, evaluate
from .test_module_export import export, parameters, with_parameters

#===============================================================================

MODELS = Path(__file__).parent / 'models'
CELLML = '{http://www.cellml.org/cellml/1.1#}'

#===============================================================================
# The kidney of BVC-RAAS6: NKE pump, Kedem-Katchalsky transport, convective Na, a renal bed port of seven variables

def test_the_kidney_module_has_the_ports_of_its_neighbours(tmp_path):
    export(MODELS / 'kidney_organ.ttl', tmp_path, 'kidney')
    (entry,) = json.loads((tmp_path / 'kidney_bondgraph_v01_modules_config.json').read_text())
    entrances = {p['port_type']: p['variables'] for p in entry['entrance_ports']}
    exits = {p['port_type']: p['variables'] for p in entry['exit_ports']}
    # one port of seven elements, in their bgf:portOrder; a source's port of only its potential
    assert entrances['renal_bed_port'] == ['u_gl', 'u_ptc', 'q_gl', 'q_ptc', 'v_aff', 'v_eff', 'v_rv']
    assert entrances['bladder_pressure_port'] == ['u_bl']
    # a resistance's and a convective flux's flows in one port, a chemostat's flow made an exit
    assert exits['urine_port'] == ['v_ure_W', 'v_ure_Na']
    assert exits['plasma_Na_flux_port'] == ['v_kid_Na'] and exits['glomerular_filtration_port'] == ['v_gl_W']
    kinds = {name: kind for name, _, _, kind in entry['variables_and_units']}
    assert kinds['u_gl'] == kinds['v_rv'] == kinds['u_bl'] == 'boundary_condition'
    assert kinds['v_ure_Na'] == kinds['v_kid_Na'] == 'variable'
    values = parameters(tmp_path)
    assert values['RT'] == 2557                     # the pump's value of the domain constant
    assert values['sigma_pt'] == 0.7 and 'sigma_pt_2' not in values     # one constant, set in two elements


def test_the_kidney_module_computes_its_flows(tmp_path):
    export(MODELS / 'kidney_organ.ttl', tmp_path, 'kidney')
    values, rates = evaluate(with_parameters(tmp_path, 'kidney'))
    RT, F = 2557.0, 96485.0
    c_ptEpi_Na, c_ptEpi_K = 0.0015/0.1, 0.00725/0.1
    c_ptc_Na, c_gl_Na, c_pt_Na = 0.00416789/0.0297465, 0.000444596/0.00317569, 0.013026/0.199755
    scale = 1 + 2.1133e-08/(2.1133e-08 + 1.13e-09)
    n_i, n_o, k_o, k_i = c_ptEpi_Na**3, c_ptc_Na**3, 0.0045**2, c_ptEpi_K**2
    nke = 3*0.00210127*scale*(310646000*n_i*k_o - n_o*k_i*__import__('math').exp(2*F*-0.04/RT)) \
        / ((1 + n_i)*(1 + n_o)*(1 + k_o)*(1 + k_i))
    assert values['v_ptEpi_NKE'] == pytest.approx(nke)
    assert values['v_kid_Na'] == pytest.approx(c_ptc_Na*0.0179727 - 0.14*0.018)
    assert rates['q_gl_Na'] == pytest.approx(0.14*0.018 - c_gl_Na*(0.0161388 + values['v_gl_W']))
    assert values['v_ure_Na'] == pytest.approx(c_pt_Na*values['v_ure_W'])
    u_pt = 1*(0.199755 - 0.1)/(0.5 - 0.199755)**2
    assert values['v_ure_W'] == pytest.approx(6.00903e-06*max(u_pt + 3.92266 - 0.480526, 0))
    assert values['FF'] == pytest.approx(values['v_gl_W']/(0.55*0.018))
    assert values['u_ptEpi_osmotic'] == pytest.approx(RT*(2*0.0015 + 0.00725 + 0.01675)/0.1)


def test_a_reference_version_keeps_its_units_names_and_parameter_metadata(tmp_path):
    '''Exported over an existing version: its names of units (where the definitions agree), its parameters' metadata
    with the bond graph's values, its other config keys and its order.'''
    reference = tmp_path / 'reference'
    export(MODELS / 'kidney_organ.ttl', reference, 'kidney')
    (entry,) = json.loads((reference / 'kidney_bondgraph_v01_modules_config.json').read_text())
    entry['required_citations_uncertain'] = ['argus2026circulatory']
    entry['variables_and_units'] = [['u_gl', 'kPa', 'access', 'boundary_condition']] + entry['variables_and_units']
    (reference / 'kidney_bondgraph_v01_modules_config.json').write_text(json.dumps([entry]))
    (reference / 'kidney_bondgraph_v01_units.cellml').write_text(
        f'<model xmlns="{CELLML[1:-1]}" name="Units"><units name="kPa"><unit prefix="kilo" units="pascal"/></units>'
        '</model>')
    path = reference / 'parameterisations' / 'default' / 'default_parameters.json'
    data = json.loads(path.read_text())
    data['K_f'].update(data_reference='Thomas 2019', sourced='yes', value=1.0)
    path.write_text(json.dumps(data))

    from bg2cellml.module_export import export_module
    import asyncio
    from bg2cellml.bondgraph.framework import get_framework
    out, warnings = tmp_path / 'out', []
    async def convert():
        framework = await get_framework()
        source = MODELS / 'kidney_organ.ttl'
        model = framework.make_bondgraph_model(source.resolve().as_uri(), source.read_text())
        cellml = model.make_cellml_model().to_xml()
        export_module(model, cellml, out, 'kidney', 'bondgraph_v01', reference=reference, warn=warnings.append)
    asyncio.run(convert())
    (new,) = json.loads((out / 'kidney_bondgraph_v01_modules_config.json').read_text())
    assert new['required_citations_uncertain'] == ['argus2026circulatory']
    units = {name: u for name, u, _, _ in new['variables_and_units']}
    assert units['u_gl'] == 'kPa' and units['u_bl'] == 'kPa'      # its own name, and the one of the same definition
    assert new['variables_and_units'][0][0] == 'u_gl'
    assert 'name="kPa"' in (out / 'kidney_bondgraph_v01_units.cellml').read_text()
    k_f = json.loads((out / 'parameterisations' / 'default' / 'default_parameters.json').read_text())['K_f']
    assert k_f['data_reference'] == 'Thomas 2019' and k_f['sourced'] == 'yes' and k_f['value'] == 0.000254357
    assert any(w.startswith('K_f:') for w in warnings)

def test_dimensionless_units_keep_the_built_in_name(tmp_path):
    '''A library's alias of dimensionless does not rename the built-in units: without a reference, and over a
    reference version that names them dimensionless, as libcuflynx's parameter files do.'''
    units = tmp_path / 'library_units.cellml'
    units.write_text(f'<model xmlns="{CELLML[1:-1]}" name="Units">'
                     '<units name="dim"><unit units="dimensionless"/></units></model>')
    reference = tmp_path / 'reference'
    export(MODELS / 'kidney_organ.ttl', reference, 'kidney')
    (entry,) = json.loads((reference / 'kidney_bondgraph_v01_modules_config.json').read_text())
    built_in = {name for name, u, _, _ in entry['variables_and_units'] if u == 'dimensionless'}
    assert {'sigma_pt', 'FF'} <= built_in
    import asyncio
    from bg2cellml.bondgraph.framework import get_framework
    from bg2cellml.module_export import export_module
    async def convert(out, reference):
        framework = await get_framework()
        source = MODELS / 'kidney_organ.ttl'
        model = framework.make_bondgraph_model(source.resolve().as_uri(), source.read_text())
        cellml = model.make_cellml_model().to_xml()
        export_module(model, cellml, out, 'kidney', 'bondgraph_v01', reference=reference, units_files=[units])
    for out, ref in ((tmp_path / 'library', None), (tmp_path / 'over_reference', reference)):
        asyncio.run(convert(out, ref))
        (entry,) = json.loads((out / 'kidney_bondgraph_v01_modules_config.json').read_text())
        names = {name: u for name, u, _, _ in entry['variables_and_units']}
        assert {name for name in built_in if names[name] != 'dimensionless'} == set(), out.name
        cellml = (out / 'kidney_bondgraph_v01_modules.cellml').read_text()
        assert '"dim"' not in cellml and '"dim"' not in (out / 'kidney_bondgraph_v01_units.cellml').read_text()
        data = json.loads((out / 'parameterisations' / 'default' / 'default_parameters.json').read_text())
        assert data['sigma_pt']['units'] == 'dimensionless'

#===============================================================================
# A vessel in SI units, whose wall resistance gives the pressure of the node a flow feeds

def test_a_resistance_fed_by_a_flow_gives_the_nodes_pressure():
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'flow_fed_node.ttl'))
    u_C = (1.2e-3 - 1e-3)/1e-7 + 500.0              # Pa, from values in L and kPa
    assert values['u_ext'] == pytest.approx(500.0)
    assert values['u_node'] == pytest.approx(u_C + 2e-4/5e-7)
    assert rates['q_C'] == pytest.approx(2e-4 - 1e-4)


def test_a_models_unit_system_converts_its_values():
    cellml = cellml_from_bgrdf(MODELS / 'flow_fed_node.ttl')
    root = etree.fromstring(cellml.encode())
    variables = {v.get('name'): v for v in root.iter(f'{CELLML}variable')}
    assert variables['q_C'].get('units') == 'metre_cubed' and float(variables['q_C'].get('initial_value')) == 1.2e-3
    assert variables['v_in'].get('units') == 'metre_cubed_per_second'
    assert float(variables['q_half'].get('initial_value')) == 1e-3
    # names are written as they are, not as Greek letters
    assert '<ci>alpha</ci>' in cellml and 'α' not in cellml


def test_an_entrance_flow_sources_port_has_the_nodes_pressure(tmp_path):
    export(MODELS / 'flow_fed_node.ttl', tmp_path, 'vessel')
    (entry,) = json.loads((tmp_path / 'vessel_bondgraph_v01_modules_config.json').read_text())
    entrances = {p['port_type']: p['variables'] for p in entry['entrance_ports']}
    assert entrances['vessel_port'] == ['v_in', 'u_in']
    values, _ = evaluate(with_parameters(tmp_path, 'vessel'))
    assert values['u_in'] == pytest.approx(values['u_node'])
    assert 'RT' not in parameters(tmp_path)           # a constant no equation uses is left out

#===============================================================================
# A vein whose wall resistance gives the inlet node's pressure, and ports that share a variable

def test_a_junction_balance_written_for_a_state_is_solved_for_the_wall_flow():
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'vein_wall.ttl'))
    G = 1e-06/0.01
    v_wall = 2e-04 - 1e-04                          # the inflow less the inertance's flow (a state)
    assert values['v_wall'] == pytest.approx(v_wall)
    assert values['u_A'] == pytest.approx((1.3e-03 - 1e-03)/1e-06 + v_wall/G)
    assert rates['q'] == pytest.approx(v_wall)
    assert rates['v'] == pytest.approx((values['u_A'] - 2000 - 1e6*1e-04)/0.01)
    # the node's pressure is written with the compliance's, not by rearranging its law
    root = etree.fromstring(cellml_from_bgrdf(MODELS / 'vein_wall.ttl').encode())
    M = '{http://www.w3.org/1998/Math/MathML}'
    (u_A,) = [a for a in root.iter(f'{M}apply') if len(a) == 3 and a[0].tag == f'{M}eq'
              and a[1].tag == f'{M}ci' and a[1].text.strip() == 'u_A']
    assert sorted(ci.text.strip() for ci in u_A[2].iter(f'{M}ci')) == ['u_C', 'u_R']


def test_a_signal_output_port_of_its_input_shares_a_variable(tmp_path):
    export(MODELS / 'vein_wall.ttl', tmp_path, 'vein')
    (entry,) = json.loads((tmp_path / 'vein_bondgraph_v01_modules_config.json').read_text())
    exits = {p['port_type']: p for p in entry['exit_ports']}
    # the inertance's flow is in the vessel_port and, under its own name, the flow_port (multi_port True by default)
    assert exits['vessel_port']['variables'] == ['v', 'u_out']
    assert exits['flow_port'] == {'port_type': 'flow_port', 'variables': ['v'], 'multi_port': 'True'}
    # bgf:multiPort "False": no multi_port
    assert exits['volume_port'] == {'port_type': 'volume_port', 'variables': ['q']}
    names = {name for name, *_ in entry['variables_and_units']}
    assert not names & {'flow_out', 'volume_out'}
    cellml = (tmp_path / 'vein_bondgraph_v01_modules.cellml').read_text()
    assert 'flow_out' not in cellml and 'volume_out' not in cellml
    assert {'v_init', 'q_init'} <= set(parameters(tmp_path))


def test_a_parameter_of_the_states_initial_name_and_value_is_its_initial_value(tmp_path):
    export(MODELS / 'vein_wall.ttl', tmp_path, 'vein')
    cellml = (tmp_path / 'vein_bondgraph_v01_modules.cellml').read_text()
    # the Constant block q_init, which r_q also uses, is the volume's initial value: one parameter q_init
    assert 'initial_value="q_init"' in cellml and 'q_init_2' not in cellml
    values = parameters(tmp_path)
    assert values['q_init'] == 1.3e-03 and 'q_init_2' not in values
    values, _ = evaluate(with_parameters(tmp_path, 'vein'))
    assert values['r_q'] == pytest.approx(1.0)


def test_a_reference_version_names_its_states_initial_values(tmp_path):
    reference = tmp_path / 'reference'
    export(MODELS / 'vein_wall.ttl', reference, 'vein')
    path = reference / 'vein_bondgraph_v01_modules.cellml'
    path.write_text(path.read_text().replace('initial_value="q_init"', 'initial_value="q_C_init"'))
    import asyncio
    from bg2cellml.bondgraph.framework import get_framework
    from bg2cellml.module_export import export_module
    async def convert():
        framework = await get_framework()
        source = MODELS / 'vein_wall.ttl'
        model = framework.make_bondgraph_model(source.resolve().as_uri(), source.read_text())
        export_module(model, model.make_cellml_model().to_xml(), tmp_path / 'out', 'vein', 'bondgraph_v01',
                      reference=reference)
    asyncio.run(convert())
    values = parameters(tmp_path / 'out')
    assert values['q_C_init'] == 1.3e-03 and 'v_init' in values
    assert 'initial_value="q_C_init"' in (tmp_path / 'out' / 'vein_bondgraph_v01_modules.cellml').read_text()

#===============================================================================
# The gut of BVC-RAAS5: a one-way pump and fluxes that laws of signal blocks set

def test_a_one_way_pump_and_fluxes_set_by_signals():
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'gut_organ_v01.ttl'))
    # Na into the epithelium at kK (q_gi - q_giEpi), not below 0, through a flux the law sets
    assert values['v_gi_Na'] == pytest.approx(max(0.003*(0.001444 - 0.0015), 0))
    # the pump runs backwards at these concentrations, so the one-way pump passes nothing
    assert values['v_giEpi_NKE'] == 0.0
    assert rates['q_giEpi_Na'] == pytest.approx(values['v_gi_Na'] - values['v_giEpi_NKE'])
    # water out by the osmotic difference, one way
    RT = 2577.0
    du = RT*0.001*0.217157/1.06921 - RT*0.001444/0.289
    assert values['v_gi_W'] == pytest.approx(max(0.00023925*du, 0))


def test_a_librarys_units_files_name_the_units(tmp_path):
    '''Without a reference version, a library's units file still names the units whose definitions agree.'''
    units = tmp_path / 'library_units.cellml'
    units.write_text(f'<model xmlns="{CELLML[1:-1]}" name="Units">'
                     '<units name="J_per_m3"><unit units="joule"/><unit units="metre" exponent="-3"/></units>'
                     '<units name="m3_per_s"><unit units="metre" exponent="3"/><unit units="second" exponent="-1"/></units>'
                     '</model>')
    import asyncio
    from bg2cellml.bondgraph.framework import get_framework
    from bg2cellml.module_export import export_module
    async def convert():
        framework = await get_framework()
        source = MODELS / 'flow_fed_node.ttl'
        model = framework.make_bondgraph_model(source.resolve().as_uri(), source.read_text())
        cellml = model.make_cellml_model().to_xml()
        export_module(model, cellml, tmp_path / 'out', 'vessel', 'bondgraph_v01', units_files=[units],
                      component_type='vessel_type')
    asyncio.run(convert())
    (entry,) = json.loads((tmp_path / 'out' / 'vessel_bondgraph_v01_modules_config.json').read_text())
    kinds = {name: u for name, u, _, _ in entry['variables_and_units']}
    assert kinds['u_node'] == 'J_per_m3' and kinds['v_in'] == 'm3_per_s'
    assert entry['component_type'] == 'vessel_type' and entry['notes'].startswith('Exported by bg2cellml')
