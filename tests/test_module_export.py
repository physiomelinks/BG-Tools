#===============================================================================
#
#  A bond graph as a libcuflynx / PhLynx module version (bg2cellml/module_export.py)
#
#===============================================================================

import asyncio
import json
import os
from pathlib import Path
import subprocess
import sys

import lxml.etree as etree
import pytest

from bg2cellml.bondgraph.framework import get_framework
from bg2cellml.module_export import export_module

from .cellml_eval import evaluate

#===============================================================================

MODELS = Path(__file__).parent / 'models'
CELLML = '{http://www.cellml.org/cellml/1.1#}'

#===============================================================================

def export(path: Path, out: Path, module_type: str) -> tuple[str, list[Path]]:
    async def convert():
        framework = await get_framework()
        model = framework.make_bondgraph_model(path.resolve().as_uri(), path.read_text())
        assert not model.has_issues, [issue.reason for issue in model.issues]
        cellml = model.make_cellml_model().to_xml()
        return cellml, export_module(model, cellml, out, module_type, 'bondgraph_v01', source=path.name)
    return asyncio.run(convert())

def parameters(out: Path) -> dict[str, float]:
    data = json.loads((out / 'parameterisations' / 'default' / 'default_parameters.json').read_text())
    assert data.pop('schema_version') == 1
    assert all(set(p) == {'value', 'units', 'data_reference', 'sourced'} for p in data.values())
    return {name: p['value'] for name, p in data.items()}

def with_parameters(out: Path, module_type: str) -> str:
    """The exported CellML with its parameters' values put back, as libcuflynx sets them."""
    root = etree.parse(str(out / f'{module_type}_bondgraph_v01_modules.cellml')).getroot()
    values = parameters(out)
    for variable in root.iter(f'{CELLML}variable'):
        init = variable.get('initial_value')
        if init in values:
            variable.set('initial_value', str(values[init]))
        elif variable.get('name') in values:
            variable.set('initial_value', str(values[variable.get('name')]))
    return etree.tostring(root).decode()

#===============================================================================

def test_a_module_version_in_libcuflynx_and_phlynx_formats(tmp_path):
    '''The files of a module version, its config entry in PhLynx's keys, every value a parameter, and one port
    for the two copies of chemostat A.'''
    _, files = export(MODELS / 'chemostat_copies.ttl', tmp_path, 'copies')
    assert sorted(f.relative_to(tmp_path).as_posix() for f in files) == [
        'copies_bondgraph_v01_modules.cellml', 'copies_bondgraph_v01_modules_config.json',
        'copies_bondgraph_v01_references.bib', 'copies_bondgraph_v01_tests.yaml',
        'copies_bondgraph_v01_units.cellml', 'copies_bondgraph_v01_verification_config.json',
        'parameterisations/default/default_parameters.json']
    (entry,) = json.loads((tmp_path / 'copies_bondgraph_v01_modules_config.json').read_text())
    assert {k: entry[k] for k in ('module_type', 'module_subtype', 'module_format', 'component_file',
                                  'component_type', 'default_parameterisation')} == {
        'module_type': 'copies', 'module_subtype': 'bondgraph_v01', 'module_format': 'cellml',
        'component_file': 'copies_bondgraph_v01_modules.cellml', 'component_type': 'copies',
        'default_parameterisation': 'default'}
    assert entry['general_ports'] == [{'port_type': 'chemostat_A_c', 'variables': ['u_A_c', 'v_A_c']}]

    root = etree.parse(str(tmp_path / 'copies_bondgraph_v01_modules.cellml')).getroot()
    (component,) = root.findall(f'{CELLML}component')
    assert component.get('name') == 'copies'
    variables = {v.get('name'): v for v in component.findall(f'{CELLML}variable')}
    assert {name for name, *_ in entry['variables_and_units']} == set(variables) - {'t'}
    assert all(v.get('public_interface') in ('in', 'out') for v in variables.values())
    kinds = {name: kind for name, _, _, kind in entry['variables_and_units']}
    values = parameters(tmp_path)
    # every value is a parameter: none is left in the CellML but a state's reference to its <state>_init
    assert {name for name, kind in kinds.items() if kind in ('constant', 'boundary_condition')} == set(values)
    assert kinds['u_A_c'] == 'boundary_condition' and values['u_A_c'] == 1000
    assert variables['q_B'].get('initial_value') == 'q_B_init' and values['q_B_init'] == 0.2
    assert entry['licence'] == 'CC0-1.0' and entry['creator'] == []
    assert entry['required_citations'] == ['argus2026circulatory']
    assert '@misc{argus2026circulatory,' in (tmp_path / 'copies_bondgraph_v01_references.bib').read_text()
    assert all(v.get('initial_value') in (None, 'q_B_init') for v in variables.values())
    assert root.find(f'{CELLML}units') is None
    assert etree.parse(str(tmp_path / 'copies_bondgraph_v01_units.cellml')).getroot().find(f'{CELLML}units') is not None

def test_the_port_flow_is_the_modules_uptake_of_the_species(tmp_path):
    '''Both copies of A take the port's potential, and the port's flow is what the module takes from its
    surroundings: B's rate of gain, here.'''
    export(MODELS / 'chemostat_copies.ttl', tmp_path, 'copies')
    values, rates = evaluate(with_parameters(tmp_path, 'copies'))
    copies = [name for name in values if name.startswith('u_Se_A')]
    assert len(copies) == 2 and all(values[c] == pytest.approx(1000.0) for c in copies)
    assert values['v_A_c'] == pytest.approx(rates['q_B'])
    assert values['v_A_c'] > 0

def test_a_product_chemostats_port_flow_is_negative_uptake(tmp_path):
    '''A chemostat a reaction produces into: the module gives the species to its surroundings, so the port's
    flow is minus the reaction's.'''
    export(MODELS / 'chemostat_product.ttl', tmp_path, 'product')
    values, rates = evaluate(with_parameters(tmp_path, 'product'))
    assert values['v_B_c'] == pytest.approx(rates['q_A_A_c'])
    assert values['v_B_c'] < 0

SIGNAL_MODELS = ['vessel_inertance', 'angII_vessel', 'baroreceptor', 'renin_secretion', 'aldosterone_transport',
                 'kedem_katchalsky']

def test_signal_port_blocks_are_entrance_and_exit_ports(tmp_path):
    '''A signal input port block is an entrance port of its bgf:portType, its value a parameter; a signal output
    port block an exit port that may feed several modules.'''
    export(MODELS / 'baroreceptor.ttl', tmp_path, 'baroreceptor')
    (entry,) = json.loads((tmp_path / 'baroreceptor_bondgraph_v01_modules_config.json').read_text())
    assert entry['entrance_ports'] == [{'port_type': 'arterial_pressure_port', 'variables': ['P_sys']}]
    assert entry['exit_ports'] == [{'port_type': 'baroreceptor_firing_rate_port', 'variables': ['f_ab'],
                                    'multi_port': 'True'}]
    kinds = {name: (units, kind) for name, units, _, kind in entry['variables_and_units']}
    assert kinds['P_sys'] == ('kPa', 'boundary_condition') and kinds['f_ab'] == ('per_second', 'variable')
    values = parameters(tmp_path)
    assert values['P_sys'] == 13 and values['P_baro_init'] == 12 and values['tau_pb'] == 2.076
    root = etree.parse(str(tmp_path / 'baroreceptor_bondgraph_v01_modules.cellml')).getroot()
    interface = {v.get('name'): v.get('public_interface') for v in root.iter(f'{CELLML}variable')}
    assert interface['P_sys'] == 'in' and interface['f_ab'] == 'out'

def _schema(name: str):
    path = Path(LIBCUFLYNX_SRC or '') / 'libcuflynx' / 'schemas' / name
    if not LIBCUFLYNX_SRC or not path.exists():
        pytest.skip('set LIBCUFLYNX_SRC to check against libcuflynx schemas')
    return json.loads(path.read_text())

@pytest.mark.parametrize('model', SIGNAL_MODELS)
def test_the_module_version_matches_libcuflynxs_schemas(model, tmp_path):
    '''The config entry and default parameterisation of each example validate against libcuflynx's schemas, and
    the required citations are keys of the version's references.'''
    jsonschema = pytest.importorskip('jsonschema')
    export(MODELS / f'{model}.ttl', tmp_path, model)
    config = json.loads((tmp_path / f'{model}_bondgraph_v01_modules_config.json').read_text())
    jsonschema.validate(config, _schema('module_config.schema.json'))
    jsonschema.validate(json.loads((tmp_path / 'parameterisations' / 'default' / 'default_parameters.json').read_text()),
                        _schema('parameters.schema.json'))
    bib = (tmp_path / f'{model}_bondgraph_v01_references.bib').read_text()
    assert all(f'{{{key},' in bib for key in config[0]['required_citations'])
    names = {name for name, *_ in config[0]['variables_and_units']}
    for port in config[0]['entrance_ports'] + config[0]['exit_ports'] + config[0]['general_ports']:
        assert set(port['variables']) <= names

#===============================================================================

LIBCUFLYNX_SRC = os.environ.get('LIBCUFLYNX_SRC')          # a circulatory_autogen checkout's src/
LIBCUFLYNX_PYTHON = os.environ.get('LIBCUFLYNX_PYTHON')    # a Python with its dependencies (and myokit)

@pytest.mark.skipif(not (LIBCUFLYNX_SRC and LIBCUFLYNX_PYTHON),
                    reason='set LIBCUFLYNX_SRC and LIBCUFLYNX_PYTHON to build the module with libcuflynx')
@pytest.mark.parametrize('model', ['chemostat_copies', 'chemostat_product', 'electrogenic_reaction', 'chemostat_named',
                                   *SIGNAL_MODELS])
def test_libcuflynx_builds_the_module_and_runs_it_as_bg2cellml(model, tmp_path):
    '''libcuflynx builds a one-instance model of the module from its default parameterisation alone, and Myokit
    runs it, from off rest, exactly as it runs bg2cellml's own CellML.'''
    module_dir = tmp_path / 'module'
    cellml, _ = export(MODELS / f'{model}.ttl', module_dir, model)
    original = tmp_path / f'{model}.cellml'
    original.write_text(cellml)
    env = dict(os.environ, PYTHONPATH=LIBCUFLYNX_SRC)
    result = subprocess.run([LIBCUFLYNX_PYTHON, str(Path(__file__).parent / 'libcuflynx_run.py'), str(module_dir),
                             model, 'bondgraph_v01', str(original), str(tmp_path / 'work')],
                            capture_output=True, text=True, env=env, timeout=600)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
    report = json.loads(result.stdout.strip().splitlines()[-1])
    assert report['moved'] > 1e-3          # the states change over the run
    assert report['difference'] <= 1e-12, report

#===============================================================================
#===============================================================================
