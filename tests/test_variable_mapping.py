#===============================================================================
#
#  The variable mapping in libcuflynx's format (bg2cellml/variable_mapping.py)
#
#===============================================================================

import asyncio
import csv
import json
import os
from pathlib import Path
import subprocess

import pytest

from bg2cellml.bondgraph.framework import get_framework
from bg2cellml.module_export import export_module
from bg2cellml.variable_mapping import (BONDGRAPH_COLUMN, COLUMNS, mapping_path, symbol, variable_mapping,
                                        write_variable_mapping)

#===============================================================================

MODELS = Path(__file__).parent / 'models'
LIBCUFLYNX_SRC = os.environ.get('LIBCUFLYNX_SRC')          # a circulatory_autogen checkout's src/
LIBCUFLYNX_PYTHON = os.environ.get('LIBCUFLYNX_PYTHON')    # a Python with its dependencies (and myokit)

def converted(name: str):
    async def load():
        framework = await get_framework()
        source = (MODELS / f'{name}.ttl').resolve()
        model = framework.make_bondgraph_model(source.as_uri(), source.read_text())
        assert not model.has_issues, [issue.reason for issue in model.issues]
        return model, model.make_cellml_model().to_xml()
    return asyncio.run(load())

def read(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with open(path, newline='') as fp:
        reader = csv.DictReader(fp)
        rows = {row['variable_name']: row for row in reader}
        return list(reader.fieldnames or []), rows

#===============================================================================

def test_the_mapping_has_libcuflynxs_columns_and_the_bond_graph_names(tmp_path):
    model, cellml = converted('chemostat_named')
    path = write_variable_mapping(variable_mapping(model, cellml), tmp_path / 'm_variable_mapping.csv')
    header, rows = read(path)
    assert header == [*COLUMNS, BONDGRAPH_COLUMN]
    by_bondgraph = {row[BONDGRAPH_COLUMN]: row for row in rows.values() if row[BONDGRAPH_COLUMN]}
    # a renamed variable: its CellML name is the key, its bond graph name the extra column, its symbol libcuflynx's
    renamed = by_bondgraph['q_B_B_c']
    assert renamed['variable_name'].endswith('/n_B_c')
    assert renamed['kind'] == 'state'
    assert renamed['latex'] == r'n_{B,c}'
    assert renamed['units']
    # a variable keeping its bond graph name is written in bond graph notation
    kept = by_bondgraph['u_f_Re_1']
    assert kept['variable_name'].endswith('/u_f_Re_1')
    assert kept['latex'] == r'u_{f,\mathrm{Re\_1}}'
    # time first, as libcuflynx orders it
    first = next(iter(rows.values()))
    assert first['kind'] == 'time' and first['latex'] == 't'
    assert {row['kind'] for row in rows.values()} <= {'time', 'state', 'parameter', 'constant', 'algebraic'}

def test_a_species_is_written_in_q_notation(tmp_path):
    model, cellml = converted('chemostat_copies')
    rows = variable_mapping(model, cellml)
    potentials = [row for row in rows if row[BONDGRAPH_COLUMN].startswith('u_Se_A')]
    assert potentials and all(row['latex'] == 'u^{A}_{c}' for row in potentials)
    # an element without a species: its identifier, less a leading repeat of the variable
    q = next(row for row in rows if row[BONDGRAPH_COLUMN] == 'q_B')
    assert q['latex'] == 'q_{B}'
    assert q['kind'] == 'state'

def test_the_mapping_sits_beside_its_cellml_as_libcuflynx_looks_for_it():
    assert mapping_path('out/glycolysis.cellml') == Path('out/glycolysis_variable_mapping.csv')
    assert mapping_path('out/glycolysis_flat.cellml') == Path('out/glycolysis_variable_mapping.csv')

def test_a_module_export_writes_its_mapping(tmp_path):
    model, cellml = converted('chemostat_product')
    files = export_module(model, cellml, tmp_path, 'bgmod', 'bondgraph_v01', source='chemostat_product.ttl',
                          variable_mapping_file=True)
    path = tmp_path / 'bgmod_bondgraph_v01_variable_mapping.csv'
    assert path in files
    header, rows = read(path)
    assert header == [*COLUMNS, BONDGRAPH_COLUMN]
    assert all(name.startswith('bgmod/') for name in rows)
    config = json.loads((tmp_path / 'bgmod_bondgraph_v01_modules_config.json').read_text())[0]
    kinds = {name.split('/', 1)[1]: row['kind'] for name, row in rows.items()}
    for name, _units, _access, role in config['variables_and_units']:
        if role == 'constant':
            assert kinds[name] == 'parameter'
    assert 'state' in kinds.values()

#===============================================================================

LIBCUFLYNX_CHECK = r'''
import json, sys
from libcuflynx.reporting.variable_mapping import default_variable_mapping, read_variable_mapping, variable_mapping
from libcuflynx.reporting.latex_names import symbol
cellml, path, names = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
ours = read_variable_mapping(path)
theirs = {row['variable_name']: row for row in default_variable_mapping(cellml)}
kept = {row['variable_name']: row['latex'] for row in variable_mapping(cellml, path)}
# libcuflynx's schema of the file, when it has one (read_variable_mapping then checks the file against it too)
try:
    from libcuflynx.reporting.variable_mapping import variable_mapping_document
    from libcuflynx.schemas import VARIABLE_MAPPING_SCHEMA
    from libcuflynx.schemas.validation import validate
    validate(VARIABLE_MAPPING_SCHEMA, variable_mapping_document(path), source=path, what='variable mapping')
    schema = True
except ImportError:
    schema = False
bondgraph = {k: r.get('bondgraph') for k, r in ours.items()} if schema else None
print(json.dumps({
    'schema': schema, 'bondgraph': bondgraph,
    'ours': {k: [r['kind'], r['component']] for k, r in ours.items()},
    'theirs': {k: [r['kind'], r['component']] for k, r in theirs.items()},
    'kept': kept, 'latex': {k: r['latex'] for k, r in ours.items()},
    'symbols': {n: symbol(n) for n in names},
}))
'''

@pytest.mark.skipif(not (LIBCUFLYNX_SRC and LIBCUFLYNX_PYTHON),
                    reason='set LIBCUFLYNX_SRC and LIBCUFLYNX_PYTHON to read the mapping with libcuflynx')
def test_libcuflynx_reads_the_mapping_and_keeps_its_symbols(tmp_path):
    model, cellml = converted('chemostat_named')
    cellml_path = tmp_path / 'chemostat_named.cellml'
    cellml_path.write_text(cellml)
    path = write_variable_mapping(variable_mapping(model, cellml), mapping_path(cellml_path))
    names = ['n_B_c', 'mu_A_c_2', 'kappa_Re_1', 'V0', 'Delta_G', 'Alpha_x', 'g_leak_Na', 'K_d']
    result = subprocess.run([LIBCUFLYNX_PYTHON, '-c', LIBCUFLYNX_CHECK, str(cellml_path), str(path), json.dumps(names)],
                            env=dict(os.environ, PYTHONPATH=LIBCUFLYNX_SRC), capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, result.stderr[-3000:]
    out = json.loads(result.stdout.strip().splitlines()[-1])
    # the same variables, kinds and components as libcuflynx's own mapping of the CellML
    assert out['ours'] == out['theirs']
    # a refresh keeps our symbols
    assert out['kept'] == out['latex']
    # the default rule is libcuflynx's
    assert {n: symbol(n) for n in names} == out['symbols']
    # with libcuflynx's schema of the mapping, the file is valid, and its bondgraph column is read
    if out['schema']:
        assert out['bondgraph']['main/n_B_c'] == 'q_B_B_c'

#===============================================================================
#===============================================================================
