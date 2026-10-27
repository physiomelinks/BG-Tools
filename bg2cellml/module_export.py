#===============================================================================
#
#  CellDL and bondgraph tools
#
#  Copyright (c) 2020 - 2025 David Brooks
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#===============================================================================

"""
A bond graph as a module version for libcuflynx (circulatory_autogen) and PhLynx.

``export_module(model, cellml, out_dir, module_type, module_subtype)`` writes the files a module version
has in a module library such as circulatory-autogen-modules (libcuflynx ``utilities/parameterisations.py``),
in the formats of libcuflynx's ``schemas/module_config.schema.json`` (PhLynx's keys)::

    <type>_<subtype>_modules.cellml         the bond graph's equations in one component, component_type
    <type>_<subtype>_units.cellml           the units they use
    <type>_<subtype>_modules_config.json    one entry: module_type, module_subtype, component_file,
                                            component_type, licence, creator, required_citations, ports,
                                            variables_and_units, default_parameterisation
    <type>_<subtype>_references.bib         the required citations' entries (by default a placeholder citing
                                            libcuflynx)
    <type>_<subtype>_tests.yaml             the version's record (reviewed: false)
    <type>_<subtype>_verification_config.json
                                            a short run of the version's states, for the library's tests
    parameterisations/<name>/<name>_parameters.json
                                            every value, libcuflynx's nested JSON (schema_version 1; units,
                                            data_reference and sourced "no")

Every variable is public. A value of the bond graph (a parameter, a potential source's potential) becomes
a parameter of the module; a state ``X`` starts from a parameter ``X_init``.

A signal input port block (bgf:SignalInputPort) is an entrance port of its port type (bgf:portType, by default
``<name>_port``) with one variable, the block's output, a boundary condition whose value is a parameter; a signal
output port block (bgf:SignalOutputPort) is an exit port (multi_port, as a signal may go to several modules) with
its output, a variable.

A chemostat (a potential source with a species) is the module's link to its surroundings, and a bond graph
may have several copies of one, one beside each reaction using it. Each species and location's copies get
one general port, port_type ``chemostat_<species>_<location>``, with variables ``[u_<species>_<location>,
v_<species>_<location>]``: the chemical potential, a boundary condition every copy takes, and the rate at
which the module takes the species from its surroundings (minus the sum of the copies' flows: a source's
flow is the flow into it). The module connected to such a port gives the potential and takes the flow.
"""

#===============================================================================

import json
from pathlib import Path
import re

import lxml.etree as etree

from bg2cellml.bondgraph.framework_support import POTENTIAL_SOURCE
from bg2cellml.variable_mapping import variable_mapping, write_variable_mapping
from bg2cellml.variable_names import naming_rules, rename_cellml, render, rule_for, variable_names

#===============================================================================

CELLML_NS = 'http://www.cellml.org/cellml/1.1#'
MATHML_NS = 'http://www.w3.org/1998/Math/MathML'
CMETA_ID = '{http://www.cellml.org/metadata/1.0#}id'
CMETA2_ID = '{http://www.cellml.org/metadata/2.0#}id'
BGF = 'https://bg-rdf.org/ontologies/bondgraph-framework#'

def C(tag: str) -> str:
    return f'{{{CELLML_NS}}}{tag}'

def M(tag: str) -> str:
    return f'{{{MATHML_NS}}}{tag}'

# The model's elements with a species, from its RDF
SOURCES = """
    SELECT DISTINCT ?uri ?species ?location
    WHERE {
        <%MODEL%> bgf:hasBondElement ?uri .
        ?uri a ?type ; bgf:hasSpecies ?species .
        OPTIONAL { ?uri bgf:hasLocation ?location }
    }"""

# The model's signal port blocks, with their port types
SIGNAL_PORTS = """
    SELECT DISTINCT ?uri ?type ?portType
    WHERE {
        <%MODEL%> bgf:hasSignalBlock ?uri .
        ?uri a ?type .
        FILTER (?type IN (bgf:SignalInputPort, bgf:SignalOutputPort))
        OPTIONAL { ?uri bgf:portType ?portType }
    }"""

NAME_PATTERN = re.compile(r'[^A-Za-z0-9_]')

# Cited when no citation is given: libcuflynx and the module library (circulatory-autogen-modules' entry)
DEFAULT_CITATION = 'argus2026circulatory'
DEFAULT_BIB = '''@misc{argus2026circulatory,
  author = {Argus, Finbar and Maso Talou, Gonzalo D. and others},
  title = {circulatory\\_autogen (libcuflynx)},
  howpublished = {GitHub},
  url = {https://github.com/physiomelinks/circulatory_autogen},
  year = {2026},
  note = {Placeholder software citation for circulatory_autogen (libcuflynx) and this module library, until a release or paper is cited; year = the version cited}
}
'''

#===============================================================================

class ModuleExportError(Exception):
    pass

#===============================================================================

def _states(component: etree._Element) -> set[str]:
    """The variables a <diff/> differentiates."""
    states = set()
    for apply in component.iter(M('apply')):
        children = list(apply)
        if children and children[0].tag == M('diff'):
            for child in children[1:]:
                if child.tag == M('ci') and child.text:
                    states.add(child.text.strip())
    return states

def _defined(component: etree._Element) -> set[str]:
    """The variables an equation defines (``<eq/>`` with a lone variable on its left)."""
    defined = set()
    for apply in component.iter(M('apply')):
        children = list(apply)
        if len(children) == 3 and children[0].tag == M('eq') and children[1].tag == M('ci') and children[1].text:
            defined.add(children[1].text.strip())
    return defined

def _ci(name: str) -> etree._Element:
    ci = etree.Element(M('ci'))
    ci.text = name
    return ci

def _equation(math: etree._Element, lhs: str, rhs: etree._Element):
    apply = etree.SubElement(math, M('apply'))
    etree.SubElement(apply, M('eq'))
    apply.append(_ci(lhs))
    apply.append(rhs)

def _signed_sum(terms: list[tuple[int, str]]) -> etree._Element:
    if len(terms) == 1 and terms[0][0] > 0:
        return _ci(terms[0][1])
    apply = etree.Element(M('apply'))
    etree.SubElement(apply, M('plus'))
    for sign, name in terms:
        if sign > 0:
            apply.append(_ci(name))
        else:
            minus = etree.SubElement(apply, M('apply'))
            etree.SubElement(minus, M('minus'))
            minus.append(_ci(name))
    return apply

def _unique(name: str, taken: set[str]) -> str:
    unique, n = name, 1
    while unique in taken:
        n += 1
        unique = f'{name}_{n}'
    taken.add(unique)
    return unique

#===============================================================================

def chemostat_groups(model) -> dict[tuple[str, str], list[tuple[str, str, int]]]:
    """
    ``{(species, location): [(potential variable, flow variable, sign), ...]}`` for the model's potential
    sources with a species: a source's flow is the flow into it, whichever way its bond points, so the
    module's uptake from its surroundings is minus its flow (sign -1).
    """

    names = variable_names(model)               # as the CellML names them
    places = {}
    for row in model.rdf_graph.query(SOURCES.replace('%MODEL%', model.uri.value)):
        species = row['species'].value
        location = row['location'].value if row.get('location') is not None else ''
        places[row['uri'].value] = (species, location)
    groups: dict[tuple[str, str], list[tuple[str, str, int]]] = {}
    for element in model.elements:
        if element.element_class != POTENTIAL_SOURCE or element.uri.value not in places:
            continue
        if element.potential is None or element.flow is None:
            continue
        sign = -1
        u, v = element.potential.symbol, element.flow.symbol
        groups.setdefault(places[element.uri.value], []).append((names.get(u, u), names.get(v, v), sign))
    return groups

def chemostat_domains(model) -> dict[tuple[str, str], str]:
    """{(species, location): the physical domain's URI} of the model's potential sources with a species."""
    places = {}
    for row in model.rdf_graph.query(SOURCES.replace('%MODEL%', model.uri.value)):
        location = row['location'].value if row.get('location') is not None else ''
        places[row['uri'].value] = (row['species'].value, location)
    return {places[e.uri.value]: e.domain.id for e in model.elements
            if e.element_class == POTENTIAL_SOURCE and e.uri.value in places}

#===============================================================================

def signal_ports(model) -> list[tuple[str, str, str]]:
    """``[(direction, CellML variable, port type), ...]`` of the model's signal port blocks: direction
    ``entrance`` (a bgf:SignalInputPort) or ``exit`` (a bgf:SignalOutputPort)."""
    names = variable_names(model)
    blocks = {block.uri.value: block for block in model.signal_blocks}
    ports = []
    for row in model.rdf_graph.query(SIGNAL_PORTS.replace('%MODEL%', model.uri.value)):
        if (block := blocks.get(row['uri'].value)) is None:
            continue
        name = names.get(block.output.symbol, block.output.symbol)
        port_type = row['portType'].value if row.get('portType') is not None else f'{name}_port'
        direction = 'entrance' if row['type'].value == f'{BGF}SignalInputPort' else 'exit'
        ports.append((direction, name, port_type))
    return sorted(ports)

def _number(value: str) -> float|int:
    number = float(value)
    return int(number) if number.is_integer() and abs(number) < 1e15 else number

def _bib_keys(bib: str) -> list[str]:
    return re.findall(r'@\w+\s*\{\s*([^,\s]+)\s*,', bib)

def export_module(model, cellml: str, out_dir: str|Path, module_type: str, module_subtype: str,
                  component_type: str|None=None, parameterisation: str='default',
                  source: str='bond graph', citations: list[str]|None=None, bib: str|None=None,
                  licence: str='CC0-1.0', creator: list[str]|None=None,
                  variable_mapping_file: bool=False) -> list[Path]:
    """
    Write the module version of ``model`` (a BondgraphModel) and its CellML (from ``make_cellml_model()``)
    into ``out_dir``; returns the files written. ``source`` names the bond graph in each parameter's
    data_reference. ``bib`` is the text of the version's references.bib and ``citations`` the keys of it the
    module requires (all its keys when not given); without them it cites libcuflynx's placeholder entry.
    ``variable_mapping_file`` also writes ``<type>_<subtype>_variable_mapping.csv`` (a module library such as
    circulatory-autogen-modules keeps no such file in a version).
    """
    for name, what in ((module_type, 'module_type'), (module_subtype, 'module_subtype')):
        if not name or NAME_PATTERN.search(name):
            raise ModuleExportError(f'{what} {name!r} must be letters, digits and _ only')
    component_type = component_type or module_type
    prefix = f'{module_type}_{module_subtype}'
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    document = etree.fromstring(cellml.encode())
    components = document.findall(C('component'))
    if len(components) != 1:
        raise ModuleExportError(f'expected the one component bg2cellml writes, found {len(components)}')
    component = components[0]
    component.set('name', component_type)
    math = component.find(M('math'))
    if math is None:
        math = etree.SubElement(component, M('math'))

    variables = {v.get('name'): v for v in component.findall(C('variable'))}
    taken = set(variables)
    states = _states(component)
    rows: list[list[str]] = []
    variables_and_units: list[list[str]] = []
    ports = []
    interfaces: dict[str, str] = {}         # each port variable's public interface

    def parameter(name: str, units: str, value: str, reference: str):
        rows.append([name, units, value, reference])

    # One port per chemostat species and location: a shared potential and the net uptake, named by the model's
    # rule for a potential source's u and v when it has one
    rules = naming_rules(model)
    domains = chemostat_domains(model)
    def port_name(var: str, place: str, species: str, location: str) -> str:
        domain = domains.get((species, location), '')
        if (pattern := rule_for(rules, POTENTIAL_SOURCE, domain, var)) is not None:
            return render(pattern, var=var, element=place, label='', species=species, location=location, name=species)
        return f'{var}_{place}'
    for (species, location), copies in sorted(chemostat_groups(model).items()):
        place = f'{species}_{location}' if location else species
        place = NAME_PATTERN.sub('_', place)
        # a copy already named as the port is renamed, so the port keeps its name
        wanted = {port_name('u', place, species, location), port_name('v', place, species, location)}
        clashes = {name: _unique(f'{name}_copy', taken) for copy in copies for name in copy[:2] if name in wanted}
        if clashes:
            rename_cellml(component, clashes)
            for old, new in clashes.items():
                variables[new] = variables.pop(old)
                taken.discard(old)
            copies = [(clashes.get(u, u), clashes.get(v, v), sign) for u, v, sign in copies]
        values = {variables[u].get('initial_value') for u, _, _ in copies if u in variables}
        values.discard(None)
        if len(values) > 1:
            raise ModuleExportError(f'the copies of chemostat {species}[{location}] have different potentials: '
                                    f'{", ".join(sorted(values))}')
        u_units = variables[copies[0][0]].get('units')
        v_units = variables[copies[0][1]].get('units')
        u_port = _unique(port_name('u', place, species, location), taken)
        v_port = _unique(port_name('v', place, species, location), taken)
        for u, _, _ in copies:
            variables[u].attrib.pop('initial_value', None)
            _equation(math, u, _ci(u_port))
        _equation(math, v_port, _signed_sum([(sign, v) for _, v, sign in copies]))
        for name, units in ((u_port, u_units), (v_port, v_units)):
            variables[name] = etree.SubElement(component, C('variable'), name=name, units=units)
        if values:
            parameter(u_port, u_units, values.pop(), f'{source}: chemostat {species}[{location}] potential')
        variables_and_units.append([u_port, u_units, 'access', 'boundary_condition'])
        variables_and_units.append([v_port, v_units, 'access', 'variable'])
        ports.append({'port_type': f'chemostat_{place}', 'variables': [u_port, v_port]})
        interfaces[u_port] = 'in'
        interfaces[v_port] = 'out'

    # Signal ports: an entrance port's variable is a boundary condition (its value a parameter), an exit port's
    # a variable the module computes
    entrance_ports = []
    exit_ports = []
    for direction, name, port_type in signal_ports(model):
        if name not in variables:
            raise ModuleExportError(f'signal port variable {name} is not in the CellML')
        units = variables[name].get('units')
        if direction == 'entrance':
            if (value := variables[name].attrib.pop('initial_value', None)) is not None:
                parameter(name, units, value, f'{source}: {name}, until the {port_type} entrance port sets it')
            variables_and_units.append([name, units, 'access', 'boundary_condition'])
            entrance_ports.append({'port_type': port_type, 'variables': [name]})
            interfaces[name] = 'in'
        else:
            variables_and_units.append([name, units, 'access', 'variable'])
            exit_ports.append({'port_type': port_type, 'variables': [name], 'multi_port': 'True'})
            interfaces[name] = 'out'

    # Every other value becomes a parameter, a state's initial value a parameter <state>_init
    port_variables = set(interfaces)
    ordered: list[etree._Element] = []
    for name, variable in list(variables.items()):
        ordered.append(variable)
        if variable.get(CMETA_ID) is not None:
            del variable.attrib[CMETA_ID]
        if variable.get(CMETA2_ID) is not None:
            del variable.attrib[CMETA2_ID]
        if name in port_variables:
            continue
        units = variable.get('units')
        value = variable.get('initial_value')
        if name == 't':
            variable.set('public_interface', 'in')
            continue
        if name in states:
            init = _unique(f'{name}_init', taken)
            variable.set('initial_value', init)
            ordered.append(etree.Element(C('variable'), name=init, units=units, public_interface='in'))
            if value is not None:
                parameter(init, units, value, f'{source}: initial {name}')
            variables_and_units.append([init, units, 'access', 'constant'])
            variable.set('public_interface', 'out')
            variables_and_units.append([name, units, 'access', 'variable'])
        elif value is not None:
            del variable.attrib['initial_value']
            variable.set('public_interface', 'in')
            parameter(name, units, value, f'{source}: {name}')
            variables_and_units.append([name, units, 'access', 'constant'])
        else:
            variable.set('public_interface', 'out')
            variables_and_units.append([name, units, 'access', 'variable'])
    for name, interface in interfaces.items():
        variables[name].set('public_interface', interface)
    # the variables come first in a component
    for variable in component.findall(C('variable')):
        component.remove(variable)
    for variable in reversed(ordered):
        component.insert(0, variable)
    undefined = [v.get('name') for v in component.findall(C('variable'))
                 if v.get('public_interface') == 'out' and v.get('name') not in states | _defined(component)
                 and v.get('name') not in port_variables]
    if undefined:
        raise ModuleExportError(f'variables neither set nor computed: {", ".join(undefined)}')

    # Units in their own file
    units_model = etree.Element(C('model'), name='Units', nsmap={None: CELLML_NS, 'cellml': CELLML_NS})
    for units in document.findall(C('units')):
        document.remove(units)
        units_model.append(units)
    document.set('name', prefix)
    for element in document.iter():
        for key in (CMETA_ID, CMETA2_ID):
            if key in element.attrib:
                del element.attrib[key]
    etree.cleanup_namespaces(document)

    files = []
    def write_xml(path: Path, root: etree._Element):
        path.write_bytes(etree.tostring(root, xml_declaration=True, encoding='UTF-8', pretty_print=True))
        files.append(path)

    write_xml(out / f'{prefix}_modules.cellml', document)
    if variable_mapping_file:
        files.append(write_variable_mapping(variable_mapping(model, etree.tostring(document).decode(),
                                                             parameters={row[0] for row in rows} - port_variables),
                                            out / f'{prefix}_variable_mapping.csv'))
    write_xml(out / f'{prefix}_units.cellml', units_model)

    if bib is None:
        bib = DEFAULT_BIB
    keys = _bib_keys(bib)
    citations = list(citations) if citations else keys
    if not citations or any(key not in keys for key in citations):
        raise ModuleExportError(f'required citations {citations} must be keys of the references ({keys})')

    config = [{
        'module_type': module_type,
        'module_subtype': module_subtype,
        'module_format': 'cellml',
        'component_file': f'{prefix}_modules.cellml',
        'component_type': component_type,
        'default_parameterisation': parameterisation,
        'licence': licence,
        'creator': list(creator or []),
        'required_citations': citations,
        'entrance_ports': entrance_ports,
        'exit_ports': exit_ports,
        'general_ports': ports,
        'variables_and_units': variables_and_units,
        'notes': f'Exported by bg2cellml from {source}.'
    }]
    config_path = out / f'{prefix}_modules_config.json'
    config_path.write_text(json.dumps(config, indent=2) + '\n')
    files.append(config_path)

    bib_path = out / f'{prefix}_references.bib'
    bib_path.write_text(bib if bib.endswith('\n') else bib + '\n')
    files.append(bib_path)

    tests_path = out / f'{prefix}_tests.yaml'
    notes = f'Exported by bg2cellml from {source}; parameter values unreviewed.'.replace("'", "''")
    tests_path.write_text(f"module_type: {module_type}\nversion: {module_subtype}\nreviewed: false\n"
                          f"notes: '{notes}'\n")
    files.append(tests_path)

    verification = {
        'module_type': module_type,
        'version': module_subtype,
        'sim_time': 10.0,
        'pre_time': 0.0,
        'dt': 0.01,
        'solver': 'CVODE_myokit',
        'solver_info': {'rtol': 1e-08, 'atol': 1e-10},
        'outputs': sorted(states)[:10] or [name for name, *_ in variables_and_units if name not in port_variables][:1],
        'validation': {parameterisation: {'baseline': {
            'status': 'pending', 'reason': f'exported by bg2cellml from {source}; no module-level data recorded'}}}
    }
    verification_path = out / f'{prefix}_verification_config.json'
    verification_path.write_text(json.dumps(verification, indent=2) + '\n')
    files.append(verification_path)

    parameters_dir = out / 'parameterisations' / parameterisation
    parameters_dir.mkdir(parents=True, exist_ok=True)
    parameters_path = parameters_dir / f'{parameterisation}_parameters.json'
    parameters: dict = {'schema_version': 1}
    for name, units, value, reference in rows:
        parameters[name] = {'value': _number(value), 'units': units, 'data_reference': reference, 'sourced': 'no'}
    parameters_path.write_text(json.dumps(parameters, indent=1) + '\n')
    files.append(parameters_path)
    return files

#===============================================================================
#===============================================================================
