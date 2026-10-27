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

A flow or potential source with a bgf:portType (and no species) is a port of that type with variables
``[flow, potential]`` (a vessel_port's order): an entrance port when its bond leaves it, else an exit port, the
flow positive into the module at an entrance and out of it at an exit. The source's value is a boundary condition
(a parameter); the other variable the module computes.

A signal input port block (bgf:SignalInputPort) is an entrance port of its port type (bgf:portType, by default
``<name>_port``) with one variable, the block's output, a boundary condition whose value is a parameter; a signal
output port block (bgf:SignalOutputPort) is an exit port (multi_port, as a signal may go to several modules) with
its output, a variable.

Any other bond element or signal block with a bgf:portType is an exit port of its quantity (a store), its output (a
block) or its flow. An element's bgf:portVariables names the variables it gives its port instead (by their template
names: ``"v"``, ``"q"``, ``"u v"``, a parameter), bgf:portDirection makes the port an ``entrance``, ``exit`` or
``general`` port and bgf:multiPort sets its multi_port. The elements with one port type and direction that have a
bgf:portOrder make one port, their variables in that order. A port variable with a value (a source's value, a
parameter, a signal input port's value) is a boundary condition and its value a parameter; a state stays a state
(from its ``_init`` parameter); any other variable the module computes.

``reference``, an existing version of the module (its directory), keeps that version's form where the bond graph
does not say: the names of its units (each variable's, when the definitions agree), its parameters' data_reference,
sourced and other keys (the values are the bond graph's; a different value is reported), its config's other keys
(e.g. required_citations_uncertain), licence, creators, citations and references, its tests.yaml and
verification_config.json, and its order of variables, parameters and ports.

A chemostat (a potential source with a species) is the module's link to its surroundings, and a bond graph
may have several copies of one, one beside each reaction using it. Each species and location's copies get
one general port, port_type ``chemostat_<species>_<location>``, with variables ``[u_<species>_<location>,
v_<species>_<location>]``: the chemical potential, a boundary condition every copy takes, and the rate at
which the module takes the species from its surroundings (minus the sum of the copies' flows: a source's
flow is the flow into it). The module connected to such a port gives the potential and takes the flow.
"""

#===============================================================================

import copy
import json
from pathlib import Path
import re

import lxml.etree as etree

from bg2cellml.bondgraph.framework_support import FLOW_SOURCE, POTENTIAL_SOURCE, QUANTITY_STORE
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

# The model's sources with a port type, and whether a power bond leaves them
SOURCE_PORTS = """
    SELECT DISTINCT ?uri ?portType ?outward
    WHERE {
        <%MODEL%> bgf:hasBondElement ?uri .
        ?uri bgf:portType ?portType .
        FILTER NOT EXISTS { ?uri bgf:hasSpecies ?species }
        BIND (EXISTS { ?bond bgf:hasSource ?uri . FILTER NOT EXISTS { ?bond a bgf:SignalBond } } AS ?outward)
    }"""

# How an element or block with a port type is a port: its variables, direction, place in a port and multi_port
PORT_OPTIONS = """
    SELECT DISTINCT ?uri ?portType ?variables ?direction ?order ?multi
    WHERE {
        { <%MODEL%> bgf:hasBondElement ?uri } UNION { <%MODEL%> bgf:hasSignalBlock ?uri }
        ?uri bgf:portType ?portType .
        OPTIONAL { ?uri bgf:portVariables ?variables }
        OPTIONAL { ?uri bgf:portDirection ?direction }
        OPTIONAL { ?uri bgf:portOrder ?order }
        OPTIONAL { ?uri bgf:multiPort ?multi }
    }"""

PORT_DIRECTIONS = ('entrance', 'exit', 'general')

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
    return [(direction, name, port_type) for direction, name, port_type, _ in _signal_port_blocks(model)]

def _signal_port_blocks(model) -> list[tuple[str, str, str, str]]:
    names = variable_names(model)
    blocks = {block.uri.value: block for block in model.signal_blocks}
    ports = []
    for row in model.rdf_graph.query(SIGNAL_PORTS.replace('%MODEL%', model.uri.value)):
        if (block := blocks.get(row['uri'].value)) is None:
            continue
        name = names.get(block.output.symbol, block.output.symbol)
        port_type = row['portType'].value if row.get('portType') is not None else f'{name}_port'
        direction = 'entrance' if row['type'].value == f'{BGF}SignalInputPort' else 'exit'
        ports.append((direction, name, port_type, block.uri.value))
    return sorted(ports)

def port_options(model) -> dict[str, dict]:
    """``{element or block URI: {variables, direction, order, multi}}`` of the model's elements with a port type, as
    their bgf:portVariables (a list of template names), bgf:portDirection, bgf:portOrder and bgf:multiPort give
    them (None where not given)."""
    options = {}
    for row in model.rdf_graph.query(PORT_OPTIONS.replace('%MODEL%', model.uri.value)):
        variables = row['variables'].value.split() if row.get('variables') is not None else None
        direction = row['direction'].value.strip().lower() if row.get('direction') is not None else None
        if direction is not None and direction not in PORT_DIRECTIONS:
            raise ModuleExportError(f'{row["uri"].value}: bgf:portDirection {direction!r} is not one of '
                                    f'{", ".join(PORT_DIRECTIONS)}')
        order = row['order'].value if row.get('order') is not None else None
        options[row['uri'].value] = {
            'port_type': row['portType'].value,
            'variables': variables,
            'direction': direction,
            'order': float(order) if order is not None else None,
            'multi': row['multi'].value if row.get('multi') is not None else None,
        }
    return options

def _number(value: str) -> float|int:
    number = float(value)
    return int(number) if number.is_integer() and abs(number) < 1e15 else number

def _bib_keys(bib: str) -> list[str]:
    return re.findall(r'@\w+\s*\{\s*([^,\s]+)\s*,', bib)

def export_module(model, cellml: str, out_dir: str|Path, module_type: str, module_subtype: str,
                  component_type: str|None=None, parameterisation: str='default',
                  source: str='bond graph', citations: list[str]|None=None, bib: str|None=None,
                  licence: str|None=None, creator: list[str]|None=None,
                  variable_mapping_file: bool=False, reference: str|Path|None=None,
                  warn=None) -> list[Path]:
    """
    Write the module version of ``model`` (a BondgraphModel) and its CellML (from ``make_cellml_model()``)
    into ``out_dir``; returns the files written. ``source`` names the bond graph in each parameter's
    data_reference. ``bib`` is the text of the version's references.bib and ``citations`` the keys of it the
    module requires (all its keys when not given); without them it cites libcuflynx's placeholder entry.
    ``variable_mapping_file`` also writes ``<type>_<subtype>_variable_mapping.csv`` (a module library such as
    circulatory-autogen-modules keeps no such file in a version). ``reference``: an existing version's directory
    whose form the export keeps (see the module's description); ``warn(message)`` is told of parameter values that
    differ from it.
    """
    for name, what in ((module_type, 'module_type'), (module_subtype, 'module_subtype')):
        if not name or NAME_PATTERN.search(name):
            raise ModuleExportError(f'{what} {name!r} must be letters, digits and _ only')
    prefix = f'{module_type}_{module_subtype}'
    ref = ReferenceVersion(reference, prefix) if reference is not None else None
    if ref is not None:
        component_type = component_type or ref.config.get('component_type')
        licence = licence or ref.config.get('licence')
        creator = creator or ref.config.get('creator')
        if bib is None and ref.bib is not None:
            bib = ref.bib
            citations = citations or ref.config.get('required_citations')
    component_type = component_type or module_type
    licence = licence or 'CC0-1.0'
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
    # a domain constant (RT, F) no equation uses is left out
    used = {ci.text.strip() for ci in component.iter(M('ci')) if ci.text}
    constants = {constant.symbol for item in [*model.elements, *model.signal_blocks] if item.domain is not None
                 for constant in item.domain.constants}
    for name in sorted(constants & set(variables) - used):
        component.remove(variables.pop(name))
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

    entrance_ports = []
    exit_ports = []
    general_ports = []
    port_lists = {'entrance': entrance_ports, 'exit': exit_ports, 'general': general_ports}
    port_states: set[str] = set()           # states in ports, which stay states
    options = port_options(model)
    groups: dict[tuple[str, str], dict] = {}

    def add_port(uri: str|None, port_type: str, direction: str, port_variables_: list[str], multi: str|None=None):
        option = options.get(uri, {}) if uri is not None else {}
        multi = option.get('multi') or multi
        if option.get('order') is not None:
            # one port of the elements with this type and direction that have a bgf:portOrder
            key = (port_type, direction)
            if key not in groups:
                groups[key] = {'port_type': port_type, 'variables': [], 'members': []}
                port_lists[direction].append(groups[key])
            group = groups[key]
            group['members'].append((option['order'], len(group['members']), port_variables_))
            if multi is not None:
                group['multi_port'] = multi
            return
        port = {'port_type': port_type, 'variables': port_variables_}
        if multi is not None:
            port['multi_port'] = multi
        port_lists[direction].append(port)

    def port_variable(name: str, port_type: str, direction: str, value_is_input: bool=True):
        """A port variable's kind: a state stays a state, one with a value is a boundary condition."""
        if name not in variables:
            raise ModuleExportError(f'port variable {name} of {port_type} is not in the CellML')
        units = variables[name].get('units')
        if name in states:
            port_states.add(name)
        elif value_is_input and (value := variables[name].attrib.pop('initial_value', None)) is not None:
            parameter(name, units, value, f'{source}: {name}, until the {port_type} {direction} port sets it')
            variables_and_units.append([name, units, 'access', 'boundary_condition'])
            interfaces[name] = 'in'
        elif name not in interfaces:
            variables_and_units.append([name, units, 'access', 'variable'])
            interfaces[name] = 'out'

    def element_variable(element, name: str) -> str|None:
        """The CellML name of an element's variable, by its template name (q, u, v or a parameter's)."""
        if name == element.domain.flow.symbol and element.flow is not None:
            symbol = element.flow.symbol
        elif name == element.domain.potential.symbol and element.potential is not None:
            symbol = element.potential.symbol
        elif name == element.domain.quantity.symbol and element.element_class == QUANTITY_STORE:
            symbol = element.variables[name].symbol if name in element.variables else None
            if symbol is None:
                symbol = element.signal_variable(None).symbol
        elif (variable := element.variables.get(name)) is not None:
            symbol = variable.symbol
        else:
            return None
        return names.get(symbol, symbol)

    def node_potential(element, u: str) -> str:
        """The potential of the node a flow source feeds, for its port: bg2cellml gives the source the potential
        conjugate to its flow out (minus the node's), so the port's variable, named as the source's, is minus it and
        the source's own is renamed."""
        own = _unique(f'{u}_source', taken)
        rename_cellml(component, {u: own})
        variables[own] = variables.pop(u)
        variables[u] = etree.SubElement(component, C('variable'), name=u, units=variables[own].get('units'))
        _equation(math, u, _signed_sum([(-1, own)]))
        return u

    # Source ports: [flow, potential], the source's value a boundary condition
    names = variable_names(model)
    sources = {row['uri'].value: (row['portType'].value, row['outward'].value in ('true', '1'))
               for row in model.rdf_graph.query(SOURCE_PORTS.replace('%MODEL%', model.uri.value))}
    for element in model.elements:
        if element.uri.value not in sources or element.element_class not in (FLOW_SOURCE, POTENTIAL_SOURCE):
            continue
        port_type, entrance = sources[element.uri.value]
        option = options.get(element.uri.value, {})
        direction = option.get('direction') or ('entrance' if entrance else 'exit')
        if option.get('variables') is not None:
            port_variables_ = []
            for var in option['variables']:
                if (name := element_variable(element, var)) is None:
                    raise ModuleExportError(f'{element.pretty_name} has no variable {var} for its {port_type} port')
                if element.element_class == FLOW_SOURCE and entrance and var == element.domain.potential.symbol:
                    name = node_potential(element, name)
                port_variables_.append(name)
                port_variable(name, port_type, direction)
            add_port(element.uri.value, port_type, direction, port_variables_)
            continue
        u = names.get(element.potential.symbol, element.potential.symbol)
        v = names.get(element.flow.symbol, element.flow.symbol)
        if element.element_class == POTENTIAL_SOURCE and entrance:
            # the source's flow is the flow into it: the module takes in minus it
            flow = _unique(f'v_into_{NAME_PATTERN.sub("_", element.symbol)}', taken)
            variables[flow] = etree.SubElement(component, C('variable'), name=flow, units=variables[v].get('units'))
            _equation(math, flow, _signed_sum([(-1, v)]))
        else:
            flow = v
        if element.element_class == FLOW_SOURCE and entrance:
            u = node_potential(element, u)
        port_variables_ = [flow, u]
        for name in port_variables_:
            port_variable(name, port_type, direction)
        add_port(element.uri.value, port_type, direction, port_variables_)

    # Signal ports: an entrance port's variable is a boundary condition (its value a parameter), an exit port's
    # a variable the module computes
    for kind, name, port_type, uri in _signal_port_blocks(model):
        if name not in variables:
            raise ModuleExportError(f'signal port variable {name} is not in the CellML')
        direction = options.get(uri, {}).get('direction') or kind
        units = variables[name].get('units')
        if kind == 'entrance':
            if (value := variables[name].attrib.pop('initial_value', None)) is not None:
                parameter(name, units, value, f'{source}: {name}, until the {port_type} {direction} port sets it')
            variables_and_units.append([name, units, 'access', 'boundary_condition'])
            interfaces[name] = 'in'
            add_port(uri, port_type, direction, [name])
        else:
            variables_and_units.append([name, units, 'access', 'variable'])
            interfaces[name] = 'out'
            add_port(uri, port_type, direction, [name], multi='True')

    # Other elements and blocks with a port type: a store's quantity, a block's output, another element's flow
    signal_port_uris = {uri for *_, uri in _signal_port_blocks(model)}
    for element in model.elements:
        uri = element.uri.value
        if uri not in options or element.element_class in (FLOW_SOURCE, POTENTIAL_SOURCE) and uri in sources:
            continue
        option = options[uri]
        port_type, direction = option['port_type'], option.get('direction') or 'exit'
        wanted = option.get('variables') or [element.domain.quantity.symbol if element.element_class == QUANTITY_STORE
                                             else element.domain.flow.symbol]
        port_variables_ = []
        for var in wanted:
            if (name := element_variable(element, var)) is None:
                raise ModuleExportError(f'{element.pretty_name} has no variable {var} for its {port_type} port')
            port_variables_.append(name)
            port_variable(name, port_type, direction)
        add_port(uri, port_type, direction, port_variables_)
    for block in model.signal_blocks:
        uri = block.uri.value
        if uri not in options or uri in signal_port_uris:
            continue
        option = options[uri]
        port_type, direction = option['port_type'], option.get('direction') or 'exit'
        port_variables_ = []
        for var in option.get('variables') or ['y']:
            variable = block.output if var == 'y' else block.variables.get(var)
            if variable is None:
                raise ModuleExportError(f'{block.pretty_name} has no variable {var} for its {port_type} port')
            name = names.get(variable.symbol, variable.symbol)
            port_variables_.append(name)
            port_variable(name, port_type, direction)
        add_port(uri, port_type, direction, port_variables_)

    # ports of several elements, their variables in bgf:portOrder
    for group in groups.values():
        for _, _, port_variables_ in sorted(group.pop('members'), key=lambda member: member[:2]):
            group['variables'].extend(port_variables_)

    # Every other value becomes a parameter, a state's initial value a parameter <state>_init
    port_variables = set(interfaces) - port_states
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
    if ref is not None:
        ref.rename_units(component, variables_and_units, rows)
        ref.order(variables_and_units, rows, port_lists)
        if warn is not None:
            for message in ref.value_differences(rows):
                warn(message)
    if undefined:
        raise ModuleExportError(f'variables neither set nor computed: {", ".join(undefined)}')

    # Units in their own file
    units_model = etree.Element(C('model'), name='Units', nsmap={None: CELLML_NS, 'cellml': CELLML_NS})
    for units in document.findall(C('units')):
        document.remove(units)
        units_model.append(units)
    if ref is not None:
        ref.units_definitions(units_model, component)
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
        'general_ports': ports + general_ports,
        'variables_and_units': variables_and_units,
        'notes': f'Exported by bg2cellml from {source}.'
    }]
    if ref is not None:
        # the reference's other keys (e.g. required_citations_uncertain), in its order
        merged = {key: config[0].get(key, value) for key, value in ref.config.items()}
        merged.update({key: value for key, value in config[0].items() if key not in merged})
        if 'notes' not in ref.config:
            merged.pop('notes', None)
        config = [merged]
    config_path = out / f'{prefix}_modules_config.json'
    config_path.write_text(json.dumps(config, indent=2) + '\n')
    files.append(config_path)

    bib_path = out / f'{prefix}_references.bib'
    bib_path.write_text(bib if bib.endswith('\n') else bib + '\n')
    files.append(bib_path)

    tests_path = out / f'{prefix}_tests.yaml'
    notes = f'Exported by bg2cellml from {source}; parameter values unreviewed.'.replace("'", "''")
    if ref is not None and ref.tests is not None:
        tests_path.write_text(ref.tests)
    else:
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
    if ref is not None and ref.verification is not None:
        verification_path.write_text(ref.verification)
    else:
        verification_path.write_text(json.dumps(verification, indent=2) + '\n')
    files.append(verification_path)

    parameters_dir = out / 'parameterisations' / parameterisation
    parameters_dir.mkdir(parents=True, exist_ok=True)
    parameters_path = parameters_dir / f'{parameterisation}_parameters.json'
    parameters: dict = {'schema_version': 1}
    for name, units, value, reference in rows:
        parameters[name] = {'value': _number(value), 'units': units, 'data_reference': reference, 'sourced': 'no'}
        if ref is not None and name in ref.parameters:
            # the reference's metadata (data_reference, sourced, ...), with the bond graph's value and units
            parameters[name] = {**ref.parameters[name], 'value': _number(value), 'units': units}
    parameters_path.write_text(json.dumps(parameters, indent=1) + '\n')
    files.append(parameters_path)
    return files

#===============================================================================

class ReferenceVersion:
    """An existing version of a module, whose form an export keeps: names of units, parameter metadata, config
    keys, references, tests.yaml, verification_config.json and order."""
    def __init__(self, directory: str|Path, prefix: str):
        d = Path(directory)
        def text(suffix: str) -> str|None:
            path = d / f'{prefix}{suffix}'
            return path.read_text() if path.exists() else None
        config = text('_modules_config.json')
        self.config: dict = json.loads(config)[0] if config else {}
        self.bib = text('_references.bib')
        self.tests = text('_tests.yaml')
        self.verification = text('_verification_config.json')
        name = self.config.get('default_parameterisation', 'default')
        path = d / 'parameterisations' / name / f'{name}_parameters.json'
        self.parameters: dict = {}
        if path.exists():
            self.parameters = {k: v for k, v in json.loads(path.read_text()).items()
                               if k != 'schema_version' and isinstance(v, dict) and 'value' in v}
        self.units_of = {row[0]: row[1] for row in self.config.get('variables_and_units', [])}
        self.__units: dict[str, etree._Element] = {}
        units_text = text('_units.cellml')
        if units_text:
            for units in etree.fromstring(units_text.encode()).iter(C('units')):
                self.__units[units.get('name')] = units

    def rename_units(self, component: etree._Element, variables_and_units: list, rows: list):
        """Name each variable's units as the reference does (its own name for the variable, else the one
        reference units of the same definition), where the definitions agree."""
        definitions = {units.get('name'): units for units in component.getparent().iter(C('units'))}
        renamed: dict[str, str] = {}
        for variable in component.findall(C('variable')):
            name, units = variable.get('name'), variable.get('units')
            wanted = self.units_of.get(name)
            if wanted is None or wanted == units:
                matches = [n for n in self.__units if _same_units(self, n, definitions, units)]
                wanted = matches[0] if len(matches) == 1 and units not in self.__units else None
            if wanted is not None and wanted != units and _same_units(self, wanted, definitions, units):
                variable.set('units', wanted)
                renamed[name] = wanted
        for row in variables_and_units:
            row[1] = renamed.get(row[0], row[1])
        for row in rows:
            row[1] = renamed.get(row[0], row[1])

    def units_definitions(self, units_model: etree._Element, component: etree._Element):
        """The units file: the reference's definitions of the names the component uses, then the export's own."""
        used = {variable.get('units') for variable in component.findall(C('variable'))}
        used |= {element.get(f'{{{CELLML_NS}}}units') for element in component.iter()
                 if element.get(f'{{{CELLML_NS}}}units') is not None}
        exported = {units.get('name'): units for units in units_model.findall(C('units'))}
        def needs(name: str, seen: set) -> set:
            # a units' definition's own units
            if name in seen:
                return seen
            seen.add(name)
            definition = self.__units.get(name, exported.get(name))
            if definition is not None:
                for unit in definition.findall(C('unit')):
                    needs(unit.get('units'), seen)
            return seen
        wanted: set = set()
        for name in used:
            needs(name, wanted)
        for units in list(units_model):
            units_model.remove(units)
        for name, units in self.__units.items():
            if name in wanted:
                units_model.append(copy.deepcopy(units))
        for name, units in exported.items():
            if name in wanted and name not in self.__units:
                units_model.append(units)

    def order(self, variables_and_units: list, rows: list, port_lists: dict):
        """The reference's order of variables, parameters and ports, then the export's own."""
        position: dict[str, int] = {}
        for n, row in enumerate(self.config.get('variables_and_units', [])):
            position.setdefault(row[0], n)
        variables_and_units.sort(key=lambda row: position.get(row[0], len(position)))
        names = list(self.parameters)
        rows.sort(key=lambda row: names.index(row[0]) if row[0] in names else len(names))
        for direction, ports in port_lists.items():
            known = [p['port_type'] for p in self.config.get(f'{direction}_ports', [])]
            ports.sort(key=lambda p: known.index(p['port_type']) if p['port_type'] in known else len(known))

    def value_differences(self, rows: list) -> list[str]:
        messages = []
        for name, _, value, _ in rows:
            if name in self.parameters:
                old = self.parameters[name].get('value')
                try:
                    differs = abs(float(value) - float(old)) > 1e-12*max(abs(float(old)), abs(float(value)), 1e-300)
                except (TypeError, ValueError):
                    differs = str(value) != str(old)
                if differs:
                    messages.append(f'{name}: the bond graph has {value}, the reference {old}')
        return messages

    def definitions(self) -> dict[str, etree._Element]:
        return dict(self.__units)

#===============================================================================

_BUILTIN = {'litre': 'liter', 'liter': 'liter', 'metre': 'meter', 'meter': 'meter', 'dimensionless': 'dimensionless'}

def _units_quantity(name: str, definitions: dict, registry, depth: int=0):
    """A pint quantity for CellML units ``name`` (defined in ``definitions``, or built in)."""
    if depth > 20:
        raise ValueError(f'units {name} are defined in terms of themselves')
    definition = definitions.get(name)
    if definition is None:
        if name in ('dimensionless', 'dim'):
            return registry.Quantity(1.0)
        return registry.Quantity(1.0, _BUILTIN.get(name, name))
    quantity = registry.Quantity(1.0)
    for unit in definition.findall(C('unit')):
        base = _units_quantity(unit.get('units'), definitions, registry, depth + 1)
        prefix = unit.get('prefix')
        if prefix:
            factor = registry.Quantity(1.0, f'{prefix}meter').to('meter').magnitude if not prefix.lstrip('-').isdigit() \
                else 10.0**int(prefix)
            base = base*factor
        base = base*float(unit.get('multiplier', 1.0))
        quantity = quantity*base**float(unit.get('exponent', 1.0))
    return quantity

def _same_units(old: ReferenceVersion, old_name: str, definitions: dict, new_name: str|None) -> bool:
    """Whether the reference's units ``old_name`` and the export's ``new_name`` (in ``definitions``) agree."""
    if new_name is None:
        return False
    import pint
    registry = _registry()
    try:
        a = _units_quantity(old_name, {**definitions, **old.definitions()}, registry)
        b = _units_quantity(new_name, definitions, registry)
        if a.dimensionality != b.dimensionality:
            return False
        ratio = (a/b).to_base_units().magnitude
        return abs(ratio - 1.0) < 1e-9
    except (pint.errors.PintError, ValueError, AttributeError, TypeError):
        return False

_REGISTRY = None
def _registry():
    global _REGISTRY
    if _REGISTRY is None:
        import pint
        _REGISTRY = pint.UnitRegistry()
    return _REGISTRY

#===============================================================================
#===============================================================================
