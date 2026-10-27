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

from typing import Sequence, TYPE_CHECKING

#===============================================================================

import networkx as nx

#===============================================================================

from ..cellml import CellMLModel
from ..mathml import Equation
from ..rdf import ResultRow, RdfGraph, Triple
from ..rdf import isBlankNode, isNamedNode, literal_as_string, NamedNode
from ..utils import Issue

from .framework_support import BondgraphElementTemplate, Variable
from .framework_support import Domain, REACTION, TRANSFORM_JUNCTION, TRANSFORM_PORT_IDS
from .model_support import BondgraphBond, BondgraphElement, BondgraphJunction, SignalBlock
from .model_support import make_element_port_id, make_symbolic_name
from .namespaces import BGF, NAMESPACES, get_curie
from .utils import clean_name, Labelled, optional_integer, pretty_name, pretty_uri

if TYPE_CHECKING:
    from .framework import BondgraphFramework

#===============================================================================
#===============================================================================

MODEL_ELEMENTS = """
    SELECT DISTINCT ?uri ?type ?label ?domain ?symbol ?species ?location
    WHERE {
        <%MODEL%> bgf:hasBondElement ?uri .
        ?uri a ?type .
        OPTIONAL { ?uri bgf:hasDomain ?domain }
        OPTIONAL { ?uri bgf:hasSymbol ?symbol }
        OPTIONAL { ?uri bgf:hasSpecies ?species }
        OPTIONAL { ?uri bgf:hasLocation ?location }
        OPTIONAL { ?uri rdfs:label ?label }
    } ORDER BY ?uri"""

MODEL_JUNCTIONS = """
    SELECT DISTINCT ?uri ?type ?label ?value ?symbol ?species ?location
    WHERE {
        <%MODEL%> bgf:hasJunctionStructure ?uri .
        ?uri a ?type .
        OPTIONAL { ?uri bgf:hasValue ?value }
        OPTIONAL { ?uri bgf:hasSymbol ?symbol }
        OPTIONAL { ?uri bgf:hasSpecies ?species }
        OPTIONAL { ?uri bgf:hasLocation ?location }
        OPTIONAL { ?uri rdfs:label ?label }
    }"""

MODEL_BONDS = """
    SELECT DISTINCT ?powerBond ?source ?target ?label ?bondCount
    WHERE {
        <%MODEL%> bgf:hasPowerBond ?powerBond .
        OPTIONAL { ?powerBond bgf:hasSource ?source }
        OPTIONAL { ?powerBond bgf:hasTarget ?target }
        OPTIONAL { ?powerBond rdfs:label ?label }
        OPTIONAL { ?powerBond bgf:bondCount ?bondCount }
        FILTER NOT EXISTS { ?powerBond a bgf:SignalBond }
    }"""

MODEL_BLOCKS = """
    SELECT DISTINCT ?uri ?type ?label ?symbol
    WHERE {
        <%MODEL%> bgf:hasSignalBlock ?uri .
        ?uri a ?type .
        OPTIONAL { ?uri bgf:hasSymbol ?symbol }
        OPTIONAL { ?uri rdfs:label ?label }
    } ORDER BY ?uri"""

SIGNAL_BONDS = """
    SELECT DISTINCT ?bond ?source ?target ?variable ?input ?modulates ?exponent
    WHERE {
        ?bond
            a bgf:SignalBond ;
            bgf:hasSource ?source ;
            bgf:hasTarget ?target .
        OPTIONAL { ?bond bgf:signalVariable ?variable }
        OPTIONAL { ?bond bgf:signalInput ?input }
        OPTIONAL { ?bond bgf:modulates ?modulates }
        OPTIONAL { ?bond bgf:modulationExponent ?exponent }
    } ORDER BY ?bond"""

#===============================================================================
#===============================================================================

BONDGRAPH_MODEL = """
    SELECT DISTINCT ?uri ?label
    WHERE {
        ?uri a bgf:BondgraphModel .
        OPTIONAL { ?uri rdfs:label ?label }
    }"""

#===============================================================================

BONDGRAPH_BONDS = """
    SELECT DISTINCT ?uri ?source ?target
    WHERE {
        ?uri
            bgf:hasSource ?source ;
            bgf:hasTarget ?target .
        FILTER NOT EXISTS { ?uri a bgf:SignalBond }
    }"""

#===============================================================================
#===============================================================================

class BondgraphModel(Labelled):   ## Component ??
    def __init__(self, framework: 'BondgraphFramework', base_iri: str, rdf_source: str, debug=False):
        self.__issues: list[Issue] = []
        self.__framework = framework
        self.__rdf_graph = RdfGraph(NAMESPACES)
        self.__debug = debug
        self.__elements: dict[str, BondgraphElement] = {}
        self.__junctions: dict[str, BondgraphJunction] = {}
        self.__blocks: dict[str, SignalBlock] = {}
        self.__symbolic_names: set[str] = set()
        self.__bonds: list[BondgraphBond] = []
        self.__graph = nx.DiGraph()
        (model_id, label) = self.__load_rdf(base_iri, rdf_source)
        if model_id is not None:
            super().__init__(model_id, label)   # pyright: ignore[reportArgumentType]
            self.__initialise()

    def report_issue(self, reason: str):
    #===================================
        self.__issues.append(Issue(reason))

    def __load_rdf(self, base_iri: str, rdf_source: str) -> tuple[str|None, str|None]:
    #=================================================================================
        try:
            self.__rdf_graph.load(base_iri, rdf_source)
        except Exception as e:
            self.report_issue(str(e))
            return (None, None)
        models: list[tuple[str, str|None]] = []
        for row in self.__rdf_graph.query(BONDGRAPH_MODEL):
            # ?uri ?label
            models.append((row['uri'].value, label.value if (label := row.get('label')) is not None else None)) # pyright: ignore[reportArgumentType]
        if len(models) == 0:
            self.report_issue('No BondgraphModels defined in RDF source')
            return (None, None)
        elif len(models) > 1:
            self.report_issue("Multiple BondgraphModels defined in RDF source -- a model's URI must be given")
            return (None, None)
        return models[0]

    def __initialise(self):
    #======================
        self.__generate_bonds()
        element_name: str|None = None
        element_type: NamedNode|None = None
        last_element_uri: str|None = None
        last_template: BondgraphElementTemplate|None = None
        symbol = None
        for row in self.__rdf_graph.query(MODEL_ELEMENTS.replace('%MODEL%', self.uri.value)):
            # ?uri ?type ?domain ?symbol ?species ?location ?label ORDER BY ?uri ?type
            if row['type'].value.startswith(NAMESPACES['bgf']):                     # pyright: ignore[reportOptionalMemberAccess]
                element_uri = row['uri'] .value                                     # pyright: ignore[reportAssignmentType]
                element_type = row['type']                                          # pyright: ignore[reportAssignmentType]
                assert element_type is not None
                if element_uri != last_element_uri:                                 # pyright: ignore[reportOptionalMemberAccess]
                    if last_element_uri is not None and last_template is None:
                        self.report_issue(f'{get_curie(element_type)} element {element_name} is not instantiated 1')
                    symbol = self.__unique_name(make_symbolic_name(row), element_uri)
                    if symbol:
                        element_name = pretty_name(symbol, element_uri)
                    else:
                        element_name = pretty_name('', element_uri)
                    element = None
                    last_template = None
                    last_element_uri = element_uri
                template = self.__framework.element_template(element_type, row.get('domain'))   # pyright: ignore[reportArgumentType]
                if template is not None:
                    if last_template is None:
                        element = BondgraphElement.for_model(self, row['uri'], template,    # pyright: ignore[reportArgumentType]
                                                             row.get('domain'),             # pyright: ignore[reportArgumentType]
                                                             symbol, literal_as_string(row.get('label')))   # pyright: ignore[reportArgumentType]
                        self.__elements[element.uri.value] = element
                        last_template = template
                    elif template != last_template:
                        self.report_issue(f'{get_curie(element_type)} element {element_name} has multiple instances')   # pyright: ignore[reportArgumentType]

        if last_element_uri is not None and last_template is None:
            assert element_type is not None
            self.report_issue(f'{get_curie(element_type)} element {element_name} is not instantiated 2')

        for row in self.__rdf_graph.query(MODEL_BLOCKS.replace('%MODEL%', self.uri.value)):
            # ?uri ?type ?label ?symbol
            block_uri = row['uri'].value                                                # pyright: ignore[reportOptionalMemberAccess]
            if block_uri in self.__blocks:
                continue
            if (template := self.__framework.block_template(row['type'])) is None:     # pyright: ignore[reportArgumentType]
                continue
            symbol = self.__unique_name(make_symbolic_name(row), block_uri)
            self.__blocks[block_uri] = SignalBlock.for_model(self, row['uri'], template,  # pyright: ignore[reportArgumentType]
                                                             symbol, literal_as_string(row.get('label')))  # pyright: ignore[reportArgumentType]

        if len(self.__elements) == 0 and len(self.__blocks) == 0:
            self.report_issue(f'Model {(pretty_uri(self.uri))} has no valid elements')
            return

        for row in self.__rdf_graph.query(MODEL_JUNCTIONS.replace('%MODEL%', self.uri.value)):
            # ?uri ?type ?value ?symbol ?species ?location ?label
            if row['type'].value.startswith(NAMESPACES['bgf']):                                 # pyright: ignore[reportOptionalMemberAccess]
                symbol = self.__unique_name(make_symbolic_name(row), row['uri'].value)
                junction = BondgraphJunction(self, row['uri'].value, row['type'], row.get('value'),     # pyright: ignore[reportArgumentType]
                                             symbol, literal_as_string(row.get('label')))               # pyright: ignore[reportArgumentType]
                self.__junctions[junction.uri.value] = junction

        for row in self.__rdf_graph.query(MODEL_BONDS.replace('%MODEL%', self.uri.value)):
            # ?powerBond ?source ?target ?label ?bondCount
            bond_id: NamedNode = row['powerBond'].value                                    # pyright: ignore[reportAssignmentType]
            source: NamedNode|None = row.get('source')                                     # pyright: ignore[reportAssignmentType]
            target: NamedNode|None = row.get('target')                                     # pyright: ignore[reportAssignmentType]
            if (source is None or isBlankNode(source)
             or target is None or isBlankNode(target)):
                self.report_issue(f'Bond {pretty_uri(bond_id)} is missing source and/or target node')
                continue
            source_id = source.value
            target_id = target.value
            issues = []
            if source_id not in self.__elements and source_id not in self.__junctions:
                issues.append(f'No element or junction for source {pretty_uri(source_id)} of bond {pretty_uri(bond_id)}')
            if target_id not in self.__elements and target_id not in self.__junctions:
                issues.append(f'No element or junction for target {pretty_uri(target_id)} of bond {pretty_uri(bond_id)}')
            if len(issues):
                for issue in issues:
                    self.report_issue(issue)
                continue

            # In terms of a BondgraphBond, the `target` is connected to the inwards port of the target node
            # and the `source` is connected to the outwards port of the source node.
            if (element := self.__elements.get(source_id)) is not None and element.element_class == REACTION:
                for port_id, port in element.power_ports.items():
                    if port.direction == BGF.OutwardPort:
                        source_id = port_id
                        break
            elif (junction := self.__junctions.get(source_id)) is not None and junction.type == TRANSFORM_JUNCTION:
                source_id = make_element_port_id(source_id, TRANSFORM_PORT_IDS[1])  # Outward port id
            if (element := self.__elements.get(target_id)) is not None and element.element_class == REACTION:
                for port_id, port in element.power_ports.items():
                    if port.direction == BGF.InwardPort:
                        target_id = port_id
                        break
            elif (junction := self.__junctions.get(target_id)) is not None and junction.type == TRANSFORM_JUNCTION:
                target_id = make_element_port_id(target_id, TRANSFORM_PORT_IDS[0])  # Inward port id
            self.__bonds.append(
                BondgraphBond(self, bond_id, source_id, target_id,
                              count=optional_integer(row.get('bondCount')),         # pyright: ignore[reportArgumentType]
                              label=literal_as_string(row.get('label'))))           # pyright: ignore[reportArgumentType]

        self.__make_bond_network()
        self.__check_and_assign_domains_to_bond_network()
        self.__assign_junction_domain_and_variables()
        self.__connect_signals()
        self.__assign_element_variables_and_equations()
        for block in self.__blocks.values():
            block.check_and_build()

        self.__equations: list[Equation] = []
        for element in self.__elements.values():
            if (cr := element.constitutive_relation) is not None:
                self.__equations.extend(cr.equations)
                for eq in cr.equations:
                    eq.provenance = 'cr'
            self.__equations.extend(element.equations)
            for eq in element.equations:
                eq.provenance = 'be'
        for junction in self.__junctions.values():
            equations = junction.build_equations(self.__graph)
            self.__equations.extend(equations)
            for eq in equations:
                eq.provenance = 'js'
        for block in self.__blocks.values():
            self.__equations.extend(block.equations)
            for eq in block.equations:
                eq.provenance = 'cr'

        if self.__debug:
            print('Elements:')
            for uri, element in self.__elements.items():
                print(' ', pretty_uri(uri))
                print('   Type:', pretty_uri(element.type))
                print('   Vars:')
                for (name, var) in element.variables.items():
                    print(f'     {name}: {var}')
                if (cr := element.constitutive_relation) is not None:
                    for eq in cr.equations:
                        print('   CR:', eq)
                for eq in element.equations:
                    print('   EQ:', eq)
            print('Junctions:')
            for uri, junction in self.__junctions.items():
                print(' ', pretty_uri(uri))
                print('   Type:', pretty_uri(junction.type))
                print('   Vars:')
                for var in junction.variables.values():
                    print(f'     {var}')
                equations = junction.equations
                for eq in equations:
                    print('   ', eq)

    def __unique_name(self, symbol: str|None, uri: str) -> str|None:
    #===============================================================
        # Variable names come from an element's symbol, species and location; two elements labelled
        # alike (e.g. copies of a chemostat) would share variables, so a later one takes its name
        # from its URI's fragment instead, numbered if that is taken too.
        if symbol is None:
            return None
        name = symbol
        if name in self.__symbolic_names:
            name = clean_name(uri.rsplit('#', 1)[-1])
            n = 1
            while name in self.__symbolic_names:
                n += 1
                name = f'{clean_name(uri.rsplit("#", 1)[-1])}_{n}'
        self.__symbolic_names.add(name)
        return name

    def __generate_bonds(self):
    #==========================
        for row in self.__rdf_graph.query(BONDGRAPH_BONDS):
            # ?bond ?source ?target
            source = row['source'] if isNamedNode(row['source']) else None
            target = row['target'] if isNamedNode(row['target']) else None
            if ((Triple(None, BGF.hasBondElement, source) in self.__rdf_graph
              or Triple(None, BGF.hasJunctionStructure, source) in self.__rdf_graph)
            and (Triple(None, BGF.hasBondElement, target) in self.__rdf_graph
              or Triple(None, BGF.hasJunctionStructure, target) in self.__rdf_graph)):
                self.__rdf_graph.add(Triple(self.uri, BGF.hasPowerBond, row['uri']))

    @property
    def elements(self) -> list[BondgraphElement]:
        return list(self.__elements.values())

    @property
    def signal_blocks(self) -> list[SignalBlock]:
        return list(self.__blocks.values())

    def __signal_source(self, source_id: str, variable_name: str|None) -> Variable|None:
    #==================================================================================
        if (block := self.__blocks.get(source_id)) is not None:
            if variable_name is not None and variable_name not in (block.template.signal_output.name    # pyright: ignore[reportOptionalMemberAccess]
                                                                   if block.template.signal_output else 'y',):
                self.report_issue(f'Signal from block {block.pretty_name} names {variable_name}, not its output')
            return block.output
        if (junction := self.__junctions.get(source_id)) is not None:
            if junction.type == TRANSFORM_JUNCTION:
                self.report_issue(f'A signal cannot come from transform node {junction.pretty_name}')
                return None
            return junction.variables.get('')
        if (element := self.__elements.get(source_id)) is not None:
            if (variable := element.signal_variable(variable_name)) is None:
                self.report_issue(f'Element {element.pretty_name} has no variable {variable_name} for a signal')
            return variable
        self.report_issue(f'Signal source {pretty_uri(source_id)} is not an element, junction or signal block')

    def __connect_signals(self):
    #===========================
        # Each signal bond sets a signal input of its target, or scales one of its parameters (bgf:modulates)
        targets: dict[str, BondgraphElement|SignalBlock] = {**self.__elements, **self.__blocks}
        for row in self.__rdf_graph.query(SIGNAL_BONDS):
            # ?bond ?source ?target ?variable ?input ?modulates ?exponent
            bond = row['bond'].value                                            # pyright: ignore[reportOptionalMemberAccess]
            if (target := targets.get(row['target'].value)) is None:            # pyright: ignore[reportOptionalMemberAccess]
                self.report_issue(f'Signal bond {pretty_uri(bond)} has no element or block as its target')
                continue
            variable_name = literal_as_string(row.get('variable'))             # pyright: ignore[reportArgumentType]
            if (signal := self.__signal_source(row['source'].value, variable_name)) is None:  # pyright: ignore[reportOptionalMemberAccess]
                continue
            if (modulates := literal_as_string(row.get('modulates'))) is not None:  # pyright: ignore[reportArgumentType]
                if not target.can_modulate(modulates):
                    self.report_issue(f'Signal bond {pretty_uri(bond)}: {target.pretty_name} has no parameter {modulates}')
                    continue
                exponent = optional_integer(row.get('exponent'), None)                       # pyright: ignore[reportArgumentType]
                if exponent is None and row.get('exponent') is not None:
                    exponent = int(float(row['exponent'].value))                             # pyright: ignore[reportOptionalMemberAccess]
                target.add_modulation(modulates, signal, exponent if exponent is not None else 1)
                continue
            inputs = target.inputs if isinstance(target, SignalBlock) else target.signal_inputs
            if (name := literal_as_string(row.get('input'))) is None:          # pyright: ignore[reportArgumentType]
                free = [n for n in inputs if n not in target.connected_inputs]
                if len(free) == 0:
                    self.report_issue(f'Signal bond {pretty_uri(bond)}: {target.pretty_name} has no free signal input')
                    continue
                name = free[0]
            elif name not in inputs:
                self.report_issue(f'Signal bond {pretty_uri(bond)}: {target.pretty_name} has no signal input {name}')
                continue
            elif name in target.connected_inputs:
                self.report_issue(f'Signal input {name} of {target.pretty_name} has two signal bonds')
                continue
            target.connect_signal(name, signal)
        # Units: a block's variables may have its signals' units
        while any([block.resolve_units() for block in self.__blocks.values()]):
            pass
        for element in self.__elements.values():
            for name, source in element.connected_inputs.items():
                declared = element.signal_inputs[name]
                if not source.units_known:
                    self.report_issue(f'Units of the signal to {name} of {element.pretty_name} are unknown')
                elif not source.units.same_scale(declared.units):
                    self.report_issue(f'Signal input {name} of {element.pretty_name} has units {source.units}, '
                                      f'not {declared.units}')
        for target in targets.values():
            for name, signal, _ in target.modulations:
                if signal.units_known and not signal.units.is_dimensionless():
                    self.report_issue(f'The signal modulating {name} of {target.pretty_name} has units '
                                      f'{signal.units}: it must be dimensionless')

    @property
    def framework(self) -> 'BondgraphFramework':
        return self.__framework

    @property
    def has_issues(self) -> bool:
        return len(self.__issues) > 0

    @property
    def issues(self) -> list[Issue]:
        return self.__issues

    @property
    def junctions(self) -> list[BondgraphJunction]:
        return list(self.__junctions.values())

    @property
    def equations(self):
        return self.__equations

    @property
    def network_graph(self) -> nx.DiGraph:
        return self.__graph.copy()

    @property
    def rdf_graph(self) -> RdfGraph:
        return self.__rdf_graph

    def __assign_element_variables_and_equations(self):
    #==================================================
        for element in self.__elements.values():
            element.assign_variables(self.__graph)
        for element in self.__elements.values():
            element.build_expressions(self.__graph)

    def __assign_junction_domain_and_variables(self):
    #================================================
        for junction in self.__junctions.values():
            if junction.type != TRANSFORM_JUNCTION:
                junction.assign_node_variables(self.__graph)
        for junction in self.__junctions.values():
            if junction.type == TRANSFORM_JUNCTION:
                junction.assign_transform_variables(self.__graph)

    # Assign junction domains from elements and check consistency
    def __check_and_assign_domains_to_bond_network(self):
    #====================================================
        seen_nodes = set()
        undirected_graph = self.__graph.to_undirected(as_view=True)

        def check_node(node: str, domain: Domain):
            if node not in seen_nodes:
                seen_nodes.add(node)
                if 'domain' not in self.__graph.nodes[node]:
                    self.__graph.nodes[node]['domain'] = domain
                    for neighbour in undirected_graph.neighbors(node):
                        check_node(neighbour, domain)
                elif domain != (node_domain := self.__graph.nodes[node]['domain']):
                    self.report_issue(f'Node {node} with domain {node_domain} incompatible with {domain}')
                    return

        for element in self.__elements.values():
            for port_id in element.power_ports.keys():
                check_node(port_id, element.domain)

    # Construct network graph of PowerBonds
    def __make_bond_network(self):
    #=============================
        for element in self.__elements.values():
            for port_id, port in element.power_ports.items():
                self.__graph.add_node(port_id, uri=port_id,
                    type=get_curie(element.type),
                    power_port=port, port_type=element.element_class,
                    element=element, label=element.symbol)
        for junction_id, junction in self.__junctions.items():
            if junction.type == TRANSFORM_JUNCTION:
                # A Transform Node has two implicit ports, with ids `0` and `1`
                for port_id in TRANSFORM_PORT_IDS:
                    port_uri_id = make_element_port_id(junction_id, port_id)
                    self.__graph.add_node(port_uri_id, uri=port_uri_id,
                        type=get_curie(junction.type), junction=junction, label=junction.symbol)
            else:
                self.__graph.add_node(junction_id, uri=junction_id,
                    type=get_curie(junction.type), junction=junction, label=junction.symbol)
        for bond in self.__bonds:
            source = bond.source_id
            target = bond.target_id
            if source not in self.__graph and target not in self.__graph:
                self.report_issue(f'No element or junction for source {pretty_uri(source)} and target {pretty_uri(target)} of bond {pretty_uri(bond.uri)}')
                continue
            elif source not in self.__graph:
                self.report_issue(f'No element or junction for source {pretty_uri(source)} of bond {pretty_uri(bond.uri)}')
                continue
            elif target not in self.__graph:
                self.report_issue(f'No element or junction for target {pretty_uri(target)} of bond {pretty_uri(bond.uri)}')
                continue

            self.__graph.add_edge(source, target, bond_count=bond.bond_count)

    def sparql_query(self, query: str) -> Sequence[ResultRow]:
    #========================================================
        return self.__rdf_graph.query(query)

#===============================================================================

    def make_cellml_model(self) -> CellMLModel:
    #==========================================
        return CellMLModel(self)

#===============================================================================
#===============================================================================
