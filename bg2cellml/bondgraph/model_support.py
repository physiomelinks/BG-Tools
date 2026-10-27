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

import re
from typing import TYPE_CHECKING

#===============================================================================

import networkx as nx
import sympy

#===============================================================================

from ..mathml import Equation, MathML
from ..rdf import Literal, NamedNode, namedNode
from ..rdf import isLiteral, isNamedNode, literal_as_string, uri_fragment
from ..units import Units, Value

from .framework_support import BondgraphElementTemplate, CompositeTemplate, ElementTemplate
from .framework_support import Domain, NamedPortVariable, PowerPort, Variable
from .framework_support import ONENODE_JUNCTION, TRANSFORM_JUNCTION, ZERONODE_JUNCTION
from .framework_support import FLOW_SOURCE, POTENTIAL_SOURCE
from .framework_support import FLOW_STORE, QUANTITY_STORE, REACTION, RESISTANCE
from .framework_support import GYRATOR_EQUATIONS, GYRATOR_NODE, TRANSFORMER_EQUATIONS, TRANSFORMER_NODE
from .framework_support import TRANSFORM_FLOW_NAME, TRANSFORM_PORT_IDS, TRANSFORM_POTENTIAL_NAME, TRANSFORM_RATIO_NAME
from .namespaces import BGF
from .utils import ModelElement, pretty_name

if TYPE_CHECKING:
    from .model import BondgraphModel

#===============================================================================

def make_element_port_id(element_id: str, port_id: str) -> str:
#==============================================================
    if port_id in [None, '']:
        return element_id
    else:
        return f'{element_id}_{port_id}'

def flow_expression(node_dict: dict, model: 'BondgraphModel') -> sympy.Expr|None:
#================================================================================
    if ('power_port' in node_dict and 'element' in node_dict
    and (expr := node_dict['element'].flow_expression) is not None):
        return expr
    return flow_symbol(node_dict, model)

def flow_symbol(node_dict: dict, model: 'BondgraphModel') -> sympy.Symbol|None:
#==============================================================================
    if 'power_port' in node_dict:       # A BondElement or TransformNode's power port
        return sympy.Symbol(node_dict['power_port'].flow.variable.symbol)
    elif 'junction' in node_dict:
        junction: BondgraphJunction = node_dict['junction']
        if junction.type == ONENODE_JUNCTION:
            return sympy.Symbol(junction.variables[''].symbol)
        elif junction.type == ZERONODE_JUNCTION:
            model.report_issue(f'Adjacent Zero Nodes to junction {junction.uri.value} must be merged')
            return
    model.report_issue(f'Unexpected bond graph node, cannot get flow: {node_dict}')

def potential_expression(node_dict: dict, model: 'BondgraphModel') -> sympy.Expr|None:
#=====================================================================================
    if ('power_port' in node_dict and 'element' in node_dict
    and (expr := node_dict['element'].potential_expression) is not None):
        return expr
    return potential_symbol(node_dict, model)

def potential_symbol(node_dict: dict, model: 'BondgraphModel') -> sympy.Symbol|None:
#===================================================================================
    if 'power_port' in node_dict:       # A BondElement or TransformNode's power port
        return sympy.Symbol(node_dict['power_port'].potential.variable.symbol)
    elif 'junction' in node_dict:
        junction: BondgraphJunction = node_dict['junction']  # pyright: ignore[reportUnknownVariableType]
        if junction.type == ZERONODE_JUNCTION:
            return sympy.Symbol(junction.variables[''].symbol)
        elif junction.type == ONENODE_JUNCTION:
            model.report_issue(f'Adjacent One Nodes to junction {junction.uri.value} must be merged')
            return
    model.report_issue(f'Unexpected bond graph node, cannot get potential: {node_dict}')

def clean_latex(latex: str) -> str:
#==================================
    latex = re.sub(r'\\[^{]*{([^}]*)}', '_\\1', latex.strip())
    latex = latex.replace('^', '_')
    latex = re.sub(r'[^a-zA-Z_0-9]', '', latex)
    return re.sub(r'^_*', '', latex)

def make_symbolic_name(result_row: dict) -> str|None:
#====================================================
    symbol = literal_as_string(result_row.get('symbol'))
    if symbol is None:
        symbol = uri_fragment(result_row['uri'].value)
    species = literal_as_string(result_row.get('species'))
    location = literal_as_string(result_row.get('location'))
    if species is not None:
        symbol += f'_{clean_latex(species)}'
    if location is not None:
        symbol += f'_{clean_latex(location)}'
    return symbol

#===============================================================================
#===============================================================================

ELEMENT_PARAMETER_VALUES = """
    SELECT DISTINCT ?name ?value ?symbol
    WHERE {
        <%ELEMENT%> bgf:parameterValue ?variable .
        ?variable
            bgf:varName ?name ;
            bgf:hasValue ?value .
        OPTIONAL { ?variable bgf:hasSymbol ?symbol }
    }"""

ELEMENT_ONE_WAY = """
    SELECT DISTINCT ?oneWay
    WHERE {
        <%ELEMENT%> bgf:oneWay ?oneWay .
    }"""

ELEMENT_STATE_VALUE = """
    SELECT DISTINCT ?value
    WHERE {
        <%ELEMENT%> bgf:hasValue ?value .
    }"""

ELEMENT_VARIABLE_VALUES = """
    SELECT DISTINCT ?name ?value ?symbol
    WHERE {
        <%ELEMENT%> bgf:variableValue ?variable .
        ?variable
            bgf:varName ?name ;
            bgf:hasValue ?value .
        OPTIONAL { ?variable bgf:hasSymbol ?symbol }
    }"""

#===============================================================================

class VariableValue:
    def __init__(self, value: Literal|NamedNode, symbol: Literal|None):
        self.__value = value
        self.__symbol = symbol.value if symbol is not None else None

    @property
    def value(self):
        return self.__value

    @property
    def symbol(self):
        return self.__symbol

#===============================================================================

class BondgraphElement(ModelElement):
    def __init__(self,  model: 'BondgraphModel', id: str, template: BondgraphElementTemplate,
                        parameter_values: dict[str, VariableValue]|None=None,
                        variable_values: dict[str, VariableValue]|None=None,
                        domain_uri: NamedNode|None=None, value: Value|MathML|None=None,
                        symbol: str|None=None, label: str|None=None):
        super().__init__(model, id, symbol=symbol, label=label)
        element_name = pretty_name(self.symbol, template.uri)
        if isinstance(template, CompositeTemplate):
            element_template = template.template
            composite = True
        else:
            element_template = template
            composite = False
        if element_template.domain is None:
            model.report_issue(f'No modelling domain for element {id} with template {element_name}/{domain_uri}')
        elif domain_uri is not None and element_template.domain.uri.value != domain_uri.value:
            model.report_issue(f'Domain mismatch for element {id} with template {element_name}/{domain_uri}')

        self.__domain = element_template.domain
        if parameter_values is not None and self.__domain is not None:
            # a value for one of the domain's constants is the model's (BondgraphModel.constant_values)
            constants = {constant.name for constant in self.__domain.constants}
            parameter_values = {name: value for name, value in parameter_values.items() if name not in constants}
        self.__element_class = element_template.element_class
        self.__type = element_template.uri.value

        if self.__element_class in [FLOW_SOURCE, POTENTIAL_SOURCE]:
            if isinstance(value, MathML):
                self.__constitutive_relation = value
            else:
                self.__constitutive_relation = None
        else:
            if element_template.constitutive_relation is None:
                model.report_issue(f'Template {element_template.uri} for element {id} has no constitutive relation')
            else:
                self.__constitutive_relation = element_template.constitutive_relation.copy()

        self.__power_ports: dict[str, PowerPort] = {}
        for port_id, port in element_template.power_ports.items():
            self.__power_ports[make_element_port_id(self.uri.value, port_id)] = port.copy(suffix=self.symbol, domain=self.__domain)

        self.__flow: Variable|None = None
        self.__potential: Variable|None = None
        self.__flow_symbol = None
        self.__potential_symbol = None
        self.__flow_expression = None
        self.__potential_expression = None
        self.__equations = []

        self.__implied_junction = None
        if composite:
            if self.__element_class == QUANTITY_STORE:
                self.__implied_junction = ZERONODE_JUNCTION
            elif self.__element_class == FLOW_STORE:
                self.__implied_junction = ONENODE_JUNCTION
            elif self.__element_class in [RESISTANCE, REACTION]:
                self.__implied_junction = ONENODE_JUNCTION
        elif self.__element_class == POTENTIAL_SOURCE:
            self.__implied_junction = ZERONODE_JUNCTION
        elif self.__element_class == FLOW_SOURCE:
            self.__implied_junction = ONENODE_JUNCTION

        self.__variables: dict[str, Variable] = {}
        self.__port_variable_names = set()
        for port in self.__power_ports.values():
            # A port, by definition, always has flow and potential??
            # Then for a composite storage element, element_flow = -sum(connected_flows)
            #                "   dissapator   "        potential = -sum(connected_potentials)
            #
            self.__variables[port.flow.name] = port.flow.variable
            self.__port_variable_names.add(port.flow.name)
            self.__variables[port.potential.name] = port.potential.variable
            self.__port_variable_names.add(port.potential.name)
        self.__variables.update({name: variable.copy(suffix=self.symbol, domain=self.__domain)
                                for name, variable in element_template.parameters.items()})
        self.__variables.update({name: variable.copy(suffix=self.symbol, domain=self.__domain)
                                for name, variable in element_template.variables.items()})

        # Defer assignment until we have the full bondgraph
        self.__variable_values: dict[str, VariableValue] = {}
        self.__missing_parameters = [name for name in element_template.parameters.keys()
                                     if parameter_values is None or name not in parameter_values]
        if parameter_values is not None:
            self.__variable_values.update(parameter_values)
        if variable_values is not None:
            for name in variable_values.keys():
                if name not in name not in self.__variables:
                    model.report_issue(f'Unknown variable {name} for element {self.pretty_name}')
                    continue
            self.__variable_values.update(variable_values)

        # Signal inputs, set by signal bonds (bgf:signalInput), and modulations, a parameter scaled by a signal
        # (bgf:modulates): (parameter, signal, exponent)
        self.__signal_inputs: dict[str, Variable] = {name: variable.copy(suffix=self.symbol, domain=self.__domain)
                                                     for name, variable in element_template.signal_inputs.items()}
        self.__connected_inputs: dict[str, Variable] = {}
        self.__modulations: list[tuple[str, Variable, int]] = []
        self.__set_by_signal: dict[str, Variable] = {}      # parameters (or a source's value) a signal sets

        self.__intrinsic_name: str|None = None
        if (intrinsic_var := element_template.intrinsic_variable) is not None:
            self.__intrinsic_name = intrinsic_var.name
            self.__variables[intrinsic_var.name] = intrinsic_var.copy(suffix=self.symbol, domain=self.__domain)
            self.__intrinsic_variable = self.__variables[intrinsic_var.name]
            if self.__element_class == FLOW_SOURCE:
                port_var = self.__power_ports[self.uri.value].flow
                port_var.variable = self.__intrinsic_variable
            elif self.__element_class == POTENTIAL_SOURCE:
                port_var = self.__power_ports[self.uri.value].potential
                port_var.variable = self.__intrinsic_variable
            if value is not None and isinstance(value, Value):
                self.__intrinsic_variable.set_value(value)
        else:
            self.__intrinsic_variable = None

        if (variable := self.__variables.get(self.__domain.flow.symbol)) is not None:
            self.__flow = variable
            self.__flow_symbol = sympy.Symbol(variable.symbol)
        if (variable := self.__variables.get(self.__domain.potential.symbol)) is not None:
            self.__potential = variable
            self.__potential_symbol = sympy.Symbol(variable.symbol)

    @classmethod
    def for_model(cls, model: 'BondgraphModel', uri: NamedNode, template: BondgraphElementTemplate,
    #==============================================================================================
                    domain_uri: NamedNode|None, symbol: str|None, label: str|None):
        parameter_values: dict[str, VariableValue] = {row['name'].value:
                                                        VariableValue(row['value'], row.get('symbol'))  # pyright: ignore[reportArgumentType]
            for row in model.sparql_query(ELEMENT_PARAMETER_VALUES.replace('%ELEMENT%', uri.value))
            # ?name ?value ?symbol
        }
        variable_values: dict[str, VariableValue] = {row['name'].value:
                                                        VariableValue(row['value'], row.get('symbol'))  # pyright: ignore[reportArgumentType]
            for row in model.sparql_query(ELEMENT_VARIABLE_VALUES.replace('%ELEMENT%', uri.value))
            # ?name ?value ?symbol
        }
        value: Value|MathML|None = None
        for row in model.sparql_query(ELEMENT_STATE_VALUE.replace('%ELEMENT%', uri.value)):
            # ?value
            if isLiteral(row['value']):
                if row['value'].datatype.value == BGF.mathml.value:   # pyright: ignore[reportAttributeAccessIssue, reportOptionalMemberAccess]
                    value = MathML.from_string(row['value'].value)    # pyright: ignore[reportOptionalMemberAccess]
                else:
                    value = Value.from_literal(row['value'])          # pyright: ignore[reportArgumentType]
                break
        element = cls(model, uri.value, template, domain_uri=domain_uri,
                      parameter_values=parameter_values, variable_values=variable_values,
                      value=value, symbol=symbol, label=label)
        for row in model.sparql_query(ELEMENT_ONE_WAY.replace('%ELEMENT%', uri.value)):
            if str(row['oneWay'].value).lower() in ('true', '1'):                   # pyright: ignore[reportOptionalMemberAccess]
                element.make_one_way()
        return element

    def make_one_way(self):
    #======================
        """Pass flow only forwards: the constitutive relation's flow becomes max(v, 0)."""
        if (self.__element_class not in [RESISTANCE, REACTION] or self.__constitutive_relation is None
         or not self.__constitutive_relation.rectify(self.__domain.flow.symbol)):
            self.report_issue(f'Element {self.pretty_name} has no flow law to make one way')

    @property
    def constitutive_relation(self) -> MathML|None:
        return self.__constitutive_relation

    @property
    def domain(self) -> Domain:
        return self.__domain

    @property
    def element_class(self) -> str:
        return self.__element_class

    @property
    def signal_inputs(self) -> dict[str, Variable]:
        '''The template's signal inputs (as declared, with their units), in order of their names.'''
        return self.__signal_inputs

    @property
    def connected_inputs(self) -> dict[str, Variable]:
        '''The variables signal bonds bring to the signal inputs.'''
        return self.__connected_inputs

    @property
    def modulations(self) -> list[tuple[str, Variable, int]]:
        return self.__modulations

    def connect_signal(self, name: str, variable: Variable):
    #======================================================
        if name in self.__signal_inputs:
            self.__connected_inputs[name] = variable
        else:
            self.__set_by_signal[name] = variable

    @property
    def set_by_signal(self) -> dict[str, Variable]:
        '''The parameters (and a source's value, v or u) signals set, in place of a value.'''
        return self.__set_by_signal

    def can_set(self, name: str) -> bool:
    #====================================
        return self.can_modulate(name)

    def settable_units(self, name: str) -> Units|None:
    #=================================================
        variable = self.__variables.get(name)
        return variable.units if variable is not None and variable.units_known else None

    def add_modulation(self, name: str, signal: Variable, exponent: int):
    #====================================================================
        self.__modulations.append((name, signal, exponent))

    def can_modulate(self, name: str) -> bool:
    #=========================================
        return (name in self.__variables and name not in self.__port_variable_names) or name == self.__intrinsic_name

    def signal_variable(self, name: str|None) -> Variable|None:
    #==========================================================
        '''The element's variable a signal bond carries: the one named, else a store's quantity, a source's value,
        and a flow store's, resistance's or reaction's flow.'''
        if name is not None:
            return self.__variables.get(name)
        if self.__element_class in [QUANTITY_STORE, FLOW_SOURCE, POTENTIAL_SOURCE]:
            return self.__intrinsic_variable
        return self.__flow

    @property
    def flow(self) -> Variable|None:
        return self.__flow

    @property
    def flow_expression(self) -> sympy.Basic|None:
        return self.__flow_expression

    @property
    def equations(self) -> list[Equation]:
        return self.__equations

    @property
    def implied_junction(self) -> str|None:
        return self.__implied_junction

    @property
    def power_ports(self) -> dict[str, PowerPort]:
        return self.__power_ports

    @property
    def potential(self) -> Variable|None:
        return self.__potential

    @property
    def potential_expression(self) -> sympy.Basic|None:
        return self.__potential_expression

    @property
    def type(self) -> str:
        return self.__type

    @property
    def variables(self) -> dict[str, Variable]:
        return self.__variables

    def __modulate_source_value(self):
    #=================================
        # A flow or potential source whose value is modulated: its value v becomes v0 (a parameter holding the
        # source's value) times the signals
        if (self.__element_class not in [FLOW_SOURCE, POTENTIAL_SOURCE] or self.__intrinsic_name is None
         or not any(name == self.__intrinsic_name for name, _, _ in self.__modulations)):
            return
        if self.__constitutive_relation is not None:
            self.report_issue(f'Source {self.pretty_name} is given by an expression and cannot be modulated')
            return
        name = self.__intrinsic_name
        if name in self.__set_by_signal:
            self.report_issue(f'The value of source {self.pretty_name} is set by a signal and cannot be modulated')
            return
        base_name = f'{name}0'
        intrinsic = self.__intrinsic_variable
        base = Variable(self, f'{base_name}_{self.symbol}', units=intrinsic.units, type=intrinsic.type)
        if intrinsic.value is not None:
            base.set_value(intrinsic.value)
            intrinsic.clear_value()
        self.__variables[base_name] = base
        self.__constitutive_relation = MathML.from_string(
            f'<math xmlns="http://www.w3.org/1998/Math/MathML"><apply><eq/><ci>{name}</ci><ci>{base_name}</ci></apply></math>')
        self.__modulations = [(base_name if n == name else n, s, e) for n, s, e in self.__modulations]

    def __set_parameters_by_signals(self):
    #=====================================
        # A parameter set by a signal is the signal's variable; a source's value set by one is defined as it
        for name in self.__missing_parameters:
            if name not in self.__set_by_signal:
                self.report_issue(f'Missing value for parameter {name} of element {self.pretty_name}')
        for name, signal in self.__set_by_signal.items():
            self.__variable_values.pop(name, None)
            if name == self.__intrinsic_name and self.__element_class in [FLOW_SOURCE, POTENTIAL_SOURCE]:
                if self.__constitutive_relation is not None:
                    self.report_issue(f'Source {self.pretty_name} is given by an expression and cannot be set by a signal')
                    continue
                self.__intrinsic_variable.clear_value()
                self.__constitutive_relation = MathML.from_string(
                    f'<math xmlns="http://www.w3.org/1998/Math/MathML"><apply><eq/><ci>{name}</ci><ci>{name}_signal</ci></apply></math>')
                self.__connected_inputs[f'{name}_signal'] = signal
            else:
                self.__variables[name] = signal

    # Substitute variable symbols into the constitutive relation
    def assign_variables(self, bond_graph: nx.DiGraph):
    #==================================================
        self.__set_parameters_by_signals()

        ## Remove variables associated with any unconnected ports
        unused_ports: list[str] = []
        for port_id, port in self.__power_ports.items():
            if bond_graph.degree(port_id) == 0:     # pyright: ignore[reportCallIssue]
                unused_ports.append(port_id)
                if self.__element_class in [RESISTANCE, REACTION]:
                    del self.__variables[port.potential.name]
                    self.__port_variable_names.remove(port.potential.name)
        for port_id in unused_ports:
            del self.__power_ports[port_id]

        for var_name, var_value in self.__variable_values.items():
            if (variable := self.__variables.get(var_name)) is None:
                self.report_issue(f'Element {self.pretty_name} has unknown name {var_name} for {self.__type}')
                continue
            if var_value.symbol is not None:
                variable.set_symbol(var_value.symbol)   ## need symbol of value[1]'s element...
            if isLiteral(var_value.value):
                variable.set_value(Value.from_literal(var_value.value))  # pyright: ignore[reportArgumentType]
            elif isNamedNode(var_value.value):
                if var_value.value.value not in bond_graph:
                    self.report_issue(f'Value for {self.pretty_name} refers to unknown element: {var_value.value}')
                    continue
                elif (element := bond_graph.nodes[var_value.value.value].get('element')) is None:
                    self.report_issue(f'Value for {self.pretty_name} is not a bond element: {var_value.value}')
                    continue
                elif element.__intrinsic_variable is None:
                    self.report_issue(f'Value for {self.pretty_name} is an element with no intrinsic variable: {var_value.value}')
                    continue
                elif variable.units != element.__intrinsic_variable.units:
                    self.report_issue(f'Units incompatible for {self.pretty_name} value: {var_value.value}')
                    continue
                else:
                    self.__variables[var_name] = element.__intrinsic_variable

        for name in self.__signal_inputs:
            if name not in self.__connected_inputs:
                self.report_issue(f'Signal input {name} of element {self.pretty_name} has no signal bond')
        self.__modulate_source_value()

        if self.__constitutive_relation is not None:
            for name, variable in self.__variables.items():
                self.__constitutive_relation.substitute(name, variable.symbol,
                                                        missing_ok=(name in self.__port_variable_names))
            for name, variable in self.__connected_inputs.items():
                try:
                    self.__constitutive_relation.substitute(name, variable.symbol)
                except ValueError as error:
                    self.report_issue(f'Signal input {name} of {self.pretty_name}: {error}')
            for name, signal, exponent in self.__modulations:
                if (variable := self.__variables.get(name)) is None:
                    self.report_issue(f'Element {self.pretty_name} has no parameter {name} to modulate')
                    continue
                try:
                    self.__constitutive_relation.scale(variable.symbol, signal.symbol, exponent)
                except ValueError as error:
                    self.report_issue(f'Modulating {name} of {self.pretty_name}: {error}')
            for eqn in self.__constitutive_relation.equations:
                if eqn.lhs == self.__flow_symbol:
                    self.__flow_expression = eqn.rhs
                elif eqn.lhs == self.__potential_symbol:
                    self.__potential_expression = eqn.rhs

    def build_expressions(self, bond_graph: nx.DiGraph):
    #===================================================
        for port_id, port in self.__power_ports.items():
            bond_inputs = []
            bond_outputs = []
            port_inputs = []
            port_outputs = []

            def update_state_equalities(node: str, input) -> sympy.Symbol|None:
                node_dict = bond_graph.nodes[node]
                edge = (node, port_id) if input else (port_id, node)
                bond_count = bond_graph.edges[edge].get('bond_count', 1)
                forward_dirn = (port.direction
                            and (input and port.direction == BGF.InwardPort
                              or not input and port.direction == BGF.OutwardPort))
                if self.__implied_junction == ONENODE_JUNCTION:
                    if (expr := potential_expression(node_dict, self.model)) is not None:
                        if bond_count != 1:
                            expr = sympy.Mul(bond_count, expr)
                        if input: bond_inputs.append(expr)
                        else: bond_outputs.append(expr)
                    # a reaction's port potentials (also of a reaction with only one bond: an outflow)
                    if self.__element_class == REACTION:
                        if (symbol := potential_symbol(node_dict, self.model)) is not None:
                            if bond_count != 1:
                                symbol = sympy.Mul(bond_count, symbol)
                            if forward_dirn: port_inputs.append(symbol)
                            else: port_outputs.append(symbol)

                elif self.__implied_junction == ZERONODE_JUNCTION:
                    if (expr := flow_expression(bond_graph.nodes[node], self.model)) is not None:
                        if bond_count != 1:
                            expr = sympy.Mul(bond_count, expr)
                        if input: bond_inputs.append(expr)
                        else: bond_outputs.append(expr)
                    if len(self.__power_ports) > 1 and self.__element_class == QUANTITY_STORE:
                        if (symbol := flow_symbol(node_dict, self.model)) is not None:
                            if bond_count != 1:
                                symbol = sympy.Mul(bond_count, symbol)
                            if forward_dirn: port_inputs.append(symbol)
                            else: port_outputs.append(symbol)

            for node in bond_graph.predecessors(port_id):
                update_state_equalities(node, True)
            for node in bond_graph.successors(port_id):
                update_state_equalities(node, False)

            bond_expr = sympy.Add(*bond_inputs, sympy.Mul(-1, sympy.Add(*bond_outputs)))
            port_expr = None
            if len(port_inputs) or len(port_outputs):
                if len(port_outputs):
                    if len(port_inputs):
                        port_expr = sympy.Add(*port_inputs, sympy.Mul(-1, sympy.Add(*port_outputs)))
                    else:
                        port_expr = sympy.Mul(-1, sympy.Add(*port_inputs))
                elif len(port_inputs):
                    port_expr = sympy.Add(*port_inputs)

            if self.__implied_junction == ONENODE_JUNCTION:
                if self.__potential_symbol is not None:
                    self.__equations.append(Equation(self.__potential_symbol, bond_expr))
                if port_expr is not None:
                    potential = sympy.Symbol(port.potential.variable.symbol)
                    self.__equations.append(Equation(potential, port_expr))
            elif self.__implied_junction == ZERONODE_JUNCTION:
                if self.__flow_symbol is not None:
                    self.__equations.append(Equation(self.__flow_symbol, bond_expr))
                if port_expr is not None:
                    flow = sympy.Symbol(port.flow.variable.symbol)
                    self.__equations.append(Equation(flow, port_expr))

#===============================================================================
#===============================================================================

class SignalBlock(ModelElement):
    '''A signal block (bgf:SignalBlock): one output, named as the block, computed from its signal inputs and
    parameters by its template's relation. Its variables' units may be its signals' (see `resolve_units`).'''
    def __init__(self, model: 'BondgraphModel', id: str, template: ElementTemplate,
                        parameter_values: dict[str, VariableValue], value: Value|None,
                        symbol: str|None=None, label: str|None=None):
        super().__init__(model, id, symbol=symbol, label=label)
        self.__template = template
        self.__type = template.uri.value
        self.__relation = template.constitutive_relation.copy() if template.constitutive_relation is not None else None
        self.__parameters = {name: variable.copy(suffix=self.symbol)
                             for name, variable in template.parameters.items()}
        self.__inputs = {name: variable.copy(suffix=self.symbol)
                         for name, variable in template.signal_inputs.items()}
        self.__connected: dict[str, Variable] = {}
        output = template.signal_output
        if output is None:
            model.report_issue(f'Signal block template {template.uri.value} has no output')
            self.__output_name = 'y'
            self.__output = Variable(self, self.symbol, free=True)
        else:
            self.__output_name = output.name
            self.__output = Variable(self, self.symbol, units=output.units if output.units_known else None,
                                     units_of=output.units_of, free=output.free)
        if value is not None:
            self.__output.set_value(value)
        for name, variable in self.__parameters.items():
            if (given := parameter_values.get(name)) is None:
                model.report_issue(f'Missing value for parameter {name} of block {self.pretty_name}')
                continue
            if given.symbol is not None:
                variable.set_symbol(given.symbol)
            if isLiteral(given.value):
                variable.set_value(Value.from_literal(given.value))      # pyright: ignore[reportArgumentType]
            else:
                model.report_issue(f'Parameter {name} of block {self.pretty_name} must be a value')
        self.__equations: list[Equation] = []
        self.__modulations: list[tuple[str, Variable, int]] = []

    @classmethod
    def for_model(cls, model: 'BondgraphModel', uri: NamedNode, template: ElementTemplate,
    #=====================================================================================
                    symbol: str|None, label: str|None):
        parameter_values: dict[str, VariableValue] = {row['name'].value:
                                                        VariableValue(row['value'], row.get('symbol'))  # pyright: ignore[reportArgumentType]
            for row in model.sparql_query(ELEMENT_PARAMETER_VALUES.replace('%ELEMENT%', uri.value))
        }
        value: Value|None = None
        for row in model.sparql_query(ELEMENT_STATE_VALUE.replace('%ELEMENT%', uri.value)):
            if isLiteral(row['value']):
                value = Value.from_literal(row['value'])          # pyright: ignore[reportArgumentType]
                break
        return cls(model, uri.value, template, parameter_values, value, symbol=symbol, label=label)

    @property
    def domain(self):
        return self.__template.domain

    @property
    def equations(self) -> list[Equation]:
        return self.__equations

    @property
    def inputs(self) -> dict[str, Variable]:
        return self.__inputs

    @property
    def connected_inputs(self) -> dict[str, Variable]:
        return self.__connected

    @property
    def output(self) -> Variable:
        return self.__output

    @property
    def template(self) -> ElementTemplate:
        return self.__template

    @property
    def type(self) -> str:
        return self.__type

    @property
    def variables(self) -> dict[str, Variable]:
        '''The block's own variables, its parameters and output (its inputs are their sources').'''
        return {**self.__parameters, self.__output_name: self.__output}

    def connect_signal(self, name: str, variable: Variable):
    #======================================================
        self.__connected[name] = variable

    def can_modulate(self, name: str) -> bool:
    #=========================================
        return name in self.__parameters

    def __local(self, name: str) -> Variable|None:
    #=============================================
        if name == self.__output_name:
            return self.__output
        if name in self.__parameters:
            return self.__parameters[name]
        return self.__connected.get(name)

    def resolve_units(self) -> bool:
    #===============================
        '''Give variables whose units are another's (bgf:unitsOf) those units, either way round, and a source of
        unknown units connected to an input of fixed units the input's. Returns whether any units were found.'''
        progress = False
        declared = [(self.__output_name, self.__output.units_of)]
        declared += [(name, variable.units_of) for name, variable in self.__parameters.items()]
        declared += [(name, variable.units_of) for name, variable in self.__inputs.items()]
        for name, units_of in declared:
            if units_of is None:
                continue
            if '*' in units_of or '/' in units_of:
                # the units of a product or quotient of two of the block's variables, once both are known
                operator = '*' if '*' in units_of else '/'
                left, right = (self.__local(n.strip()) for n in units_of.split(operator, 1))
                target = self.__local(name)
                if (target is not None and not target.units_known and left is not None and right is not None
                and left.units_known and right.units_known):
                    target.set_units(left.units.times(right.units) if operator == '*' else left.units.over(right.units))
                    progress = True
                continue
            a, b = self.__local(name), self.__local(units_of)
            if a is None or b is None:
                continue
            if a.units_known and not b.units_known:
                b.set_units(a.units)
                progress = True
            elif b.units_known and not a.units_known:
                a.set_units(b.units)
                progress = True
        for name, variable in self.__inputs.items():
            source = self.__connected.get(name)
            if source is not None and variable.units_known and not source.units_known:
                source.set_units(variable.units)
                progress = True
        return progress

    def check_and_build(self):
    #=========================
        for name, variable in self.__inputs.items():
            if (source := self.__connected.get(name)) is None:
                self.report_issue(f'Signal input {name} of block {self.pretty_name} has no signal bond')
            elif variable.units_known and source.units_known and not source.units.same_scale(variable.units):
                self.report_issue(f'Signal input {name} of block {self.pretty_name} has units {source.units}, '
                                  f'not {variable.units}')
        for name, variable in self.variables.items():
            if not variable.units_known:
                self.report_issue(f'Units of {name} of block {self.pretty_name} are unknown: give its value units')
            elif (issue := variable.settle_value()) is not None:
                self.report_issue(f'Block {self.pretty_name}: {issue}')
        if self.__relation is None:
            return
        try:
            self.__relation.substitute(self.__output_name, self.__output.symbol, missing_ok=True)
            for name, variable in self.__parameters.items():
                self.__relation.substitute(name, variable.symbol)
            for name, variable in self.__connected.items():
                self.__relation.substitute(name, variable.symbol)
        except ValueError as error:
            self.report_issue(f'Block {self.pretty_name}: {error}')
            return
        for name, signal, exponent in self.__modulations:
            try:
                self.__relation.scale(self.__parameters[name].symbol, signal.symbol, exponent)
            except ValueError as error:
                self.report_issue(f'Modulating {name} of block {self.pretty_name}: {error}')
        self.__equations = self.__relation.equations

    @property
    def modulations(self) -> list[tuple[str, Variable, int]]:
        return self.__modulations

    def add_modulation(self, name: str, signal: Variable, exponent: int):
    #====================================================================
        self.__modulations.append((name, signal, exponent))

#===============================================================================
#===============================================================================

class BondgraphBond(ModelElement):
    def __init__(self, model: 'BondgraphModel', id: str,
                        source_id: str, target_id: str,
                        count: int|None=None, label: str|None=None):
        super().__init__(model, id, label=label)
        self.__source_id = source_id
        self.__target_id = target_id
        self.__bond_count = count if count is not None else 1

    @property
    def bond_count(self) -> int:
        return self.__bond_count

    @property
    def source_id(self) -> str|None:
        return self.__source_id

    @property
    def target_id(self) -> str|None:
        return self.__target_id

#===============================================================================
#===============================================================================

class BondgraphJunction(ModelElement):
    def __init__(self, model: 'BondgraphModel', id: str, type: NamedNode,
            value: Literal|None, symbol: str|None, label: str|None):
        super().__init__(model, id, symbol=symbol, label=label)
        self.__type = type.value
        self.__junction = model.framework.junction(self.__type)
        # A Transformer or Gyrator node is a transform node that says which it is (a plain one is a
        # transformer within a physical domain and a gyrator between two)
        self.__transform_kind: str|None = None
        if self.__type in (TRANSFORMER_NODE, GYRATOR_NODE):
            self.__transform_kind = self.__type
            self.__type = TRANSFORM_JUNCTION
        if self.__junction is None:
            model.report_issue(f'Unknown Junction {self.__type} for node {id}')
        self.__transform_relation: MathML|None = None
        self.__value = value
        self.__variables: dict[str, Variable] = {}
        self.__equations: list[Equation] = []

    @property
    def equations(self) -> list[Equation]:
        return self.__equations

    @property
    def type(self) -> str:
        return self.__type

    @property
    def variables(self) -> dict[str, Variable]:
        return self.__variables

    def __get_domain(self, attributes: dict) -> Domain|None:
    #=======================================================
        if (domain := attributes.get('domain')) is None:
            self.report_issue(f'Cannot find domain for junction {self.pretty_name}. Are there bonds to it?')
            return
        return domain

    def assign_node_variables(self, bond_graph: nx.DiGraph):
    #=======================================================
        assert self.__type != TRANSFORM_JUNCTION
        if (domain := self.__get_domain(bond_graph.nodes[self.uri.value])) is not None:
            if self.__type == ONENODE_JUNCTION:
                self.__variables[''] = Variable(self, self.symbol, units=domain.flow.units, value=self.__value, type=BGF.FlowVariable)
            elif self.__type == ZERONODE_JUNCTION:
                self.__variables[''] = Variable(self, self.symbol, units=domain.potential.units, value=self.__value, type=BGF.PotentialVariable)

    def assign_transform_variables(self, bond_graph: nx.DiGraph):
    #============================================================
        assert self.__type == TRANSFORM_JUNCTION
        domains = []
        graph = bond_graph.to_undirected(as_view=True)
        self.__variables[TRANSFORM_RATIO_NAME] = Variable(self, self.symbol, value=self.__value)
        for port_id in TRANSFORM_PORT_IDS:
            port_uri_id = make_element_port_id(self.uri.value, port_id)
            # Domains have been assigned to the transform junction's power ports as part of
            # `check_and_assign_domains_to_bond_network()` when building the model's graph
            if (domain := self.__get_domain(bond_graph.nodes[port_uri_id])) is None:
                self.report_issue(f'Cannot determine physical domain for port of {self.pretty_name}')
                return
            domains.append(domain)
            neighbours = list(graph.neighbors(port_uri_id))
            power_port = None
            flow_name = f'{TRANSFORM_FLOW_NAME}_{self.symbol}_{port_id}'
            flow_var_name = f'{TRANSFORM_FLOW_NAME}_{port_id}'
            potential_name = f'{TRANSFORM_POTENTIAL_NAME}_{self.symbol}_{port_id}'
            potential_var_name = f'{TRANSFORM_POTENTIAL_NAME}_{port_id}'
            if len(neighbours):
                neighbour = graph.nodes[neighbours[0]]
                junction_type = None
                variable = None
                if 'element' in neighbour:
                    # Should be a composite element
                    bond_element: BondgraphElement = neighbour['element']
                    junction_type = bond_element.implied_junction
                    if junction_type == ONENODE_JUNCTION:
                        variable = bond_element.flow
                    elif junction_type == ZERONODE_JUNCTION:
                        variable = bond_element.potential
                elif 'junction' in neighbour:
                    junction: BondgraphJunction = neighbour['junction']
                    junction_type = junction.type
                    variable = junction.variables['']
                if variable is not None:
                    port_uri = namedNode(port_uri_id)
                    if junction_type == ONENODE_JUNCTION:
                        self.__variables[flow_var_name] = variable
                        potential = Variable(self, potential_name, units=domain.potential.units)
                        self.__variables[potential_var_name] = potential
                        power_port = PowerPort(port_uri, NamedPortVariable(flow_name, variable),    # pyright: ignore[reportArgumentType]
                                                         NamedPortVariable(potential_name, potential))
                    elif junction_type == ZERONODE_JUNCTION:
                        flow = Variable(self, flow_name, units=domain.flow.units)
                        self.__variables[flow_var_name] = flow
                        self.__variables[potential_var_name] = variable
                        power_port = PowerPort(port_uri, NamedPortVariable(flow_name, flow),        # pyright: ignore[reportArgumentType]
                                                         NamedPortVariable(potential_name, variable))
            if power_port is not None:
                bond_graph.nodes[port_uri_id]['power_port'] = power_port
                bond_graph.nodes[port_uri_id]['port_type'] = TRANSFORM_JUNCTION
        if self.__transform_kind == TRANSFORMER_NODE \
        or (self.__transform_kind is None and domains[0] == domains[1]):
            self.__transform_relation = TRANSFORMER_EQUATIONS.copy()
        else:
            self.__transform_relation = GYRATOR_EQUATIONS.copy()

    def __own_lhs(self, equation: Equation) -> Equation:
    #===================================================
        # A transform node's equations are written for its port variables (u_1 = RATIO u_0, v_1 = v_0/RATIO);
        # where a port takes a connected element's variable (a store's potential, a reaction's flow), that
        # variable is already defined there, so the equation is solved for the node's own variable instead,
        # keeping the CellML causal (each variable defined once) whichever way the bonds are drawn.
        own = {sympy.Symbol(v.symbol) for v in self.__variables.values()   # its ports' potentials and flows
               if self.symbol and v.symbol.startswith((f'u_{self.symbol}_', f'v_{self.symbol}_'))}
        lhs = equation.lhs
        if not isinstance(lhs, sympy.Symbol) or lhs in own:
            return equation
        unknowns = [s for s in equation.rhs.free_symbols if s in own]
        if len(unknowns) != 1:
            return equation
        solution = sympy.solve(sympy.Eq(lhs, equation.rhs), unknowns[0])
        if len(solution) != 1:
            return equation
        return Equation(unknowns[0], solution[0], equation.provenance)

    def build_equations(self, bond_graph: nx.DiGraph) -> list[Equation]:
    #===================================================================
        ## is this where we multiply by bondCount??
        if self.__type == TRANSFORM_JUNCTION:
            if self.__transform_relation is not None:
                for name, variable in self.__variables.items():
                    self.__transform_relation.substitute(name, variable.symbol)
                self.__equations = [self.__own_lhs(eq) for eq in self.__transform_relation.equations]
        elif bond_graph.degree[self.uri.value] > 1:
            # we are connected to several nodes
            inputs = []
            outputs = []
            equation_lhs = []   # Use a list to pass value from `update_state_equalities` function
            equal_value: list[str] = []

            def update_state_equalities(node: str, input):
                node_dict = bond_graph.nodes[node]
                edge = (node, self.uri.value) if input else (self.uri.value, node)
                bond_count = bond_graph.edges[edge].get('bond_count', 1)
                if self.__type == ONENODE_JUNCTION:
                    if (symbol := potential_symbol(bond_graph.nodes[node], self.model)) is not None:
                        if bond_count != 1:
                            symbol = sympy.Mul(bond_count, symbol)
                        if input: inputs.append(symbol)
                        else: outputs.append(symbol)
                    if 'element' in node_dict and node_dict['element'].implied_junction is None:
                        equal_value.append(node_dict['power_port'].flow.variable.symbol)
                elif self.__type == ZERONODE_JUNCTION:
                    if (symbol := flow_symbol(bond_graph.nodes[node], self.model)) is not None:
                        if bond_count != 1:
                            symbol = sympy.Mul(bond_count, symbol)
                        if input: inputs.append(symbol)
                        else: outputs.append(symbol)
                    if 'element' in node_dict and node_dict['element'].implied_junction is None:
                        equal_value.append(node_dict['power_port'].potential.variable.symbol)
                if len(equation_lhs) == 0 and (port := node_dict.get('power_port')) is not None:
                    if self.__type == ONENODE_JUNCTION:
                        if node_dict.get('port_type') in [RESISTANCE, REACTION]:
                            equation_lhs.append(sympy.Symbol(port.potential.variable.symbol))
                    elif self.__type == ZERONODE_JUNCTION:
                        if (node_dict.get('port_type') == QUANTITY_STORE and 'element' in node_dict
                        and node_dict['element'].implied_junction is None):
                            equation_lhs.append(sympy.Symbol(port.flow.variable.symbol))

            for node in bond_graph.predecessors(self.uri.value):
                update_state_equalities(node, True)
            for node in bond_graph.successors(self.uri.value):
                update_state_equalities(node, False)

            if len(inputs) or len(outputs):
                lhs = equation_lhs[0] if len(equation_lhs) else None
                if len(outputs):
                    if lhs in inputs:
                        assert lhs is not None
                        inputs.remove(lhs)
                        if len(inputs):
                            self.__equations.append(Equation(lhs, sympy.Add(*outputs, sympy.Mul(-1, sympy.Add(*inputs)))))
                        else:
                            self.__equations.append(Equation(lhs, sympy.Add(*outputs)))
                    else:
                        if lhs in outputs:
                            assert lhs is not None
                            outputs.remove(lhs)
                        else:
                            lhs = outputs.pop()
                        if len(outputs):
                            self.__equations.append(Equation(lhs, sympy.Add(*inputs, sympy.Mul(-1, sympy.Add(*outputs)))))
                        else:
                            self.__equations.append(Equation(lhs, sympy.Add(*inputs)))
                elif len(inputs) > 1:
                    if lhs in inputs:
                        assert lhs is not None
                        inputs.remove(lhs)
                    else:
                        lhs = inputs.pop()
                    self.__equations.append(Equation(lhs, sympy.Mul(-1, sympy.Add(*inputs))))

            if len(equal_value):
                # The first junction variable represents the flow/potential of the node itself
                junction_symbol = sympy.Symbol(self.__variables[''].symbol)
                for value in equal_value:
                    # Filter out known equality between an implied junction and an element's port
                    if (symbol := sympy.Symbol(value)) != junction_symbol:
                        self.__equations.append(Equation(junction_symbol, sympy.Symbol(value)))

        return self.__equations

#===============================================================================
#===============================================================================
