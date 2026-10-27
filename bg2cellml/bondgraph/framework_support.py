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

from dataclasses import dataclass
from typing import Self, TYPE_CHECKING

#===============================================================================

from ..mathml import MathML
from ..rdf import literal_as_string, uri_fragment
from ..rdf import isLiteral, Literal, literal, NamedNode, namedNode, RdfGraph
from ..units import Units, Value
from ..utils import Issue

from .namespaces import BGF, CDT
from .utils import clean_name, Labelled, ModelElement, optional_integer

#===============================================================================
#===============================================================================

FLOW_SOURCE        = BGF.FlowSource.value
FLOW_STORE         = BGF.FlowStore.value
POTENTIAL_SOURCE   = BGF.PotentialSource.value
QUANTITY_STORE     = BGF.QuantityStore.value
REACTION           = BGF.Reaction.value
RESISTANCE         = BGF.Resistance.value
SIGNAL_BLOCK       = BGF.SignalBlock.value

ONENODE_JUNCTION   = BGF.OneNode.value
TRANSFORM_JUNCTION = BGF.TransformNode.value  # Can act as a transformer of gyrator
TRANSFORMER_NODE   = BGF.TransformerNode.value  # A transform node that is a transformer, also between domains
GYRATOR_NODE       = BGF.GyratorNode.value      # ... a gyrator, also within a domain
ZERONODE_JUNCTION  = BGF.ZeroNode.value

#===============================================================================
#===============================================================================

# Variable of integration

VOI_SYMBOL = 't'
VOI_UCUMUNIT: Literal = literal('s', datatype=CDT.ucumunit) # pyright: ignore[reportArgumentType, reportAssignmentType]

#===============================================================================
#===============================================================================

TRANSFORM_RATIO_NAME     = 'RATIO'

TRANSFORM_FLOW_NAME      = 'v'
TRANSFORM_POTENTIAL_NAME = 'u'

TRANSFORM_PORT_IDS       = ['0', '1']

#===============================================================================

GYRATOR_EQUATIONS = MathML.from_string(f"""
    <math xmlns="http://www.w3.org/1998/Math/MathML">
        <apply>
            <eq/>
            <ci>{TRANSFORM_POTENTIAL_NAME}_{TRANSFORM_PORT_IDS[0]}</ci>
            <apply>
                <times/>
                <ci>{TRANSFORM_RATIO_NAME}</ci>
                <ci>{TRANSFORM_FLOW_NAME}_{TRANSFORM_PORT_IDS[1]}</ci>
            </apply>
        </apply>
        <apply>
            <eq/>
            <ci>{TRANSFORM_POTENTIAL_NAME}_{TRANSFORM_PORT_IDS[1]}</ci>
            <apply>
                <times/>
                <ci>{TRANSFORM_RATIO_NAME}</ci>
                <ci>{TRANSFORM_FLOW_NAME}_{TRANSFORM_PORT_IDS[0]}</ci>
            </apply>
        </apply>
    </math>
""")

#===============================================================================

TRANSFORMER_EQUATIONS = MathML.from_string(f"""
    <math xmlns="http://www.w3.org/1998/Math/MathML">
        <apply>
            <eq/>
            <ci>{TRANSFORM_POTENTIAL_NAME}_{TRANSFORM_PORT_IDS[1]}</ci>
            <apply>
                <times/>
                <ci>{TRANSFORM_RATIO_NAME}</ci>
                <ci>{TRANSFORM_POTENTIAL_NAME}_{TRANSFORM_PORT_IDS[0]}</ci>
            </apply>
        </apply>
        <apply>
            <eq/>
            <ci>{TRANSFORM_FLOW_NAME}_{TRANSFORM_PORT_IDS[1]}</ci>
            <apply>
                <divide/>
                <ci>{TRANSFORM_FLOW_NAME}_{TRANSFORM_PORT_IDS[0]}</ci>
                <ci>{TRANSFORM_RATIO_NAME}</ci>
            </apply>
        </apply>
    </math>
""")

#===============================================================================
#===============================================================================

class Variable:
    def __init__(self, element: Labelled|ModelElement|None, name: str,
                        units: Literal|Units|None=None,
                        value: Literal|None=None,
                        type: NamedNode|str|None=None,
                        units_of: str|None=None, free: bool=False):
        # A variable may have no units yet: those of another variable of its template (`units_of`), known once a
        # model's signals are connected, or (`free`) those of its value or of the signal connected to it
        self.__element = element
        self.__name = clean_name(name)
        self.__symbol: str|None = None
        self.__type = type.value if isinstance(type, NamedNode) else type
        self.__units_of = units_of
        self.__free = free
        if isLiteral(units):
            self.__units = Units.from_ucum(units)   # pyright: ignore[reportArgumentType]
        else:
            self.__units: Units|None = units   # pyright: ignore[reportAttributeAccessIssue]
        if value is not None:
            self.__value = Value.from_literal(value)
            if self.__units is None and self.__value.units is not None:
                self.__units = self.__value.units
        else:
            self.__value = None
        if self.__units is None:
            if units_of is None and not free:
                self.__report_issue(f'Variable {name} for {self.element_uri} has no Units specified')
            return

        if self.__value is not None:
            if self.__value.units is None:
                self.__value.set_units(self.__units)
            elif not self.__units.is_compatible_with(self.__value.units):
                self.__report_issue(f'Value for variable {name} has incompatible units ({self.__value.units} != {self.__units})')

    def __str__(self):
        return f'{self.symbol} ({self.__value if self.__value is not None else self.__units}) [{self.type}]'

    def __report_issue(self, issue: str):
        if isinstance(self.__element, ModelElement):
            self.__element.report_issue(issue)
        else:
            raise Issue(issue)

    @property
    def element_uri(self) -> NamedNode|None:
        return self.__element.uri if self.__element is not None else None

    @property
    def name(self):
        return self.__name

    @property
    def symbol(self) -> str:
        return self.__symbol if self.__symbol is not None else self.__name

    @property
    def type(self):
        return self.__type

    @property
    def value(self):
        return self.__value

    @property
    def units(self) -> Units:
        return self.__units         # type: ignore

    @property
    def units_known(self) -> bool:
        return self.__units is not None

    @property
    def units_of(self) -> str|None:
        return self.__units_of

    @property
    def free(self) -> bool:
        return self.__free

    def set_units(self, units: Units):
    #=================================
        self.__units = units

    def clear_value(self):
    #=====================
        self.__value = None

    def settle_value(self) -> str|None:
    #==================================
        # Once its units are known, a value given before them (or in other units) is expressed in them: a
        # dimensionless number is taken as it is, compatible units are converted. Returns an issue, if any.
        if self.__value is None or self.__units is None:
            return None
        value_units = self.__value.units
        number = self.__value.value
        if value_units is not None and value_units != self.__units and not value_units.is_dimensionless():
            if not value_units.is_compatible_with(self.__units):
                return f'Value of {self.__name} has units {value_units}, not compatible with {self.__units}'
            number = value_units.convert(number, self.__units)
        self.__value = Value(number, self.__units)
        return None

    def copy(self, suffix: str|None=None, domain: 'Domain|None'=None) -> 'Variable':
    #===============================================================================
        if suffix is None:
            name = self.__name
        else:
            suffix = clean_name(suffix)
            if domain is None:
                name = f'{self.__name}_{suffix}'
            else:
                domain_symbols = domain.intrinsic_symbols
                suffix_parts = suffix.split('_')
                if self.name in domain_symbols and self.name == suffix_parts[0]:
                    name = suffix
                else:
                    name = f'{self.__name}_{suffix}'
        copy = Variable(self.__element, name, units=self.__units, type=self.__type,
                        units_of=self.__units_of, free=self.__free)
        copy.__value = self.__value.copy() if self.__value is not None else None
        return copy

    def set_symbol(self, symbol: str):
    #=================================
        self.__symbol = symbol

    def set_value(self, value: Value):
    #=================================
        self.__value = value
        if self.__units is None:
            # a free variable takes its value's units; one with bgf:unitsOf settles its value once they are known
            if self.__units_of is None and value.units is not None:
                self.__units = value.units
            return
        if self.__value.units is None:
            self.__value.set_units(self.__units)                        # type: ignore
        elif not self.__units.is_compatible_with(self.__value.units):   # type: ignore
            self.__report_issue(
                f'Value for variable {self.__name} has incompatible units ({self.__units} != {self.__value.units})') # type: ignore

#===============================================================================

VOI_VARIABLE = Variable(None, VOI_SYMBOL, units=VOI_UCUMUNIT, type=BGF.VariableOfIntergration)

#===============================================================================
#===============================================================================

DOMAIN_CONSTANTS = """
    SELECT DISTINCT ?name ?value
    WHERE {
        <%DOMAIN_URI%>
            a bgf:PhysicalDomain ;
            bgf:hasConstant [
                bgf:varName ?name ;
                bgf:hasValue ?value
            ] .
    }"""

#===============================================================================

class Domain(Labelled):
    def __init__(self, id: str, label: str|None,
                    flow_name: str, flow_units: Literal,
                    potential_name: str, potential_units: Literal,
                    quantity_name: str, quantity_units: Literal):
        super().__init__(id, label)
        self.__flow = Variable(self, flow_name, units=flow_units, type=BGF.FlowVariable)
        self.__potential = Variable(self, potential_name, units=potential_units, type=BGF.PotentialVariable)
        self.__quantity = Variable(self, quantity_name, units=quantity_units, type=BGF.QuantityVariable)
        self.__intrinsic_symbols = [
            self.__flow.symbol,
            self.__potential.symbol,
            self.__quantity.symbol
        ]
        self.__constants: list[Variable] = []

    @classmethod
    def from_rdf_graph(cls, graph: RdfGraph,
                    uri: NamedNode, label: Literal|None,
                    flow_name: Literal, flow_units: Literal,
                    potential_name: Literal, potential_units: Literal,
                    quantity_name: Literal, quantity_units: Literal) -> Self:
        self = cls(uri.value, literal_as_string(label),
                        flow_name.value, flow_units,
                        potential_name.value, potential_units,
                        quantity_name.value, quantity_units)
        self.__add_constants(graph)
        return self

    def __eq__(self, other):
        return self.uri.value == other.uri.value

    @property
    def constants(self):
        return self.__constants

    @property
    def flow(self):
        return self.__flow

    @property
    def intrinsic_symbols(self):
        return self.__intrinsic_symbols

    @property
    def name(self) -> str:
        return self.symbol

    @property
    def potential(self):
        return self.__potential

    @property
    def quantity(self):
        return self.__quantity

    def __add_constants(self, graph: RdfGraph):
    #==========================================
        self.__constants.extend([Variable(self, row['name'].value, value=row['value'])  # pyright: ignore[reportArgumentType]
                                for row in graph.query(
                                    # ?name ?value
                                    DOMAIN_CONSTANTS.replace('%DOMAIN_URI%', self.uri.value))])

#===============================================================================

SIGNAL_DOMAIN_CONSTANTS = """
    SELECT DISTINCT ?name ?value
    WHERE {
        <%DOMAIN_URI%>
            a bgf:SignalDomain ;
            bgf:hasConstant [
                bgf:varName ?name ;
                bgf:hasValue ?value
            ] .
    }"""

class SignalDomain(Labelled):
    # The domain of signal blocks: no power variables, only constants (e.g. RT)
    def __init__(self, id: str, label: str|None):
        super().__init__(id, label)
        self.__constants: list[Variable] = []

    @classmethod
    def from_rdf_graph(cls, graph: RdfGraph, uri: NamedNode, label: Literal|None) -> Self:
        self = cls(uri.value, literal_as_string(label))
        self.__constants.extend([Variable(self, row['name'].value, value=row['value'])  # pyright: ignore[reportArgumentType]
                                for row in graph.query(SIGNAL_DOMAIN_CONSTANTS.replace('%DOMAIN_URI%', self.uri.value))])
        return self

    def __eq__(self, other):
        return self.uri.value == other.uri.value

    @property
    def constants(self):
        return self.__constants

    @property
    def intrinsic_symbols(self) -> list[str]:
        return []

#===============================================================================
#===============================================================================

@dataclass
class NamedPortVariable:
    name: str
    variable: Variable

    def __str__(self):
        return f'(name: {self.name}, variable: {self.variable})'

#===============================================================================

class PowerPort:
    def __init__(self, uri: NamedNode, flow: NamedPortVariable, potential: NamedPortVariable,
                       direction: NamedNode|None=None, bond_count: int=1):
        self.__uri = uri
        self.__flow = flow
        self.__potential = potential
        self.__direction = direction
        self.__bond_count = bond_count

    def __str__(self):
        dirn = uri_fragment(self.__direction) if self.__direction is not None else None
        return f'PP: {uri_fragment(self.__uri)}, potential: {self.__potential}, flow: {self.__flow}, dirn: {dirn}'

    @property
    def direction(self) -> NamedNode|None:
        return self.__direction

    @property
    def flow(self) -> NamedPortVariable:
        return self.__flow

    @property
    def potential(self) -> NamedPortVariable:
        return self.__potential

    def copy(self, suffix: str|None=None, domain: Domain|None=None) -> 'PowerPort':
    #========================================================================================
        return PowerPort(self.__uri,
            NamedPortVariable(name=self.__flow.name, variable=self.__flow.variable.copy(suffix=suffix, domain=domain)),
            NamedPortVariable(name=self.__potential.name, variable=self.__potential.variable.copy(suffix=suffix, domain=domain)),
            direction = self.__direction, bond_count=self.__bond_count)

#===============================================================================
#===============================================================================

TEMPLATE_PARAMETERS = """
    SELECT DISTINCT ?name ?units ?value ?unitsOf
    WHERE {
        <%TEMPLATE_URI%> bgf:hasParameter ?variable .
        ?variable bgf:varName ?name .
        OPTIONAL { ?variable bgf:hasUnits ?units }
        OPTIONAL { ?variable bgf:hasValue ?value }
        OPTIONAL { ?variable bgf:unitsOf ?unitsOf }
    }"""

TEMPLATE_VARIABLES = """
    SELECT DISTINCT ?name ?units ?value ?unitsOf
    WHERE {
        <%TEMPLATE_URI%> bgf:hasVariable ?variable .
        ?variable bgf:varName ?name .
        OPTIONAL { ?variable bgf:hasUnits ?units }
        OPTIONAL { ?variable bgf:hasValue ?value }
        OPTIONAL { ?variable bgf:unitsOf ?unitsOf }
    }"""

TEMPLATE_SIGNAL_INPUTS = """
    SELECT DISTINCT ?name ?units ?unitsOf
    WHERE {
        <%TEMPLATE_URI%> bgf:hasSignalInput ?variable .
        ?variable bgf:varName ?name .
        OPTIONAL { ?variable bgf:hasUnits ?units }
        OPTIONAL { ?variable bgf:unitsOf ?unitsOf }
    } ORDER BY ?name"""

TEMPLATE_SIGNAL_OUTPUT = """
    SELECT DISTINCT ?name ?units ?unitsOf
    WHERE {
        <%TEMPLATE_URI%> bgf:hasSignalOutput ?variable .
        ?variable bgf:varName ?name .
        OPTIONAL { ?variable bgf:hasUnits ?units }
        OPTIONAL { ?variable bgf:unitsOf ?unitsOf }
    }"""

#===============================================================================

TEMPLATE_PORT_BONDS = """
    SELECT DISTINCT ?portId ?direction ?bondCount
    WHERE {
        {
            { <%TEMPLATE_URI%> bgf:hasPort ?portId .
            }
        UNION {
            <%TEMPLATE_URI%> bgf:hasPort ?port .
            ?port
                bgf:portId ?portId ;
                bgf:direction ?direction .
            OPTIONAL { ?port bgf:bondCount ?bondCount }
            }
        }
    } ORDER BY ?portId"""

#===============================================================================

class ElementTemplate(Labelled):
    def __init__(self, id: str, element_class: NamedNode,
                    label: str|None, domain: 'Domain|SignalDomain', relation: str|Literal):
        super().__init__(id, label)
        self.__element_class = element_class.value
        self.__domain = domain
        if self.__element_class in [FLOW_SOURCE, POTENTIAL_SOURCE]:
            self.__relation = None
        else:
            mathml = None
            if isLiteral(relation):
                if relation.datatype.value == BGF.mathml.value: # pyright: ignore[reportAttributeAccessIssue]
                    mathml = relation.value                     # pyright: ignore[reportAttributeAccessIssue]
                else:
                    # Do we insist on datatyping? Default to MathML ??
                    mathml = relation
            if mathml is None:
                raise Issue(f'BondElement {id} has no constitutive relation')
            elif isLiteral(mathml):
                mathml = mathml.value                           # pyright: ignore[reportAttributeAccessIssue]
            try:
                self.__relation = MathML.from_string(mathml)    # pyright: ignore[reportArgumentType]
            except ValueError as error:
                raise Issue(f'{self.uri}: {error}')
        self.__power_ports: dict[str, PowerPort] = {}
        self.__parameters: dict[str, Variable] = {}
        self.__variables: dict[str, Variable] = {}
        self.__signal_inputs: dict[str, Variable] = {}
        self.__signal_output: Variable|None = None
        self.__intrinsic_variable: Variable|None = None

    @classmethod
    def from_rdf_graph(cls, graph: RdfGraph, uri: NamedNode, element_class: NamedNode,
                        label: Literal|None, domain: 'Domain|SignalDomain', relation: Literal) -> Self:
        self = cls(uri.value, element_class, literal_as_string(label), domain, relation)
        if self.__element_class != SIGNAL_BLOCK:
            self.__add_ports(graph)
        self.__add_variables(graph)
        self.__check_names()
        return self

    @property
    def constitutive_relation(self) -> MathML|None:
        return self.__relation

    @property
    def domain(self) -> Domain:
        return self.__domain        # pyright: ignore[reportReturnType]

    @property
    def is_signal_block(self) -> bool:
        return self.__element_class == SIGNAL_BLOCK

    @property
    def signal_inputs(self) -> dict[str, Variable]:
        '''The template's signal inputs, in order of their names.'''
        return self.__signal_inputs

    @property
    def signal_output(self) -> Variable|None:
        return self.__signal_output

    @property
    def element_class(self) -> str:
        return self.__element_class

    @property
    def intrinsic_variable(self) -> Variable|None:
        return self.__intrinsic_variable

    @property
    def parameters(self) -> dict[str, Variable]:
        return self.__parameters

    @property
    def power_ports(self) -> dict[str, PowerPort]:
        return self.__power_ports

    @property
    def variables(self) -> dict[str, Variable]:
        return self.__variables

    def __add_ports(self, graph: RdfGraph):
    #======================================
        for row in graph.query(TEMPLATE_PORT_BONDS.replace('%TEMPLATE_URI%', self.uri.value)):
            # ?portId ?direction ?bondCount
            if isLiteral(row['portId']):
                id = row['portId'].value
                suffix = f'_{id}'
                flow_var = self.__port_name_variable(self.domain.flow)
                potential_var = self.__port_name_variable(self.domain.potential, suffix)
                self.__power_ports[id] = PowerPort(namedNode(f'{self.uri.value}{suffix}'),  # pyright: ignore[reportArgumentType]
                                                    flow_var, potential_var,
                                                    direction=row.get('direction'),         # pyright: ignore[reportArgumentType]
                                                    bond_count=optional_integer(row.get('bondCount'), 1)) # pyright: ignore[reportArgumentType]
        if len(self.__power_ports) == 0:
            self.__power_ports = {'': PowerPort(self.uri,
                                    self.__port_name_variable(self.domain.flow),
                                    self.__port_name_variable(self.domain.potential)
                                 )}

    def __port_name_variable(self, domain_variable: Variable, suffix: str='') -> NamedPortVariable:
    #==============================================================================================
        port_var_name = f'{domain_variable.name}{suffix}'
        return NamedPortVariable(name=port_var_name,
                                variable=Variable(self, port_var_name, units=domain_variable.units, type=domain_variable.type))

    def __add_variables(self, graph: RdfGraph):
    #==========================================
        for row in graph.query(TEMPLATE_PARAMETERS.replace('%TEMPLATE_URI%', self.uri.value, True)):
            # ?name ?units ?value
            var_name = row['name'].value             # pyright: ignore[reportOptionalMemberAccess]
            if var_name in self.__domain.intrinsic_symbols:
                raise Issue(f'Cannot specify domain symbol {var_name} as a variable for {self.uri}')
            self.__parameters[var_name] = self.__template_variable(row)
        for row in graph.query(TEMPLATE_VARIABLES.replace('%TEMPLATE_URI%', self.uri.value, True)):
            # ?name ?units ?value
            var_name = row['name'].value             # pyright: ignore[reportOptionalMemberAccess]
            if var_name in self.__domain.intrinsic_symbols:
                raise Issue(f'Cannot specify domain symbol {var_name} as a variable for {self.uri}')
            self.__variables[var_name] = self.__template_variable(row)
        for row in graph.query(TEMPLATE_SIGNAL_INPUTS.replace('%TEMPLATE_URI%', self.uri.value, True)):
            # ?name ?units ?unitsOf
            self.__signal_inputs[row['name'].value] = self.__template_variable(row)  # pyright: ignore[reportOptionalMemberAccess]
        for row in graph.query(TEMPLATE_SIGNAL_OUTPUT.replace('%TEMPLATE_URI%', self.uri.value, True)):
            # ?name ?units ?unitsOf
            self.__signal_output = self.__template_variable(row)
        # A variable that is intrinsic to the element's class
        # Values of intrinsic variables are set by bgf:hasValue
        if self.__element_class == QUANTITY_STORE:
            self.__intrinsic_variable = self.__domain.quantity.copy()
        if self.__element_class == FLOW_STORE:
            self.__intrinsic_variable = self.__domain.flow.copy()
        elif self.__element_class == POTENTIAL_SOURCE:
            self.__intrinsic_variable = self.__domain.potential.copy()
        elif self.__element_class == FLOW_SOURCE:
            self.__intrinsic_variable = self.__domain.flow.copy()

    def __template_variable(self, row) -> Variable:
    #===============================================
        units_of = row['unitsOf'].value if row.get('unitsOf') is not None else None
        return Variable(self, row['name'].value, units=row.get('units'), value=row.get('value'),
                        units_of=units_of, free=(row.get('units') is None and units_of is None))

    def __check_names(self):
    #=======================
        names = []
        def add_name(name: str, unique=True):
            if name not in names:
                names.append(name)
            elif unique:
                raise Issue(f'Duplicate name `{name}` for {self.uri}')
        for name in self.__parameters.keys():
            add_name(name)
        for name in self.__variables.keys():
            add_name(name)
        eqn_names = self.__relation.variables if self.__relation is not None else []
        if len(names) > len(eqn_names):
            raise Issue(f"{self.uri} has variables that are not in it's constitutive relation")
        for name in self.__signal_inputs.keys():
            add_name(name)
        if self.__signal_output is not None:
            add_name(self.__signal_output.name)
        for port in self.__power_ports.values():
            if port.flow is not None:
                add_name(port.flow.name, False)
            if port.potential is not None:
                add_name(port.potential.name, False)
        names.extend([c.name for c in self.__domain.constants])
        names.extend(self.__domain.intrinsic_symbols)
        names.append(VOI_VARIABLE.name)
        for name in eqn_names:
            if name not in names:
                raise Issue(f'Constitutive relation of {self.uri} has undeclared name {name}')

#===============================================================================
#===============================================================================

class JunctionStructure(Labelled):
    def __init__(self, id: str, label: str|None):
        super().__init__(id, label)

#===============================================================================
#===============================================================================

class CompositeElement(Labelled):
    def __init__(self, id: str, template: ElementTemplate, junction: JunctionStructure, label: str|None):
        super().__init__(id, label)
        self.__template = template
        self.__junction = junction

    @property
    def junction(self):
        return self.__junction

    @property
    def template(self):
        return self.__template

#===============================================================================

class CompositeTemplate(Labelled):
    def __init__(self, id: str, template: ElementTemplate, label: str|None):
        super().__init__(id, label)
        self.__template = template

    @property
    def template(self):
        return self.__template

#===============================================================================
#===============================================================================

type BondgraphElementTemplate = ElementTemplate | CompositeTemplate

#===============================================================================
#===============================================================================
