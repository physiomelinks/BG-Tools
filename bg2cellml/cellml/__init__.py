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

from collections import defaultdict
import copy
from datetime import datetime, timezone
from pathlib import Path
import re
from typing import TYPE_CHECKING

#===============================================================================

import lxml.etree as etree
import sympy

#===============================================================================

from ..annotations import model_annotations, ModelAnnotations
from ..causality import assign_causality, CausalityError
from ..bondgraph.framework_support import QUANTITY_STORE, Variable, VOI_VARIABLE
from ..bondgraph.model_support import BondgraphElement, BondgraphJunction
from ..bondgraph.utils import clean_name, ModelElement
from ..mathml import Equation, MATHML_NS
from ..rdf import literal, namedNode, RdfGraph, Triple, uri_fragment
from ..rdf.namespace import Namespace, RDF, RDFS, XSD
from ..units import Units
from ..utils import XMLNamespace
from ..version import __version__

if TYPE_CHECKING:
    from ..bondgraph.model import BondgraphModel

#===============================================================================

CELLML_MODEL_URI = 'file:///CELLML_MODEL_URI'
CELLML_MODEL_NS = Namespace(f'{CELLML_MODEL_URI}#')

BQMODEL = Namespace('http://biomodels.net/model-qualifiers/')

METADATA_NAMESPACES = {
    'bgf': 'https://bg-rdf.org/ontologies/bondgraph-framework#',
    'bqmodel': str(BQMODEL),
    'rdfs': str(RDFS),
}

# OMEX Metadata 1.2 composite annotations (section 2.3.7): local resources are in the annotation file's own
# namespace (`local:`), physical properties are OPB terms
ANNOTATION_URI = 'file:///CELLML_ANNOTATION_URI'
LOCAL_NS = Namespace(f'{ANNOTATION_URI}#')

BGF = Namespace(METADATA_NAMESPACES['bgf'])
BQBIOL = Namespace('http://biomodels.net/biology-qualifiers/')
OPB = Namespace('https://identifiers.org/opb:')
SEMSIM = Namespace('http://bime.uw.edu/semsim/')

OPB_CHEMICAL_MOLAR_FLOW_RATE = OPB.OPB_00592
OPB_MOLAR_AMOUNT_OF_CHEMICAL = OPB.OPB_00425

ANNOTATION_NAMESPACES = {
    'bqbiol': str(BQBIOL),
    'local': str(LOCAL_NS),
    'opb': str(OPB),
    'semsim': str(SEMSIM),
}

def local_name(*parts: str) -> str:
    return '_'.join(re.sub(r'[^A-Za-z0-9_]', '_', part) for part in parts if part)

#===============================================================================

CELLML_NS = XMLNamespace('http://www.cellml.org/cellml/1.1#')
CMETA_NS = XMLNamespace('http://www.cellml.org/metadata/2.0#')

#===============================================================================

def cellml_element(tag: str, *args, **attributes) -> etree.Element:
#==================================================================
    if 'id' in attributes:
        attributes[CMETA_NS('id')] = attributes.pop('id')
    return etree.Element(CELLML_NS(tag), *args, **attributes)

def cellml_subelement(parent: etree.Element, tag: str, *args, **attributes) -> etree.Element:
#============================================================================================
    if 'id' in attributes:
        attributes[CMETA_NS('id')] = attributes.pop('id')
    return etree.SubElement(parent, CELLML_NS(tag), *args, **attributes)

def header_comment(source_uri: str) -> etree.Comment:
#====================================================
    utc = datetime.now(timezone.utc)
    return etree.Comment(f'''
  This CellML file was generated at {utc.isoformat()}

  by [BG-Tools](https://github.com/CellDL/BG-Tools), version {__version__}

  from {source_uri}
''')

#===============================================================================

CELLML_UNITS = [
    'ampere', 'farad', 'katal', 'lux', 'pascal', 'tesla',
    'becquerel', 'gram', 'kelvin', 'meter', 'radian', 'volt',
    'candela', 'gray', 'kilogram', 'metre', 'second', 'watt',
    'celsius', 'henry', 'liter', 'mole', 'siemens', 'weber',
    'coulomb', 'hertz', 'litre', 'newton', 'sievert',
    'dimensionless', 'joule', 'lumen', 'ohm', 'steradian',
]

#===============================================================================

DIMENSIONLESS_UNITS_NAME = 'dim'
CELLML_UNITS_ATTRIB = CELLML_NS('units')

DIMENSIONLESS_UNIT_DEFINITION = [
    f'<units name="{DIMENSIONLESS_UNITS_NAME}">',
    '<unit units="dimensionless"/>',
    '</units>'
]

#===============================================================================

def symbol_sort_key(symbol: str) -> str:
    return (symbol[2:] + symbol[0:2]) if symbol[0:2] in ['u_', 'v_'] else symbol

#===============================================================================

class CellMLComponent:
    def __init__(self, name: str, parent: etree.Element):
        self.__id = name
        self.__name = name
        self.__element = cellml_subelement(parent, 'component', name=name, id=self.__id)
        self.__bg_elements = defaultdict(list[str])

    @property
    def bg_elements(self):
    #=====================
        return self.__bg_elements

    @property
    def element(self):
    #=================
        return self.__element

    @property
    def id(self):
    #============
        return self.__id

    def add_dimensionless_attrib(self):
    #====================================
        for element in self.__element.findall(f'.//{MATHML_NS.cn}'):
            element.attrib[CELLML_UNITS_ATTRIB] = DIMENSIONLESS_UNITS_NAME

    def add_element(self, element: etree.Element):
    #=============================================
        self.__element.append(element)

    def add_variable(self, variable: 'CellMLVariable', bg_node_id: str|None=None):
    #=============================================================================
        self.__element.append(variable.get_element())
        if bg_node_id is not None:
            self.__bg_elements[bg_node_id].append(f'{self.__name}/{variable.symbol}')

#===============================================================================

class CellMLVariable:
    def __init__(self, component: CellMLComponent, variable: Variable):
        self.__id = f'{component.id}-{variable.symbol}'
        self.__symbol = variable.symbol
        self.__units = variable.units.name
        if variable.value is not None:
            self.__initial_value = variable.value.value
            value_units = variable.value.units
            if (value_units is not None and variable.units_known and value_units != variable.units
            and not value_units.is_dimensionless() and value_units.is_compatible_with(variable.units)):
                self.__initial_value = float(f'{value_units.convert(variable.value.value, variable.units):.15g}')
        else:
            self.__initial_value = None

    @property
    def id(self):
    #============
        return self.__id

    def get_element(self) -> etree.Element:
    #======================================
        element = cellml_element('variable', name=self.__symbol, units=self.__units, id=self.__id)
        if self.__initial_value is not None:
            element.attrib['initial_value'] = f'{self.__initial_value}'
        return element

    @property
    def symbol(self):
    #================
        return self.__symbol

#===============================================================================

class CellMLModel:
    def __init__(self, model: 'BondgraphModel'):
        name = uri_fragment(model.uri).rsplit('.')[0]
        self.__name = f'BG_{clean_name(name)}'
        self.__cellml = cellml_element('model', name=self.__name,
                                        nsmap={
                                            None: str(CELLML_NS),
                                            'cellml': str(CELLML_NS),
                                            'cmeta': str(CMETA_NS)
                                        })
        self.__components: list[CellMLComponent] = []
        self.__components.append(CellMLComponent('main', self.__cellml))
        self.__first_component_element = self.__components[0].element
        self.__current_component = self.__components[0]

        self.__model = model
        self.__elements = {element.uri.value: element for element in model.elements}
        self.__known_units: set[str] = set()
        self.__known_fixed: set[str] = set()
        self.__known_variables: dict[str, tuple[Variable, ModelElement]] = {}
        self.__variable_ids: dict[str, str] = {}        # bond graph symbol -> the CellML variable's cmeta:id
        self.__has_annotations = False

        self.__metadata = RdfGraph()
        self.__metadata.add(Triple(namedNode(CELLML_MODEL_URI), BQMODEL.isDescribedBy, model.uri))

        self.__add_unit_xml(DIMENSIONLESS_UNIT_DEFINITION)  ## Only if <cn> in MathML??
        self.__add_fixed(VOI_VARIABLE)       # only if VOI in some element's CR??

        for element in model.elements:
            self.__add_element(element)
        for block in model.signal_blocks:
            for constant in block.domain.constants:
                self.__add_fixed(constant)
            for variable in block.variables.values():
                if variable.units_known:
                    self.__add_variable(variable, block)
        for junction in model.junctions:
            self.__add_junction_variables(junction)
        self.__output_variable_definitions()
        self.__add_annotations(model_annotations(model))
        self.__equations_to_mathml()
        for component in self.__components:
            component.add_dimensionless_attrib()

    @property
    def name(self):
    #==============
        return self.__name

    def __add_element(self, element: BondgraphElement):
    #==================================================
        for constant in element.domain.constants:
            self.__add_fixed(constant)
        for variable in element.variables.values():
            self.__add_variable(variable, element)

    def __add_fixed(self, variable: Variable):
    #===========================================
        if variable.symbol not in self.__known_fixed:
            # a domain constant the model gives its own value
            variable = self.__model.constant_values.get(variable.symbol, variable)
            self.__add_units(variable.units)
            cellml_variable = CellMLVariable(self.__current_component, variable)
            self.__current_component.add_variable(cellml_variable)
            self.__known_fixed.add(variable.symbol)

    def __add_junction_variables(self, junction: BondgraphJunction):
    #===============================================================
        for variable in junction.variables.values():
            self.__add_variable(variable, junction)

    def __add_metadata(self, cellml_variable: CellMLVariable, variable: Variable, model_element: ModelElement):
    #==========================================================================================================
        if variable.type is not None:
            self.__metadata.add(Triple(CELLML_MODEL_NS(cellml_variable.id), BQMODEL.isDerivedFrom, namedNode(model_element.id)))
            self.__metadata.add(Triple(CELLML_MODEL_NS(cellml_variable.id), RDF.type, variable.type))
        if model_element.label:
            # This won't create duplicate statements
            self.__metadata.add(Triple(namedNode(model_element.id), RDFS.label, literal(model_element.label)))

    def __add_units(self, units: Units):
    #===================================
        elements = self.__elements_from_units(units)
        if len(elements):
            for element in elements:
                self.__add_unit_xml(element)

    def __add_unit_xml(self, unit_xml: list[str]):
    #=============================================
        if len(unit_xml):
            units_element = etree.fromstring(''.join(unit_xml))
            self.__first_component_element.addprevious(units_element)

    def __add_variable(self, variable: Variable, element: ModelElement):
    #==================================================================
        if variable.symbol not in self.__known_variables:
            # variables/component
            self.__known_variables[variable.symbol] = (variable, element)

    def __output_variable_definitions(self):
    #=======================================
        for symbol in sorted(self.__known_variables.keys(), key=symbol_sort_key):
            variable, model_element = self.__known_variables[symbol]
            self.__add_units(variable.units)
            cellml_variable = CellMLVariable(self.__current_component, variable)
            self.__current_component.add_variable(cellml_variable, model_element.id)
            self.__variable_ids[symbol] = cellml_variable.id
            self.__add_metadata(cellml_variable, variable, model_element)

    def __element_variable_id(self, element_uri: str, variable_type: str) -> str|None:
    #=================================================================================
        # the cmeta:id of the element's variable of the type (it stays when naming rules rename the variable)
        element = self.__elements.get(element_uri)
        if element is not None:
            for variable in element.variables.values():
                if variable.type is not None and variable.type.value == variable_type:
                    return self.__variable_ids.get(variable.symbol)
        return None

    def __add_annotations(self, annotations: ModelAnnotations):
    #===========================================================
        # OMEX Metadata 1.2 composite annotations, for the reactions and species that are annotated or have
        # annotated participants
        add = lambda s, p, o: self.__metadata.add(Triple(s, p, o))
        entities = set()
        for reaction in annotations.reactions:
            if not (reaction.annotated or any(annotations.entity(p).annotated for p in reaction.participants)):
                continue
            element = local_name(uri_fragment(reaction.element))
            process = LOCAL_NS(f'process_{element}')
            if (flux := self.__element_variable_id(reaction.element, BGF.FlowVariable.value)) is not None:
                add(CELLML_MODEL_NS(flux), BQBIOL.isVersionOf, OPB_CHEMICAL_MOLAR_FLOW_RATE)
                add(CELLML_MODEL_NS(flux), BQBIOL.isPropertyOf, process)
            for term in reaction.rhea + reaction.sbo:
                add(process, BQBIOL.isVersionOf, namedNode(term))
            for term in reaction.rhea_parts:
                add(process, BQBIOL.hasPart, namedNode(term))
            counts = {'left': 0, 'right': 0}
            for participant in reaction.participants:
                counts[participant.side] += 1
                role, predicate = (('source', SEMSIM.hasSourceParticipant) if participant.side == 'left'
                              else ('sink', SEMSIM.hasSinkParticipant))
                node = LOCAL_NS(f'{role}_{element}_{counts[participant.side]}')
                add(process, predicate, node)
                add(node, SEMSIM.hasMultiplier, literal(f'{float(participant.coefficient)}', XSD.double))
                add(node, SEMSIM.hasPhysicalEntityReference,
                    LOCAL_NS(local_name('entity', participant.species, participant.location)))
                entities.add((participant.species, participant.location))
            for k, mediator in enumerate(reaction.mediators, 1):
                node = LOCAL_NS(f'mediator_{element}' if len(reaction.mediators) == 1 else f'mediator_{element}_{k}')
                add(process, SEMSIM.hasMediatorParticipant, node)
                add(node, SEMSIM.hasPhysicalEntityReference, namedNode(mediator))
            self.__has_annotations = True
        for key, entity in annotations.entities.items():
            if not (entity.annotated or key in entities):
                continue
            node = LOCAL_NS(local_name('entity', entity.species, entity.location))
            for term in entity.chebi:
                add(node, BQBIOL('is'), namedNode(term))
            for term in entity.go:
                add(node, BQBIOL.isPartOf, namedNode(term))
            for element_uri in entity.elements:
                element = self.__elements.get(element_uri)
                if (element is not None and element.element_class == QUANTITY_STORE
                and (amount := self.__element_variable_id(element_uri, BGF.QuantityVariable.value)) is not None):
                    add(CELLML_MODEL_NS(amount), BQBIOL.isVersionOf, OPB_MOLAR_AMOUNT_OF_CHEMICAL)
                    add(CELLML_MODEL_NS(amount), BQBIOL.isPropertyOf, node)
            self.__has_annotations = True

    def __elements_from_units(self, units: Units) -> list[list[str]]:
    #================================================================
        result = []
        def elements_from_units(units: Units) -> list[str]:
            if units.name in self.__known_units or units.name in CELLML_UNITS:
                return []
            elements = []
            elements.append(f'<units xmlns="{CELLML_NS}" name="{units.name}">')
            for item in units.base_items():
                if item[0] not in self.__known_units:
                    item_elements = elements_from_units(Units(item[0]))
                    result.append(item_elements)
                name = Units.normalise_name(item[0])
                if item[1] == 0: elements.append(f'<unit units="{name}"/>')
                else: elements.append(f'<unit units="{name}" exponent="{item[1]}"/>')
            elements.append('</units>')
            self.__known_units.add(units.name)
            return elements
        result.append(elements_from_units(units))
        return result

    def __output_equations(self, equations: list[Equation], description: str):
    #=========================================================================
        if len(equations):
            self.__current_component.add_element(etree.Comment(f' {description}'))
            for equation in sorted(equations, key=lambda eq: str(eq.lhs)):
                self.__current_component.add_element(equation.mathml_equation())

    def __equations_to_mathml(self):
    #===============================
        equations = self.__model.equations
        # where the bond graph's causality differs from the equations' left-hand sides, solve them for it
        known = {symbol for symbol, (variable, _) in self.__known_variables.items() if variable.value is not None}
        known |= self.__known_fixed
        known |= {str(eq.lhs.expr) for eq in equations if isinstance(eq.lhs, sympy.Derivative)}
        try:
            equations = assign_causality(equations, known)
        except CausalityError as error:
            self.__model.report_issue(f'Causality: {error}')
        element_odes: list[Equation] = []
        element_algebraics: list[Equation] = []
        junction_algebraics: list[Equation] = []
        for equation in equations:
            if equation.provenance == 'cr':
                if isinstance(equation.lhs, sympy.Symbol):
                    element_algebraics.append(equation)
                elif isinstance(equation.lhs, sympy.Derivative):
                    element_odes.append(equation)
            elif equation.provenance == 'be':
                element_algebraics.append(equation)
            else:
                junction_algebraics.append(equation)
        self.__output_equations(element_odes, 'Element ODEs')
        self.__output_equations(element_algebraics, 'Element algebraics')
        self.__output_equations(junction_algebraics, 'Junction algebraics')

    def annotation(self) -> list[dict[str, dict[str, list[str]]]]:
    #=============================================================
        from ..variable_names import variable_names
        names = variable_names(self.__model)
        def renamed(path: str) -> str:
            component, _, name = path.rpartition('/')
            return f'{component}/{names.get(name, name)}' if component else names.get(name, name)
        bg_vars = defaultdict(list[str])
        for component in self.__components:
            for bg_element, variables in component.bg_elements.items():
                bg_vars[bg_element].extend(renamed(v) for v in variables)
        return [
            { id: { 'variables': vars } }
                for id, vars in bg_vars.items()
        ]

    async def metadata(self, cellml_uri: str|None, annotation_uri: str|None=None) -> str:
    #====================================================================================
        # `annotation_uri` is where the metadata is written, the base of its `local:` resources; by default the
        # CellML's URI with a .ttl suffix
        namespaces = METADATA_NAMESPACES.copy()
        if self.__has_annotations:
            namespaces.update(ANNOTATION_NAMESPACES)
        namespaces['model'] = str(CELLML_MODEL_NS)
        # the source document's namespace (a BG-RDF model's URI has a fragment, a CellDL one's none)
        namespaces['diagram'] = f'{self.__model.uri.value.split("#", 1)[0]}#'
        metadata = await self.__metadata.serialise(namespaces)
        if cellml_uri is not None:
            if annotation_uri is None:
                annotation_uri = f'{cellml_uri.removesuffix(".cellml")}.ttl'
            metadata = metadata.replace(CELLML_MODEL_URI, cellml_uri)
        if annotation_uri is not None:
            metadata = metadata.replace(ANNOTATION_URI, annotation_uri)
        return metadata

    def to_xml(self) -> str:
    #=======================
        # the model's own names for its variables (bgf:NamingRule, bgf:VariableName), if it gives any
        from ..variable_names import rename_cellml, variable_names
        cellml = rename_cellml(copy.deepcopy(self.__cellml), variable_names(self.__model))
        cellml.addprevious(header_comment(self.__model.uri.value))
        cellml_tree = etree.ElementTree(cellml)
        return etree.tostring(cellml_tree,
            encoding='unicode', inclusive_ns_prefixes=['cellml'],
            pretty_print=True)

#===============================================================================
#===============================================================================
