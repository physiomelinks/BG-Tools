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
A model's variables named its own way in the CellML (BG-RDF's bgf:NamingRule and bgf:VariableName).

By default a variable is named after its element: ``q_G6P_c``, ``u_q_G6P_c``, ``v_Re_HK_c``. A model may give
rules, each naming one variable (``bgf:forVariable``: the template's q, u, v, u_f, K, kappa, ...) of the elements
of one class (``bgf:forElementClass``, e.g. bgf:QuantityStore; every class when absent) in one physical domain
(``bgf:forDomain``, e.g. bgf:Chemical; every domain when absent) by a pattern (``bgf:namePattern``), and an element
may name one of its variables explicitly (``bgf:hasVariableName``). An explicit name wins over rules; of the rules,
the one for the element's class and domain wins, then its class, then its domain, then neither.

A pattern's placeholders: ``{species}``, ``{location}``, ``{name}`` (the species, else the element's identifier),
``{element}`` (the bond graph's name for the element, e.g. q_G6P_c), ``{label}`` (its symbol, e.g. q) and ``{var}``.
A name keeps letters, digits and _ (any other character becomes _), drops repeated and surrounding _, starts with
``{var}_`` when it would start with a digit, and gets _2, _3, ... when it is already taken. The CellML editor's
Variables panel writes the same rules (src/plugins/bondgraph/naming.ts).
"""

#===============================================================================

import re
from typing import TYPE_CHECKING

import lxml.etree as etree

if TYPE_CHECKING:
    from .bondgraph import BondgraphModel

#===============================================================================

NAMING_RULES = """
    SELECT DISTINCT ?rule ?class ?domain ?variable ?pattern
    WHERE {
        <%MODEL%> bgf:hasNamingRule ?rule .
        ?rule bgf:forVariable ?variable ;
              bgf:namePattern ?pattern .
        OPTIONAL { ?rule bgf:forElementClass ?class }
        OPTIONAL { ?rule bgf:forDomain ?domain }
    }"""

VARIABLE_NAMES = """
    SELECT DISTINCT ?element ?variable ?name
    WHERE {
        <%MODEL%> bgf:hasBondElement ?element .
        ?element bgf:hasVariableName ?given .
        ?given bgf:forVariable ?variable ;
               bgf:cellmlName ?name .
    }"""

ELEMENTS = """
    SELECT DISTINCT ?uri ?symbol ?species ?location
    WHERE {
        <%MODEL%> bgf:hasBondElement ?uri .
        OPTIONAL { ?uri bgf:hasSymbol ?symbol }
        OPTIONAL { ?uri bgf:hasSpecies ?species }
        OPTIONAL { ?uri bgf:hasLocation ?location }
    }"""

CELLML_NS = 'http://www.cellml.org/cellml/1.1#'
MATHML_NS = 'http://www.w3.org/1998/Math/MathML'

#===============================================================================

def valid_name(name: str, var: str) -> str:
    """A CellML name: letters, digits and _, no repeated or surrounding _, not starting with a digit."""
    name = re.sub(r'_+', '_', re.sub(r'[^A-Za-z0-9_]', '_', name)).strip('_')
    if not name:
        return var
    if name[0].isdigit():
        name = f'{var}_{name}'
    return name

def render(pattern: str, *, var: str, element: str, label: str, species: str, location: str, name: str) -> str:
    """A pattern with its placeholders filled; an unknown placeholder is left out."""
    values = {'var': var, 'element': element, 'label': label, 'species': species, 'location': location, 'name': name}
    return valid_name(re.sub(r'\{(\w+)\}', lambda m: values.get(m.group(1), ''), pattern), var)

def _value(row, key: str) -> str:
    term = row.get(key)
    return term.value if term is not None else ''

def has_naming(model: 'BondgraphModel') -> bool:
    uri = model.uri.value
    return bool(list(model.rdf_graph.query(NAMING_RULES.replace('%MODEL%', uri)))
             or list(model.rdf_graph.query(VARIABLE_NAMES.replace('%MODEL%', uri))))

def naming_rules(model: 'BondgraphModel') -> list[dict[str, str]]:
    """The model's rules: {class, domain, variable, pattern}, class and domain '' when the rule names none."""
    return [{'class': _value(row, 'class'), 'domain': _value(row, 'domain'), 'variable': _value(row, 'variable'),
             'pattern': _value(row, 'pattern')}
            for row in model.rdf_graph.query(NAMING_RULES.replace('%MODEL%', model.uri.value))]

def rule_for(rules: list[dict[str, str]], element_class: str, domain: str, var: str) -> str|None:
    """The pattern of the most specific rule for the variable: class and domain, then class, then domain, then neither."""
    fits = [r for r in rules if r['variable'] == var and r['class'] in ('', element_class) and r['domain'] in ('', domain)]
    if not fits:
        return None
    return max(fits, key=lambda r: (r['class'] != '', r['domain'] != ''))['pattern']

def variable_table(model: 'BondgraphModel') -> list[dict[str, str]]:
    """Every element variable: its element, class, template name, bond graph name and CellML name, and the
    element's species, location and identifier (name)."""
    uri = model.uri.value
    rules = naming_rules(model)
    given: dict[tuple[str, str], str] = {}
    for row in model.rdf_graph.query(VARIABLE_NAMES.replace('%MODEL%', uri)):
        given[(_value(row, 'element'), _value(row, 'variable'))] = _value(row, 'name')
    context: dict[str, dict[str, str]] = {}
    for row in model.rdf_graph.query(ELEMENTS.replace('%MODEL%', uri)):
        context[row['uri'].value] = {'label': _value(row, 'symbol'), 'species': _value(row, 'species'),
                                     'location': _value(row, 'location')}

    # names not of elements' variables (time, the domains' constants) are taken
    taken: set[str] = {'t'}
    for element in model.elements:
        taken.update(constant.symbol for constant in element.domain.constants)
    rows = []
    assigned: dict[str, str] = {}               # bond graph name -> CellML name, of variables already named
    # in the order of the elements' identifiers, so a name taken twice gets _2 on the later one, as the editor shows
    for element in sorted(model.elements, key=lambda e: e.uri.value.rsplit('#', 1)[-1]):
        ctx = context.get(element.uri.value, {'label': '', 'species': '', 'location': ''})
        identifier = element.uri.value.rsplit('#', 1)[-1]
        for var, variable in element.variables.items():
            old = variable.symbol
            if (element.uri.value, var) in given:
                wanted = valid_name(given[(element.uri.value, var)], var)
            elif (pattern := rule_for(rules, element.element_class, element.domain.id, var)) is not None:
                wanted = render(pattern, var=var, element=element.symbol or identifier, label=ctx['label'],
                                species=ctx['species'], location=ctx['location'], name=ctx['species'] or identifier)
            else:
                wanted = old
            if old in assigned and wanted == old:
                # a variable several elements share (e.g. a parameter a signal sets) keeps its one name
                new = assigned[old]
            else:
                new, n = wanted, 1
                while new in taken:
                    n += 1
                    new = f'{wanted}_{n}'
                taken.add(new)
                assigned.setdefault(old, new)
            rows.append({'element': element.uri.value, 'class': element.element_class.rsplit('#', 1)[-1],
                         'variable': var, 'bondgraph': old, 'cellml': new, 'species': ctx['species'],
                         'location': ctx['location'], 'name': identifier})
    return rows

def variable_names(model: 'BondgraphModel') -> dict[str, str]:
    """The renamed variables, bond graph name -> CellML name (empty when the model names none its own way)."""
    if not has_naming(model):
        return {}
    return {row['bondgraph']: row['cellml'] for row in variable_table(model) if row['bondgraph'] != row['cellml']}

def rename_cellml(root: etree._Element, names: dict[str, str]) -> etree._Element:
    """Rename variables in a CellML document (their declarations, their uses in MathML, initial values naming them)."""
    if not names:
        return root
    for variable in root.iter(f'{{{CELLML_NS}}}variable'):
        if (name := variable.get('name')) in names:
            variable.set('name', names[name])
        if (initial := variable.get('initial_value')) in names:
            variable.set('initial_value', names[initial])
    for ci in root.iter(f'{{{MATHML_NS}}}ci'):
        if ci.text is not None and ci.text.strip() in names:
            ci.text = names[ci.text.strip()]
    return root

#===============================================================================
#===============================================================================
