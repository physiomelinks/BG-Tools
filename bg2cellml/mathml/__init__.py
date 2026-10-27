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
from copy import deepcopy
from dataclasses import dataclass

#===============================================================================

import lxml.etree as etree
import sympy
from sympy.printing.mathml import MathMLContentPrinter

#===============================================================================

from ..utils import etree_from_string, XMLNamespace

#===============================================================================

MATHML_NS = XMLNamespace('http://www.w3.org/1998/Math/MathML')

#===============================================================================
#===============================================================================

def sympy_from_mathml(node: etree.Element) -> sympy.Expr:
#========================================================
    if node.tag == MATHML_NS.bvar:
        return sympy_from_mathml(node.getchildren()[0])
    elif node.tag == MATHML_NS.ci:
        return sympy.Symbol(node.text)
    elif node.tag == MATHML_NS.cn:
        return sympy.Float(float(node.text))
    elif node.tag == MATHML_NS.pi:
        return sympy.pi
    elif node.tag != MATHML_NS.apply:
        raise ValueError(f'Expected MathML `apply`, got `{node.tag}`')
    children = node.getchildren()
    operator = children[0].tag
    if operator == MATHML_NS.times:
        terms = [sympy_from_mathml(child) for child in children[1:]
                    if not isinstance(child, etree._Comment)]
        return sympy.Mul(*terms)
    elif operator == MATHML_NS.divide:
        return sympy.Mul(sympy_from_mathml(children[1]),
                         sympy.Pow(sympy_from_mathml(children[2]), sympy.Integer(-1)))
    elif operator == MATHML_NS.plus:
        terms = [sympy_from_mathml(child) for child in children[1:]
                    if not isinstance(child, etree._Comment)]
        return sympy.Add(*terms)
    elif operator == MATHML_NS.minus:
        if len(children) == 2:
            return sympy.Mul(sympy.Integer(-1), sympy_from_mathml(children[1]))
        else:
            return sympy.Add(sympy_from_mathml(children[1]),
                             sympy.Mul(sympy.Integer(-1), sympy_from_mathml(children[2])))
    elif operator == MATHML_NS.ln:
        return sympy.log(sympy_from_mathml(children[1]))
    elif operator == MATHML_NS.exp:
        return sympy.exp(sympy_from_mathml(children[1]))
    elif operator == MATHML_NS.power:
        return sympy.Pow(sympy_from_mathml(children[1]), sympy_from_mathml(children[2]))
    elif operator == MATHML_NS.diff:
        return sympy.Derivative(sympy_from_mathml(children[2]), sympy_from_mathml(children[1]), evaluate=False)
    elif operator == MATHML_NS.cos:
        return sympy.cos(sympy_from_mathml(children[1]))
    elif operator == MATHML_NS.sin:
        return sympy.sin(sympy_from_mathml(children[1]))
    elif operator == MATHML_NS.abs:
        return sympy.Abs(sympy_from_mathml(children[1]))
    elif operator in (MATHML_NS.max, MATHML_NS.min):
        args = [sympy_from_mathml(child) for child in children[1:] if not isinstance(child, etree._Comment)]
        return (sympy.Max if operator == MATHML_NS.max else sympy.Min)(*args)
    elif operator == MATHML_NS.root:
        # <root/> with an optional <degree> (a square root without one)
        degree = sympy.Integer(2)
        args = [child for child in children[1:] if not isinstance(child, etree._Comment)]
        if args and args[0].tag == MATHML_NS.degree:
            degree = sympy_from_mathml(args[0].getchildren()[0])
            args = args[1:]
        return sympy.Pow(sympy_from_mathml(args[0]), sympy.Integer(1)/degree)
    else:
        raise ValueError(f'Unsupported MathML operator: {operator}')

#===============================================================================

sympy_mathml_printer = MathMLContentPrinter({'disable_split_super_sub': True})

def sympy_to_mathml(sympy_object) -> str:
#========================================
    mathml = ['<math xmlns="http://www.w3.org/1998/Math/MathML">']
    mathml.append(sympy_mathml_printer.doprint(sympy_object))
    mathml.append('</math>')
    return ''.join(mathml)

#===============================================================================
#===============================================================================

@dataclass
class Equation:
    lhs: sympy.Symbol | sympy.Derivative
    rhs: sympy.Basic
    provenance: str|None = None

    def __str__(self):
        return f'{self.lhs} = {self.rhs}'

    def sympy_equation(self):
    #========================
        return sympy.Eq(self.lhs, self.rhs)

    def mathml_equation(self) -> etree.Element:
    #==========================================
        return etree.fromstring(
            sympy_to_mathml(
                self.sympy_equation()))

#===============================================================================
#===============================================================================

class MathML:
    def __init__(self, mathml: etree.Element):
        self.__mathml = mathml
        self.__variables = defaultdict(list)
        for element in self.__mathml.findall(f'.//{MATHML_NS.ci}'):
            self.__variables[element.text].append(element)
        self.__equations: list[Equation] = []
        self.__build_equations()

    @classmethod
    def from_string(cls, formulae: str) -> 'MathML':
        return cls(etree_from_string(formulae))

    @property
    def equations(self) -> list[Equation]:
        return self.__equations

    @property
    def variables(self) -> list[str]:
        return list(self.__variables.keys())

    def __str__(self):
    #=================
        return etree.tostring(self.__mathml, encoding='unicode', pretty_print=True)

    def copy(self) -> 'MathML':
    #==========================
        return MathML(deepcopy(self.__mathml))

    def scale(self, symbol: str, factor: str, exponent: int=1):
    #===========================================================
        '''Replace each use of `symbol` by `symbol*factor^exponent` (`symbol/factor` when the exponent is -1).'''
        if symbol not in self.__variables:
            raise ValueError(f'Variable {symbol} not in formulae, cannot scale it')
        for element in self.__variables[symbol]:
            product = etree.Element(MATHML_NS.apply)
            etree.SubElement(product, MATHML_NS.divide if exponent == -1 else MATHML_NS.times)
            element.addprevious(product)
            product.append(element)         # moves the variable into the product
            signal = etree.SubElement(product, MATHML_NS.ci)
            signal.text = factor
            if exponent not in (1, -1):
                power = etree.Element(MATHML_NS.apply)
                etree.SubElement(power, MATHML_NS.power)
                product.replace(signal, power)
                power.append(signal)
                etree.SubElement(power, MATHML_NS.cn).text = str(exponent)
        self.__variables = defaultdict(list)
        for element in self.__mathml.findall(f'.//{MATHML_NS.ci}'):
            self.__variables[element.text].append(element)
        self.__build_equations()

    def substitute(self, name: str, symbol: str, missing_ok=False):
    #==============================================================
        if name not in self.__variables:
            if missing_ok:
                return
            raise ValueError(f'Variable {name} not in formulae, cannot substitute it')
        elif name == symbol:
            return
        elif symbol in self.__variables:
            raise ValueError(f'Symbol {symbol} is already in formulae, cannot substitute to it')
        for element in self.__variables[name]:
            element.text = symbol
        self.__variables[symbol] = self.__variables[name]
        del self.__variables[name]
        self.__build_equations()

    def __build_equations(self):
    #===========================
        self.__equations: list[Equation] = []
        for equation in self.__mathml.findall(f'.//{MATHML_NS.apply}/{MATHML_NS.eq}'):
            lhs = sympy_from_mathml(equation.getnext())
            rhs = sympy_from_mathml(equation.getnext().getnext())
            if isinstance(rhs, sympy.Derivative):
                if isinstance(lhs, sympy.Derivative):
                    raise ValueError(f'Unsupported form of ODE equation: {equation.xml}')
                self.__equations.append(Equation(rhs, lhs))
            elif isinstance(lhs, sympy.Symbol) or isinstance(lhs, sympy.Derivative):
                self.__equations.append(Equation(lhs, rhs))
            else:
                raise ValueError(f'Unsupported form of equation: {equation.xml}')

#===============================================================================
