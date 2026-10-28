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

from typing import Self

#===============================================================================

import re

import pint
from pint import Quantity
from ucumvert import PintUcumRegistry

#===============================================================================

from ..bondgraph.namespaces import CDT
from ..rdf import isLiteral, Literal, NamedNode

#===============================================================================

ucum_registry = PintUcumRegistry()

PREFERRED_BASE_ITEMS = {
    'kilopascal': [
        ('joule', 0),
        ('litre', -1),
    ]
}

SUBSTITUTIONS = {
    'kilopascal': 'kPa',
    'liter': 'litre',
    'meter': 'metre'
}

#===============================================================================

class Units:
    def __init__(self, units: str|pint.Unit):
        if isinstance(units, str):
            single_unit = ucum_registry(units)
            units = single_unit.u               # type: ignore
        self.__units: pint.Unit = units         # type: ignore
        self.__name = Units.normalise_name(str(self.__units))

    @classmethod
    def from_ucum(cls, ucum_units: Literal|str) -> Self:
        if isLiteral(ucum_units):
            if (ucum_units.datatype is not None                     # pyright: ignore[reportAttributeAccessIssue]
            and ucum_units.datatype.value != CDT.ucumunit.value):   # pyright: ignore[reportAttributeAccessIssue]
                raise TypeError(f'Units value has unexpected datatype: {ucum_units.datatype}')  # pyright: ignore[reportAttributeAccessIssue]
            return cls(ucum_registry.from_ucum(ucum_units.value).u) # pyright: ignore[reportAttributeAccessIssue]
        return cls(ucum_registry.from_ucum(ucum_units).u)

    @staticmethod
    def normalise_name(name: str) -> str:
    #====================================
        name = (name.replace(' * ', '_')
                    .replace(' / ', '_per_')
                    .replace(' ** 2', '_squared')
                    .replace(' ** 3', '_cubed')
                    .replace('1_per_', 'per_'))
        # any other power, e.g. second ** -1.5 -> second_pow_m1_5
        name = re.sub(r' \*\* (-?)([0-9.]+)', lambda m: f'_pow_{"m" if m[1] else ""}{m[2].replace(".", "_")}', name)
        for fullname, replacement in SUBSTITUTIONS.items():
            name = name.replace(fullname, replacement)
        return name

    def __eq__(self, other):
        return self.__units == other.__units

    def __str__(self):
        return str(self.__units)

    @property
    def name(self):
        return self.__name

    def base_items(self) -> list[tuple[str, int]]:
    #=============================================
        unit_quantity: Quantity = pint.Quantity(1, self.__units)
        return PREFERRED_BASE_ITEMS.get(str(self),
                                        unit_quantity.unit_items())
        ## if not PREFERRED_BASE_ITEMS start by going through
        ## self.__pint_units.unit_items() and only go to_base_units() if not a known
        ## CELLML_UNIT...

    def si_definition(self) -> tuple[float, list[tuple[str, int]]]|None:
    #===================================================================
        """
        A single unit that is not made of others (a prefixed or non-SI one: nanogram, millilitre, hour) as a
        multiplier of SI base units, e.g. (1e-12, [('kilogram', 1)]); None for any other units.
        """
        items = list(self.base_items())
        if len(items) != 1 or items[0][1] != 1 or items[0][0] != str(self.__units):
            return None
        base = pint.Quantity(1, self.__units).to_base_units()
        multiplier = float(f'{base.magnitude:.15g}')
        return (multiplier, [(Units.normalise_name(name), int(exponent)) for name, exponent in base.unit_items()])

    def times(self, other: 'Units') -> 'Units':
    #==========================================
        return Units(self.__units * other.__units)

    def over(self, other: 'Units') -> 'Units':
    #=========================================
        return Units(self.__units / other.__units)

    def is_compatible_with(self, other: Self):
    #=========================================
        return self.__units.is_compatible_with(other.__units)

    def is_dimensionless(self) -> bool:
    #==================================
        return bool(self.__units.dimensionless)

    def convert(self, value: float, other: 'Units') -> float:
    #=========================================================
        return float(pint.Quantity(value, self.__units).to(other.__units).magnitude)

    def same_scale(self, other: 'Units') -> bool:
    #=============================================
        # equal, or compatible with a factor of 1 (kPa and J/L)
        if self == other:
            return True
        return self.is_compatible_with(other) and abs(self.convert(1.0, other) - 1.0) < 1e-12

#===============================================================================

class Value:
    def __init__(self, value: float, units: Units|None):
        self.__units = units
        self.__value = value

    @classmethod
    def from_literal(cls, literal_value: Literal) -> Self:
        if literal_value.datatype.value == CDT.ucum.value:  # pyright: ignore[reportAttributeAccessIssue]
            parts = str(literal_value.value).split()
            value = float(parts[0])
            units = Units.from_ucum(parts[1])
        elif literal_value.datatype is None:                # pyright: ignore[reportAttributeAccessIssue]
            value = float(literal_value.value)
            units = None
        else:
            raise TypeError(f'Literal value has unexpected datatype: {literal_value.datatype}') # pyright: ignore[reportAttributeAccessIssue]
        return cls(value, units)

    def __str__(self):
        return f'{self.__value} {self.__units}'

    @property
    def units(self):
        return self.__units

    @property
    def value(self) -> float:
        return self.__value

    def copy(self) -> 'Value':
    #=========================
        return Value(self.__value, self.__units)

    def set_units(self, units: Units):
    #=================================
        if self.__units is None:
            self.__units = units
        elif units != self.__units:
            raise ValueError(f'Can not reassign Units of a Value ({self.__units} != {units})')

    def set_value(self, value: float):
    #================================
        self.__value = value

#===============================================================================

class Quantity:
    def __init__(self, uri: NamedNode, units: Literal, label: Literal|None=None, variable: Literal|None=None):
        self.__uri = uri
        self.__units = Units.from_ucum(units)
        self.__label = str(label) if label is not None else str(uri)
        self.__variable = str(variable) if variable is not None else self.__label

    @property
    def label(self) -> str|None:
        return self.__label

    @property
    def units(self):
        return self.__units

    @property
    def uri(self):
        return self.__uri

    @property
    def variable(self):
        return self.__variable

#===============================================================================
#===============================================================================
