#===============================================================================
#
#  CellDL and bondgraph tools
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#===============================================================================

"""
Causality: each algebraic equation computes one variable.

A bond graph's equations are written with a fixed variable on each left-hand side: a resistance's flow from its
potentials, a junction's potential or flow from its bonds'. Where the bond graph's causality differs, e.g. a
resistance between a node a flow source feeds and a compliance (its flow is given, so it must give the node's
pressure), one variable has two equations and another none. ``assign_causality`` then matches each equation to a
variable it computes (keeping its own left-hand side where it can, by augmenting paths from the over-determined
equations to the variables no equation computes) and solves each equation whose variable changed for it.
"""

#===============================================================================

import sympy

from .mathml import Equation

#===============================================================================

class CausalityError(Exception):
    pass

def assign_causality(equations: list[Equation], known: set[str]) -> list[Equation]:
    """The equations, each solved for the variable it computes. ``known``: the names of variables no equation
    computes (parameters, states, time). Equations are returned unchanged when every variable has one equation."""
    algebraic = [n for n, eq in enumerate(equations) if isinstance(eq.lhs, sympy.Symbol)]
    lhs_count: dict[sympy.Symbol, int] = {}
    for n in algebraic:
        lhs_count[equations[n].lhs] = lhs_count.get(equations[n].lhs, 0) + 1
    unknowns_of: dict[int, set[sympy.Symbol]] = {}
    unknown: set[sympy.Symbol] = set()
    for n in algebraic:
        eq = equations[n]
        symbols = {s for s in eq.rhs.free_symbols | {eq.lhs} if isinstance(s, sympy.Symbol) and s.name not in known}
        unknowns_of[n] = symbols
        unknown |= symbols
    for n, eq in enumerate(equations):
        if n not in unknowns_of:
            unknown |= {s for s in eq.rhs.free_symbols if isinstance(s, sympy.Symbol) and s.name not in known}
    undefined = {s for s in unknown if s not in lhs_count}
    if not undefined or all(count == 1 for count in lhs_count.values()):
        # nothing to change (or nothing an equation could be solved for instead)
        return equations

    # the matching: an equation keeps its left-hand side unless another already has it
    match: dict[int, sympy.Symbol] = {}            # equation -> the variable it computes
    owner: dict[sympy.Symbol, int] = {}            # variable -> its equation
    unmatched = []
    for n in algebraic:
        lhs = equations[n].lhs
        if lhs not in owner:
            match[n], owner[lhs] = lhs, n
        else:
            unmatched.append(n)

    def augment(n: int, seen: set[int]) -> bool:
        # give equation n a variable: a free one, or one whose equation can take another
        candidates = sorted(unknowns_of[n], key=lambda s: s.name)
        for symbol in candidates:
            if symbol not in owner:
                match[n], owner[symbol] = symbol, n
                return True
        for symbol in candidates:
            other = owner[symbol]
            if other in seen:
                continue
            seen.add(other)
            if augment(other, seen):
                match[n], owner[symbol] = symbol, n
                return True
        return False

    for n in unmatched:
        if not augment(n, {n}):
            raise CausalityError(f'No variable left for {equations[n].lhs} = {equations[n].rhs}')

    result = list(equations)
    for n, symbol in match.items():
        eq = equations[n]
        if symbol == eq.lhs:
            continue
        solutions = sympy.solve(sympy.Eq(eq.lhs, eq.rhs), symbol, dict=False)
        if len(solutions) != 1:
            raise CausalityError(f'Cannot solve {eq.lhs} = {eq.rhs} for {symbol}')
        result[n] = Equation(symbol, solutions[0], eq.provenance)
    return result

#===============================================================================
#===============================================================================
