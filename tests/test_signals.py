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

"""Signals (bgf:SignalBond), signal blocks (bgf:SignalBlock) and the elements for vessel, organ and control modules:
the CellML bg2cellml makes for the example models in tests/models, evaluated at their initial values against the
relations worked by hand, and the issues a wrongly connected signal gives.

    pytest tests/test_signals.py
"""

import asyncio
import math
from pathlib import Path

import pytest

from bg2cellml.bondgraph.framework import get_framework

from .cellml_eval import cellml_from_bgrdf, evaluate

MODELS = Path(__file__).parent / 'models'
RT = 2578.73058

#===============================================================================

def test_a_vessel_with_inertances_a_viscous_wall_and_a_one_way_outlet():
    '''Inlet segment dv/dt = (u - R v)/I, wall u = (q - q_us)/C + R_v dq/dt, outlet inertance dv/dt = u/I into a
    stiffening store u = E (q - q_hat)/(q_bar - q)^2, emptied through v = k max(u - u_0, 0).'''
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'vessel_inertance.ttl'))
    v_in, v_out = 0.1, 0.05
    u_v = (1.2 - 0.5)/0.1 + 0.1*(v_in - v_out)
    u_b = 0.1*(0.2 - 0.05)/(0.6 - 0.2)**2
    assert values['u_q_v'] == pytest.approx(u_v)
    assert rates['v_I_in'] == pytest.approx((13 - u_v - 2*v_in)/0.01)
    assert rates['v_I_out'] == pytest.approx((u_v - u_b)/0.05)
    assert values['u_q_b'] == pytest.approx(u_b)
    assert values['v_R_void'] == 0.0            # u_b is below the opposing pressure, so nothing flows
    assert rates['q_b'] == pytest.approx(v_out) and rates['q_v'] == pytest.approx(v_in - v_out)

def test_angiotensin_scales_a_resistance_and_a_compliance():
    '''ang_scale = 1 + alpha q/(q + q_h) divides the efferent conductance and the peritubular compliance; the
    glomerular 0-junction's pressure is lagged into an output port.'''
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'angII_vessel.ttl'))
    scale = 1 + 0.5*2e-11/(2e-11 + 2e-11)
    u_gc = (0.6 - 0.1)/0.05
    u_pc = scale*(0.4 - 0.2)/0.1
    assert values['ang_scale'] == pytest.approx(scale)
    assert values['u_gc'] == pytest.approx(u_gc) and values['u_q_pc'] == pytest.approx(u_pc)
    assert values['v_R_eff'] == pytest.approx(0.3*(u_gc - u_pc)/scale)
    assert rates['q_C_gc'] == pytest.approx(0.4*(13 - u_gc) - values['v_R_eff'])
    assert rates['P_gc_f'] == pytest.approx((u_gc - 8)/5) and values['P_gc'] == pytest.approx(8)

def test_the_baroreceptor_lags_pressure_and_maps_it_to_a_firing_rate():
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'baroreceptor.ttl'))
    assert rates['P_baro'] == pytest.approx((13 - 12)/2.076)
    assert values['f_ab'] == pytest.approx((2.52 + 47.78)/2)

def test_renin_secretion_and_the_michaelis_menten_cascade():
    '''Renin's secretion k f_P(P) (q/q_0)^-B, a modulated flow source; Ang I formed at kcat a E x/(x + Km) with
    E renin's amount, a signal, and x = K q of angiotensinogen.'''
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'renin_secretion.ttl'))
    f_P = 1 + 4/(1 + math.exp((12 - 11.5)/0.8))
    assert values['f_P'] == pytest.approx(f_P) and values['fb_ATII'] == pytest.approx(1.0)
    assert values['v_renin'] == pytest.approx(1e-15*f_P)
    x = 3.3e5*3e-6
    v1 = 1.0*1.0*1.3e-12*x/(x + 1.1)
    assert values['v_v1'] == pytest.approx(v1)
    assert rates['q_AGT_AGT'] == pytest.approx(6e-11 - v1)
    assert rates['q_renin_renin'] == pytest.approx(1e-15*f_P - 7.7e-16*(1e12*1.3e-12 - math.exp(-80000/RT)))

def test_aldosterone_secretion_lags_and_scales_a_sodium_pump():
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'aldosterone_transport.ttl'))
    assert values['aldo_target'] == pytest.approx(0.3 + 0.7*1.0*1.0)
    assert rates['aldo_secretion'] == pytest.approx(0.0, abs=1e-15)
    assert values['v_aldo_sec'] == pytest.approx(2.5e-12)
    scale = 1 + 5e-9/(5e-9 + 5e-9)
    x_f, x_r = 100*0.0015, 0.3*0.4
    assert values['v_NKE'] == pytest.approx(scale*1e-4*(x_f - x_r)/(1 + (x_f + x_r)/2))

def test_kedem_katchalsky_water_and_solute_fluxes():
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'kedem_katchalsky.ttl'))
    pi_l, pi_c = 2*RT*0.014/0.1, 2*RT*0.0006/0.05
    c_l, c_c = 0.014/0.1, 0.0006/0.05
    jv = 1e-4*((2.0 - 1.0) - 0.9*(pi_l - pi_c))
    js = 1e-3*(c_l - c_c) + 0.1*(c_l + c_c)/2*jv
    assert values['Pi_lumen'] == pytest.approx(pi_l) and values['c_Na_cell'] == pytest.approx(c_c)
    assert values['v_Jv'] == pytest.approx(jv)
    assert values['v_Js'] == pytest.approx(js)
    assert values['v_Jc'] == pytest.approx(c_c*jv)
    assert rates['q_W_cell'] == pytest.approx(jv)
    assert rates['q_Na_cell_Na_cell'] == pytest.approx(js - c_c*jv)

#===============================================================================

def issues(turtle: str) -> list[str]:
    async def load():
        framework = await get_framework()
        model = framework.make_bondgraph_model('file:///signals.ttl', turtle)
        return model, [issue.reason for issue in model.issues]
    return asyncio.run(load())

PREFIXES = '''@prefix bgf: <https://bg-rdf.org/ontologies/bondgraph-framework#> .
@prefix cdt: <https://w3id.org/cdt/> .
@prefix : <#> .
'''

STORES = PREFIXES + '''
:m a bgf:BondgraphModel ; bgf:hasBondElement :q_1, :R, :q_2 .
:q_1 a bgf:HydraulicComplianceNode ; bgf:hasValue "1 L"^^cdt:ucum ;
    bgf:parameterValue [ bgf:varName "C" ; bgf:hasValue "1 L2/J"^^cdt:ucum ] ,
                       [ bgf:varName "q_us" ; bgf:hasValue "0 L"^^cdt:ucum ] .
:q_2 a bgf:HydraulicComplianceNode ; bgf:hasValue "0.5 L"^^cdt:ucum ;
    bgf:parameterValue [ bgf:varName "C" ; bgf:hasValue "1 L2/J"^^cdt:ucum ] ,
                       [ bgf:varName "q_us" ; bgf:hasValue "0 L"^^cdt:ucum ] .
:R a bgf:HydraulicResistanceNode ; bgf:parameterValue [ bgf:varName "k" ; bgf:hasValue "1 L2.s-1.J-1"^^cdt:ucum ] .
:b1 bgf:hasSource :q_1 ; bgf:hasTarget :R .
:b2 bgf:hasSource :R ; bgf:hasTarget :q_2 .
'''

def test_a_signal_bond_between_elements_is_not_a_power_bond():
    '''A signal from q_1 to q_2 adds no bond to the bond graph's network.'''
    model, found = issues(STORES + ''':s a bgf:SignalBond ; bgf:hasSource :q_1 ; bgf:hasTarget :R ;
        bgf:signalVariable "u" ; bgf:modulates "k" .''')
    assert any('dimensionless' in issue for issue in found), found
    edges = {(a.rsplit('#')[-1], b.rsplit('#')[-1]) for a, b in model.network_graph.edges}
    assert edges == {('q_1', 'R'), ('R', 'q_2')}

def test_a_modulation_by_a_dimensionless_signal_scales_the_parameter():
    model, found = issues(STORES + ''':m bgf:hasSignalBlock :p .
        :p a bgf:PowerLaw ; bgf:parameterValue [ bgf:varName "u_ref" ; bgf:hasValue "0.5 L"^^cdt:ucum ] ,
                                               [ bgf:varName "n" ; bgf:hasValue "2 1"^^cdt:ucum ] .
        :s1 a bgf:SignalBond ; bgf:hasSource :q_1 ; bgf:hasTarget :p .
        :s2 a bgf:SignalBond ; bgf:hasSource :p ; bgf:hasTarget :R ; bgf:modulates "k" .''')
    assert found == []
    values, _ = evaluate(model.make_cellml_model().to_xml())
    assert values['v_R'] == pytest.approx(1.0*(1.0/0.5)**2*(1.0 - 0.5))

def test_an_unconnected_signal_input_is_an_issue():
    _, found = issues(STORES + ''':m bgf:hasSignalBlock :lag .
        :lag a bgf:FirstOrderLag ; bgf:hasValue "1 J/L"^^cdt:ucum ;
            bgf:parameterValue [ bgf:varName "tau" ; bgf:hasValue "1 s"^^cdt:ucum ] .''')
    assert any('Signal input u of block' in issue for issue in found), found

def test_a_signal_of_the_wrong_units_is_an_issue():
    '''A water volume (L) into a Concentration block's solute amount (mol).'''
    _, found = issues(STORES + ''':m bgf:hasSignalBlock :c .
        :c a bgf:Concentration .
        :s1 a bgf:SignalBond ; bgf:hasSource :q_1 ; bgf:hasTarget :c ; bgf:signalInput "q_s" .
        :s2 a bgf:SignalBond ; bgf:hasSource :q_2 ; bgf:hasTarget :c ; bgf:signalInput "q_w" .''')
    assert any('has units' in issue and 'q_s' in issue for issue in found), found

def test_a_parameter_in_other_units_is_converted_to_its_signals():
    '''A lag of a pressure in J/L whose initial value is given in Pa.'''
    model, found = issues(STORES + ''':m bgf:hasSignalBlock :lag .
        :lag a bgf:FirstOrderLag ; bgf:hasSymbol "u_f" ; bgf:hasValue "500 Pa"^^cdt:ucum ;
            bgf:parameterValue [ bgf:varName "tau" ; bgf:hasValue "2 s"^^cdt:ucum ] .
        :s a bgf:SignalBond ; bgf:hasSource :q_1 ; bgf:hasTarget :lag ; bgf:signalVariable "u" .''')
    assert found == []
    values, rates = evaluate(model.make_cellml_model().to_xml())
    assert values['u_f'] == pytest.approx(0.5)
    assert rates['u_f'] == pytest.approx((1.0 - 0.5)/2)

#===============================================================================
#===============================================================================
