"""Tests for TMT quant node
"""
import os

import pytest

from TMT_quant_0_0_1 import fix_crosstalk

ISOTOPE_OVERLAP = """
Mass-Tag,126,127L,127H,128L,128H,129L,129H,130L,130H,131L,131H
126,0.92048137463,0.00417403803988,0.0587149294683,0.0,0.00384616691329,0.0,0.0,0.0,0.0,0.0,0.0
127L,0.00619612281462,0.900667547057,0.0221106415412,0.056756557396,0.0,0.00336298002748,0.0,0.0,0.0,0.0,0.0
127H,0.00849174724589,0.0,0.919886648804,0.0042812036856,0.0514567152392,0.0,0.00353336477389,0.0,0.0,0.0,0.0
128L,0.0,0.00819531631876,0.00557034275631,0.904200326705,0.0238942113867,0.0479658854406,0.0,0.00274980930335,0.0,0.0,0.0
128H,0.0,0.0,0.0148826706851,0.0,0.920354911702,0.00407441326769,0.0416243485347,0.0,0.003048713424,0.0,0.0
129L,0.0,0.0,0.0,0.0148088852717,0.00563317274792,0.924006221484,0.0,0.0407297385733,0.0,0.00280137752717,0.0
129H,0.0,0.0,0.0,0.00209994968005,0.0219584782296,0.0,0.926572835987,0.00343018728874,0.0352243843091,0.0,0.00220633892969
130L,0.0,0.0,0.0,0.00184737739711,0.0,0.0226891566799,0.00637595993929,0.908083634915,0.0201297814059,0.0327456482015,0.0
130H,0.0,0.0,0.0,0.0,0.00226685881517,0.0,0.0270021540993,0.0,0.929046590085,0.0037428792119,0.0277780298204
131L,0.0,0.0,0.0,0.0,0.0,0.00214048204202,0.0,0.0278155481012,0.00600576057087,0.943689752776,0.0114510264682
131H,0.0,0.0,0.0,0.0,0.0,0.0,0.00233662045211,0.0,0.0256287420956,0.0,0.958053930181
"""


def test_fix_crosstalk():
    breakpoint()
    pass


def test_get_intensities():
    assert 1 == 2


def test_extract_TMT_signals():
    pass


def test_get_isotope_envelope():
    pass


def test_calc_isotope_envelopes():
    pass


def test_calculate_S2I():
    pass


def test_correct_S2I():
    pass


def test_calculate_P2T():
    pass