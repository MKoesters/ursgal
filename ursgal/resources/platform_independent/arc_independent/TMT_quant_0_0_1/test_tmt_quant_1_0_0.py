"""Tests for TMT quant node
"""
import os
from io import StringIO

import numpy as np
import pandas as pd
import pytest
from pyqms.chemical_composition import ChemicalComposition

from TMT_quant_0_0_1 import (
    calc_isotope_envelopes,
    calculate_P2T,
    calculate_S2I,
    fix_crosstalk,
    get_intensities,
    interpolate_qc,
)

PROTON = 1.007_276_466_583

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


class Spectrum:
    def __init__(self, peaks=None, noise_level=1):
        self._peaks = peaks
        self._noise_level = noise_level

    def peaks(self, type=None):
        return self._peaks

    def estimated_noise_level(self, mode="median"):
        return self._noise_level


def test_fix_crosstalk():
    """Test crosstalk fix.

    Use precomputed impurities and assert fix sets all values to 100.0
    """
    TMT_channels = np.array(
        [
            93.51692447,
            91.30369014,
            102.11652333,
            98.39943001,
            102.9410515,
            100.42391389,
            100.74452838,
            98.28089182,
            101.90839719,
            98.29796577,
            99.94893254,
        ]
    )
    matrix = StringIO(ISOTOPE_OVERLAP)
    matrix = pd.read_csv(matrix)
    impurity_matrix = matrix.drop(columns="Mass-Tag").values
    impurity_matrix_inversed = np.linalg.inv(impurity_matrix.T)
    fixed = fix_crosstalk(TMT_channels, impurity_matrix_inversed)
    assert np.allclose(fixed, 100)


@pytest.mark.parametrize(
    "mz,peaks,unit,tol,expected_result",
    [
        (
            [126],
            np.array([(126.002, 100), (126.003, 100), (125.998, 100)]),
            "da",
            0.003,
            2,
        ),
        (
            [126],
            np.array([(126.002, 100), (126.003, 100), (125.998, 100)]),
            "da",
            0.002,
            2,
        ),
        (
            [126],
            np.array([(126.002, 100), (126.003, 100), (125.998, 100)]),
            "da",
            0.0031,
            3,
        ),
        (
            [126],
            np.array([(126.002, 100), (126.003, 100), (125.998, 100)]),
            "da",
            0.0,
            0,
        ),
        (
            [126],
            np.array([(126.002, 100), (126.003, 100), (125.998, 100)]),
            "ppm",
            20e-5,
            3,
        ),
        (
            [126],
            np.array([(126.002, 100), (126.003, 100), (125.998, 100)]),
            "ppm",
            15.87396e-6,
            2,
        ),
        (
            [126],
            np.array([(126.002, 100), (126.003, 100), (125.998, 100)]),
            "ppm",
            1e-5,
            0,
        ),
        (
            [126],
            np.array([(126.002, 100), (126.003, 100), (125.998, 100)]),
            "ppm",
            0,
            0,
        ),
    ],
    ids=[
        "tol=0.0030 unit=da expected_peaks=2",
        "tol=0.0020 unit=da expected_peaks=2",
        "tol=0.0031 unit=da expected_peaks=3",
        "tol=0.0000 unit=da expected_peaks=0",
        "tol=20e-5  unit=ppm expected_peaks=3",
        "tol=15e-6  unit=ppm expected_peaks=2",
        "tol=1e-5   unit=ppm expected_peaks=0",
        "tol=0      unit=ppm expected_peaks=0",
    ],
)
def test_get_intensities(mz, peaks, unit, tol, expected_result):
    """Test extracting intensities for given masses from peak list.
    
    Args:
        mz (list): mz to extract
        peaks (np.array): peak list
        unit (str): tolerance unit
        tol (float): tolerance
        expected_result (int): number of extracted peaks
    """
    if unit == "da":
        extracted_peaks = get_intensities(
            mz, peaks, tolerance_unit=unit, tolerance_da=tol
        )
    else:
        extracted_peaks = get_intensities(
            mz, peaks, tolerance_unit=unit, tolerance_ppm=tol
        )

    test_result = len(extracted_peaks)
    assert expected_result == test_result


@pytest.mark.parametrize(
    "peptides,charges,expected_results",
    [
        (["ELVISLIVES"], [1], (1, 1100.633 + PROTON)),
        # test some more cases
    ],
    ids=["Peptide=ELVISLIVES Charge=1"],
)
def test_calc_isotope_envelopes(peptides, charges, expected_results):
    """Test calculation of isotope envelopes.
    
    Args:
        peptides (list): list of peptide sequences
        charges (list): list of charge states
        expected_results (tuple): tuple with monoisotopic mass
    """
    lib = calc_isotope_envelopes(peptides, charges)
    assert len(lib.keys()) == expected_results[0]
    peptide_cc = ChemicalComposition(peptides[0]).hill_notation_unimod()
    assert [peptide_cc] == list(lib.keys())
    arr = np.array(lib[peptide_cc]["env"][(("N", "0.000"),)][charges[0]]["mz"])
    assert (abs(arr - expected_results[1]) < expected_results[1] * 5e-6).any()


@pytest.mark.parametrize(
    "ms1_spec,peptide,charge,window,result",
    [
        (
            Spectrum(
                peaks=np.array(
                    [
                        (1101.6401678803702, 100),
                        (1102.6431931587235, 400),
                        (1101.8401678803702, 500),
                    ]
                )
            ),
            "ELVISLIVES",
            1,
            (1100.633 + PROTON - 2, 1100.633 + PROTON + 2),
            0.5,
        ),
        (
            Spectrum(
                peaks=np.array([(1101.6401678803702, 100), (1102.6431931587235, 400)])
            ),
            "ELVISLIVES",
            1,
            (1100.633 + PROTON - 2, 1100.633 + PROTON + 2),
            1.0,
        ),
        (
            Spectrum(peaks=np.array([(1101.8401678803702, 500)])),
            "ELVISLIVES",
            1,
            (1100.633 + PROTON - 2, 1100.633 + PROTON + 2),
            0.0,
        ),
    ],
    ids=[
        "Isotopes=2 prec Int=500 total Int=100",
        "Isotopes=2 prec Int=500 total Int=500",
        "Isotopes=2 prec Int=0 total Int=500",
    ],
)
def test_calculate_S2I(ms1_spec, peptide, charge, window, result):
    """Test signal 2 intensity calculation
    
    Args:
        ms1_spec (Spectrum): spectrum class with peaks() and estimated_noise_level methods
        peptide (str): peptide sequence
        charge (int): Description
        window (tuple): lower and upper isolation window border
        result (float): expeted S2I
    """
    cc_obj = ChemicalComposition()
    lib = lib = calc_isotope_envelopes([peptide], [charge])
    s2i = calculate_S2I(ms1_spec, peptide, charge, window, lib, cc_obj)
    assert np.isclose(s2i, result)


@pytest.mark.parametrize(
    "ms1_spec,peptide,charge,window,result",
    [
        (
            Spectrum(
                peaks=np.array(
                    [
                        (1101.6401678803702, 100),
                        (1102.6431931587235, 400),
                        (1101.8401678803702, 500),
                    ]
                ),
                noise_level=1,
            ),
            "ELVISLIVES",
            1,
            (1100.633 + PROTON - 2, 1100.633 + PROTON + 2),
            500,
        ),
        (
            Spectrum(
                peaks=np.array(
                    [
                        (1101.6401678803702, 100),
                        (1102.6431931587235, 400),
                        (1101.8401678803702, 500),
                    ]
                ),
                noise_level=2,
            ),
            "ELVISLIVES",
            1,
            (1100.633 + PROTON - 2, 1100.633 + PROTON + 2),
            250,
        ),
        (
            Spectrum(
                peaks=np.array(
                    [
                        (1101.6401678803702, 100),
                        (1102.6431931587235, 400),
                        (1101.8401678803702, 500),
                    ]
                ),
                noise_level=0.5,
            ),
            "ELVISLIVES",
            1,
            (1100.633 + PROTON - 2, 1100.633 + PROTON + 2),
            1000,
        ),
        (
            Spectrum(peaks=np.array([(1101.8401678803702, 500)]), noise_level=1),
            "ELVISLIVES",
            1,
            (1100.633 + PROTON - 2, 1100.633 + PROTON + 2),
            0.0,
        ),
    ],
    ids=[
        "2 Isotopes noise level=1",
        "2 Isotopes noise level=2",
        "2 Isotopes noise level=0.5",
        "0 Isotopes noise level=1",
    ],
)
def test_calculate_P2T(ms1_spec, peptide, charge, window, result):
    cc_obj = ChemicalComposition()
    lib = lib = calc_isotope_envelopes([peptide], [charge])
    p2t = calculate_P2T(ms1_spec, peptide, charge, window, lib, cc_obj)
    assert np.isclose(p2t, result)


@pytest.mark.parametrize(
    "input_data,output",
    [
        ((0.5, 0.5, 1, 3, 2), 0.5),  # same S2I always leeds to same output
        ((100, 100, 1, 5, 2), 100),
        ((0.25, 0.75, 1, 3, 2), 0.5),  # different S2I, ms2 same distance to e and l
        ((0.75, 0.25, 1, 3, 2), 0.5),  # and the other way around :)
        ((-0.75, -0.25, 1, 3, 2), -0.5),  # try negative values
        ((-0.25, -0.75, 1, 3, 2), -0.5),  # different S2I, ms2 same distance to e and l
        ((0.25, 0.75, 1, 1e10, 2), 0.25),  # ms2 way closer to early scan
        ((0.25, 0.75, 1, 1e10, 1e10 - 1), 0.75),  # ms2 way closer to early scan
    ],
    ids=[
        "Same S2I",
        "Same S2I high value",
        "Different S2I equidistance",
        "Different S2I equidistance reverse",
        "Negative S2I",
        "Negative S2I reverse",
        "MS2 close to early scan",
        "MS2 close to late scan",
    ],
)
def test_interpolate_qc(input_data, output):
    qc_e, qc_l, RT_e, RT_l, RT_ms2 = input_data
    interpolated_qc = interpolate_qc(qc_e, qc_l, RT_e, RT_l, RT_ms2)
    # assert interpolated_qc == output
    assert np.isclose(interpolated_qc, output)


if __name__ == "__main__":
    pytest.main()
