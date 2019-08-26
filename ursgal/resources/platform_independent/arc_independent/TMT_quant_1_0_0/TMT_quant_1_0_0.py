#!/usr/bin/env python3
"""Ursgal main resource for TMT quantification.

Takes an mzML as input, removes crosstalk of TMT channels, performs normalization
and computes S2I and P2T
"""
import os
import pickle
import time
from collections import OrderedDict as ODict
from io import StringIO

from scipy import linalg
import numpy as np
import pandas as pd
import pymzml
import pyqms

PROTON = 1.007276466583


def fix_crosstalk_linalg_solve(channels, matrix):
    corrected_intensities = linalg.solve(matrix, channels)
    return corrected_intensities


def fix_crosstalk_dot(channels, impurity_matrix_inversed, params=None):
    """Correct crosstalk between TMT channels.

    Args:
        ms2_spectrum (pymzml.spec.Spectrum): input ms2 spectrum
        params (dict, optional): additional params

    """
    corrected_intensities = impurity_matrix_inversed.dot(channels)
    return corrected_intensities


def get_intensities(
    masses,
    all_peaks,
    tolerance_unit="da",
    tolerance_ppm=5e-6,
    tolerance_da=0.002,
    offset=0,
):
    """Extract given masses with intensities from spectrum.

    Args:
        masses (list): list of masses/mz values to extract from spectrum
        all_peaks (np.ndarray): mz, i array
        tolerance_unit (str, optional): tolerane unit
        tolerance_ppm (float, optional): use tolerance in ppm (relative)
        tolerance_da (float, optional): use tolerance in da (absolute)

    """
    signals = {}
    if tolerance_unit.lower() == "ppm":
        tolerance = (all_peaks * tolerance_ppm)[:, 0]
    elif tolerance_unit.lower() == "da":
        tolerance = tolerance_da
    else:
        raise Exception("Dont know unit: {tolerance_unit}. Please use either ppm or da")
    for i, rep_ion in enumerate(masses):
        idx = abs(all_peaks[:, 0] - rep_ion) < tolerance
        peak = all_peaks[idx]
        for p in peak:
            signals[p[0]] = p[1]
    return signals


def extract_reporter_signals(ms2_spectrum, reporter_ions, tolerance, unit, offset):
    """Extract signals of given reporter ions.

    Args:
        ms2_spectrum (pymzml.spec.Spectrum): input MS2 spectrum
        reporter_ions (dict): dict mapping reporter ion name to mass
        tolerance (TYPE): tolerance in dalton

    Returns:
        list: sorted list of reporter intensities

    """
    all_peaks = ms2_spectrum.peaks("centroided")
    all_peaks[:, 0] += all_peaks[:, 0] * offset
    if not len(all_peaks) == 0:
        sliced_peaks = all_peaks
        signals = np.array([0 for x in range(len(reporter_ions))])
        if len(sliced_peaks) > 0:
            for i, (triv_name, rep_ion_mz) in enumerate(sorted(reporter_ions.items())):
                if unit.lower() == "ppm":
                    tolerance = tolerance * rep_ion_mz * 5e-6
                idx = abs(sliced_peaks[:, 0] - rep_ion_mz) <= tolerance
                test = sliced_peaks[idx]
                if len(test) == 1:
                    signals[i] = test[0][1]
                else:
                    signals[i] = 0
        else:
            signals = []
    else:
        signals = []
    return signals


def calculate_S2I(ms1_spectrum, peptide, charge, isolation_window_borders, offset):
    """Calculate Signal to Intensity.

    For that, the sum of all isotope intensities in the isolation windows is divided by
    the sum of all intensities in the isolation window.

    Args:
        ms1_spectrum (pymzml.spec.Spectrum): input_spectrum
        peptide (str): peptide sequence with mods in unimod style
        charge (int): charge of the precursor peptide
        isolation_window_borders (tuple): upper and lower precurso isolation border

    """
    lower_border, upper_border = isolation_window_borders
    mz = ms1_spectrum.peaks("centroided")
    mz[:, 0] += mz[:, 0] * offset
    all_ions_in_iso_border = mz[lower_border < mz[:, 0]]
    all_ions_in_iso_border = all_ions_in_iso_border[
        upper_border > all_ions_in_iso_border[:, 0]
    ]
    isotope_mz = [peptide, peptide + (PROTON / charge)]
    isotope_int = get_intensities(
        isotope_mz, all_ions_in_iso_border, tolerance_unit="ppm", offset=offset
    )

    isotope_int_sum = sum(isotope_int.values())
    all_ions_int_sum = sum(all_ions_in_iso_border[:, 1]) + 0.000001
    S2I = isotope_int_sum / all_ions_int_sum
    return S2I


def calculate_P2T(ms1_spectrum, peptide, charge, isolation_window_borders, offset):
    """Calculate Precursor to Threshold.

    Divide the sum of all isotope intensities in isolation window by the noise level.

    Args:
        ms1_spectrum (pymzml.spec.Spectrum): input_spectrum
        peptide (str): peptide sequence with mods in unimod style
        charge (int): charge of the precursor peptide
        isolation_window_borders (tuple): upper and lower precurso isolation border

    """
    lower_border, upper_border = isolation_window_borders
    mz = ms1_spectrum.peaks("centroided")
    mz[:, 0] += mz[:, 0] * offset
    all_ions_in_iso_border = mz[lower_border < mz[:, 0]]
    all_ions_in_iso_border = all_ions_in_iso_border[
        upper_border > all_ions_in_iso_border[:, 0]
    ]
    isotope_mz = [peptide, peptide + (PROTON / charge)]
    isotope_int = get_intensities(
        isotope_mz, all_ions_in_iso_border, tolerance_unit="ppm", offset=offset
    )

    isotope_int_sum = sum(isotope_int.values())
    noise_level = ms1_spectrum.estimated_noise_level(mode="median")
    P2T = isotope_int_sum / noise_level
    return P2T


def correct_S2I(raw_channels, S2I):
    """Correct Signal 2 Intensity.

    To be implemented

    Args:
        raw_channels (np.ndarray): array with measured channel intensities
        S2I (float): Signal 2 Intensity

    Returns:
        np.ndarray: array with S2I corrected intensities

    """
    normalized_channels = raw_channels / (raw_channels.sum() + 1e-5)
    interfering_signal = sum(raw_channels * (1 - S2I))
    interfering_signal_per_reporter = normalized_channels * interfering_signal
    no_inference_channels = raw_channels - interfering_signal_per_reporter
    corrected_channels = no_inference_channels / (raw_channels.sum() + 1e-5)
    return corrected_channels


def interpolate_qc(qc_e, qc_l, RT_e, RT_l, RT_ms2):
    """Calculate time weighted linear combination to extrapolate quality control values.

    Args:
        qc_e (float): quality score of early scan
        qc_l (float): quality score of late scan
        RT_e (float): Retention time of early scan
        RT_l (float): Retention time of late scan
        RT_ms2 (float): Retention time of MS2 scan

    Returns:
        float: Interpolated QC value

    """
    return (RT_ms2 - RT_e) * ((qc_l - qc_e) / (RT_l - RT_e)) + qc_e


def main(mzml_file, output_file, param_dict):
    """
    Perform quantification based on mzML input file.

    Take an mzML as input, remove coalescence of TMT channels, perform normalization
        and compute S2I and S2N.

    Args:
        mzml_file (str): path to input mzml
        output_file (str): path to output csv
        param_dict (dict): dict with all node related parameters

    """
    param_dict["reporter_ion_mzs"] = ODict(
        list(param_dict["reporter_ion_mzs"].values())[0]
    )
    trivial_names = list(param_dict["reporter_ion_mzs"].keys())

    reporter_ion_tolerance = list(param_dict["reporter_ion_tolerance"].values())[0]
    reporter_ion_tolerance_unit = list(
        param_dict["reporter_ion_tolerance_unit"].values()
    )[0]

    ppm_offset = list(param_dict["machine_offset_in_ppm"].values())[0]
    offset = ppm_offset * 1e-6

    impurity_data = StringIO(list(param_dict["impurity_matrix"].values())[0])
    impurity_data = pd.read_csv(impurity_data)
    impurity_matrix = impurity_data.drop(columns="Mass-Tag").values
    impurity_matrix_transposed = impurity_matrix.T
    impurity_matrix_inversed = np.linalg.inv(impurity_matrix_transposed)

    reader = pymzml.run.Reader(mzml_file)
    previous_ms2_idents = []

    all_lines = {}
    all_data = {}
    tmt_data = {}
    S2I_data = {}
    P2T_data = {}

    for i, spec in enumerate(reader):
        if i % 500 == 0:
            print(
                "Process spec {i} with MS level {level}".format(
                    i=i, level=spec.ms_level
                ),
                end="\r",
            )
        if spec.ms_level == 1:
            ms1_spec = spec
            if len(previous_ms2_idents) > 0:
                for ident_info in previous_ms2_idents:
                    s2i_spec_data = S2I_data[ident_info["ID"]]
                    p2t_spec_data = P2T_data[ident_info["ID"]]
                    S2I = calculate_S2I(
                        ms1_spec,
                        ident_info["mz"],
                        ident_info["charge"],
                        ident_info["isolation_window_borders"],
                        offset,
                    )
                    s2i_spec_data["S2I_after"] = S2I
                    s2i_spec_data["RT_after"] = ms1_spec.scan_time_in_minutes()
                    P2T = calculate_P2T(
                        ms1_spec,
                        ident_info["mz"],
                        ident_info["charge"],
                        ident_info["isolation_window_borders"],
                        offset,
                    )
                    p2t_spec_data["P2T_after"] = P2T
                    p2t_spec_data["RT_after"] = ms1_spec.scan_time_in_minutes()
                    interpol_S2I = interpolate_qc(
                        s2i_spec_data["S2I_before"],
                        s2i_spec_data["S2I_after"],
                        s2i_spec_data["RT_before"],
                        s2i_spec_data["RT_after"],
                        s2i_spec_data["RT_MS2"],
                    )
                    interpol_P2T = interpolate_qc(
                        p2t_spec_data["P2T_before"],
                        p2t_spec_data["P2T_after"],
                        p2t_spec_data["RT_before"],
                        p2t_spec_data["RT_after"],
                        p2t_spec_data["RT_MS2"],
                    )
                    s2i_corrected_channels = correct_S2I(
                        tmt_data[ident_info["ID"]]["raw"], interpol_S2I
                    )
                    # Implement normalization over peptides
                    for i, channel in enumerate(tmt_data[ident_info["ID"]]["raw"]):
                        # breakpoint()
                        nd = {**all_lines[ident_info["ID"]]}
                        nd["Quant Ion"] = list(trivial_names)[i]
                        nd["rel Intensity"] = tmt_data[ident_info["ID"]]["normalized"][
                            i
                        ]
                        nd["abs Intensity"] = tmt_data[ident_info["ID"]]["raw"][i]
                        nd["corr Intensity"] = s2i_corrected_channels[i]
                        nd["P2T"] = interpol_P2T
                        nd["S2I"] = interpol_S2I
                        all_data[ident_info["ID"]].append(nd)
                previous_ms2_idents = []
        elif spec.ms_level == 2:
            ms2_spec = spec
            raw_channels = extract_reporter_signals(
                ms2_spec,
                param_dict["reporter_ion_mzs"],
                reporter_ion_tolerance,
                reporter_ion_tolerance_unit,
                offset,
            )
            # fixed_channels = fix_crosstalk_dot(raw_channels, impurity_matrix_transposed)
            fixed_channels = fix_crosstalk_linalg_solve(
                raw_channels, impurity_matrix_inversed
            )
            # if all(fixed_channels == 0):
            #     continue
            # assert np.allclose(fixed_channels, fixed_channels2), f'{fixed_channels}\n{fixed_channels2}'
            normalized_channels = fixed_channels / fixed_channels.sum()
            if spec.ID not in tmt_data:
                tmt_data[spec.ID] = {"raw": None, "fixed": None}
            tmt_data[spec.ID]["raw"] = fixed_channels
            tmt_data[spec.ID]["normalized"] = normalized_channels
            if spec.ID not in all_data:
                all_data[spec.ID] = []

            iso_mz = spec["isolation window target m/z"]
            lower_border = iso_mz - spec["isolation window lower offset"]
            upper_border = iso_mz + spec["isolation window upper offset"]
            charge = spec["charge state"]
            ident_info = {
                "mz": iso_mz,
                "charge": charge,
                "isolation_window_borders": (lower_border, upper_border),
                "RT": ms2_spec.scan_time_in_minutes(),
                "ID": ms2_spec.ID,
            }
            previous_ms2_idents.append(ident_info)
            S2I = calculate_S2I(
                ms1_spec,
                ident_info["mz"],
                ident_info["charge"],
                ident_info["isolation_window_borders"],
                offset,
            )
            if ms2_spec.ID not in S2I_data:
                S2I_data[ms2_spec.ID] = {
                    "S2I_before": S2I,
                    "S2I_after": None,
                    "RT_before": ms1_spec.scan_time_in_minutes(),
                    "RT_after": None,
                    "RT_MS2": ms2_spec.scan_time_in_minutes(),
                }
            P2T = calculate_P2T(
                ms1_spec,
                ident_info["mz"],
                ident_info["charge"],
                ident_info["isolation_window_borders"],
                offset,
            )
            if ms2_spec.ID not in P2T_data:
                P2T_data[ms2_spec.ID] = {
                    "P2T_before": P2T,
                    "P2T_after": None,
                    "RT_before": ms1_spec.scan_time_in_minutes(),
                    "RT_after": None,
                    "RT_MS2": ms2_spec.scan_time_in_minutes(),
                }
            line = ODict(
                {
                    "Spectrum ID": ms2_spec.ID,
                    "Charge": charge,
                    "RT": ms2_spec.scan_time_in_minutes(),
                    "Quant Ion": list(trivial_names),
                    "rel Intensity": normalized_channels,
                    "abs Intensity": fixed_channels,
                    "S2I": None,
                    "P2T": None,
                }
            )
            all_lines[ms2_spec.ID] = line
    flat = [item for list in all_data.values() for item in list]
    final_df = pd.DataFrame(flat)
    float_cols = {"rel Intensity": 5, "abs Intensity": 5, "RT": 3, "S2I": 5, "P2T": 5}
    final_df = final_df.round(float_cols)
    # sometimes there are -0 and 0, making it hard to compare files using hashes
    final_df["rel Intensity"] += 0.0
    final_df["abs Intensity"] += 0.0
    final_df.to_csv(output_file, index=False)
