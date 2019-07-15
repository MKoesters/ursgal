#!/usr/bin/env python3
"""Ursgal main resource for TMT quantification.

Takes an mzML as input, removes coalescence of TMT channels, performs normalization
and computes S2I and S2N 
"""
import os
import pickle
from io import StringIO
from pprint import pprint

import numpy as np
import pandas as pd
import pymzml
import pyqms


# @profile
def fix_crosstalk(channels, impurity_matrix_inversed, params=None):
    """Correct crosstalk between TMT channels.
    
    Args:
        ms2_spectrum (pymzml.spec.Spectrum): input ms2 spectrum
        params (dict, optional): additional params
    """
    corrected_intensities = impurity_matrix_inversed.dot(channels.T)
    return  corrected_intensities


# @profile
def get_intensities(
    masses,
    all_peaks,
    tolerance_unit='da',
    tolerance_ppm=5e-6,
    tolerance_da=0.002
):

    """Extract given masses with intensities from spectrum.
    
    Args:
        masses (list): list of masses/mz values to extract from spectrum
        spectrum (pymzml.spec.Spectrum): Description
    """
    signals = {}
    if tolerance_unit == 'ppm':
        tolerance = (all_peaks * tolerance_ppm)[:,0]
    else:
        tolerance = tolerance_da
    # breakpoint()
    for i, rep_ion in enumerate(masses):
        idx = abs(all_peaks[:,0] - rep_ion) < tolerance
        peak = all_peaks[idx]
        for p in peak:
            signals[p[0]] = p[1]
    return signals


# @profile
def extract_TMT_signals(ms2_spectrum, reporter_ions, tolerance, unit):
    """Summary
    
    Args:
        ms2_spectrum (pymzml.spec.Spectrum): input MS2 spectrum
        TMT_type (str): Type of TMT used (6Plex, 8Plex, 10/11plex, 20plex)
    """
    all_peaks = ms2_spectrum.peaks('centroided')
    if not len(all_peaks) == 0:
        sliced_peaks = all_peaks
        signals = np.array([0 for x in range(len(reporter_ions))])
        if len(sliced_peaks) > 0:
            for i, rep_ion in enumerate(reporter_ions):
                idx = abs(sliced_peaks[:,0] - rep_ion) < tolerance
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


# @profile
def calc_isotope_envelopes(peptides, charges):
    # breakpoint()
    cc = pyqms.chemical_composition.ChemicalComposition()
    cc2pep = {}
    for pep in peptides:
        cc.use(pep)
        cc2pep[cc.hill_notation_unimod()] = pep
    print('Start calculating isotopes')
    lib = pyqms.IsotopologueLibrary(molecules=list(peptides), charges=list(charges), verbose=False)
    print('Finished calculating isotopes')
    # for cc, pep in cc2pep.items():
    #     if cc in lib:
    #         lib[pep] = lib[cc]
    #         del lib[cc]
    return lib


# @profile
def calculate_S2I(ms1_spectrum, peptide, charge, isolation_window_borders, isotope_lib, cc_obj):
    """Calculate Signal to Intensity.
    
    For that, the sum of all isotope intensities in the isolation windows is divided by
    the sum of all intensities in the isolation window.
    
    Args:
        ms1_spectrum (pymzml.spec.Spectrum): input_spectrum
        peptide (str): peptide sequence with mods in unimod style
        charge (int): charge of the precursor peptide
        isolation_window_borders (tuple): upper and lower precurso isolation border
    
    Returns:
        float: signal 2 intensity
    """
    lower_border, upper_border = isolation_window_borders
    mz = ms1_spectrum.peaks('centroided')
    all_ions_in_iso_border = mz[lower_border < mz[:,0]]
    all_ions_in_iso_border = all_ions_in_iso_border[upper_border > all_ions_in_iso_border[:,0]]
    cc_obj.use(peptide)
    formula = cc_obj.hill_notation_unimod()
    isotope_mz = isotope_lib[formula]['env'][(('N', '0.000'),)][charge]['mz']
    isotope_int = get_intensities(
        isotope_mz,
        all_ions_in_iso_border,
        tolerance_unit='ppm'
    )
    isotope_int_sum = sum(isotope_int.values())
    all_ions_int_sum = sum(all_ions_in_iso_border[:,1]) + 0.000001
    S2I = isotope_int_sum / all_ions_int_sum
    return S2I


# @profile
def calculate_P2T(ms1_spectrum, peptide, charge, isolation_window_borders, isotope_lib, cc_obj):
    """Calculate Precursor to Threshold.
    
    Divide the sum of all isotope intensities in isolation window by the noise level.
    
    Args:
        ms1_spectrum (pymzml.spec.Spectrum): input_spectrum
        peptide (str): peptide sequence with mods in unimod style
        charge (int): charge of the precursor peptide
        isolation_window_borders (tuple): upper and lower precurso isolation border
    """
    lower_border, upper_border = isolation_window_borders
    mz = ms1_spectrum.peaks('centroided')
    all_ions_in_iso_border = mz[lower_border < mz[:,0]]
    all_ions_in_iso_border = all_ions_in_iso_border[upper_border > all_ions_in_iso_border[:,0]]
    cc_obj.use(peptide)
    formula = cc_obj.hill_notation_unimod()
    isotope_mz = isotope_lib[formula]['env'][(('N', '0.000'),)][charge]['mz']
    isotope_int = get_intensities(
        isotope_mz,
        all_ions_in_iso_border,
        tolerance_unit='ppm',
    )

    isotope_int_sum = sum(isotope_int.values())
    noise_level = ms1_spectrum.estimated_noise_level(mode='median')
    P2T = isotope_int_sum / noise_level
    return P2T


# @profile
def correct_S2I(ms1_spectrum, params=None):
    """Correct 
    
    Args:
        spectrum (pymzml.spec.Spectrum): input spectrum
        params (dict, optional): additional params
    """
    pass


# @profile
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


# @profile
def prepare_idents_df(csv_path):
    """Prepare idents df by adding a column with peptide sequence and mods in unimod style.
    
    Args:
        csv_path (str): path to ursgal unified idents csv.
    
    Returns:
        pandas.DataFrame: pandas dataframe with `full_seq` column
    """
    validated_idents = pd.read_csv(csv_path)
    validated_idents.sort_values('Spectrum ID', inplace=True)
    validated_idents['Modifications'] = validated_idents['Modifications'].fillna('')
    validated_idents['full_seq'] = validated_idents['Sequence'].astype(str) + '#' + validated_idents['Modifications'].astype(str)
    validated_idents['full_seq'] = validated_idents['full_seq'].str.rstrip('#$')
    validated_idents['full_seq'] = validated_idents['full_seq'].dropna()
    return validated_idents


# @profile
def main(mzml_file, output_file, param_dict):
    """
    Takes an mzML as input, removes coalescence of TMT channels, performs normalization
    and computes S2I and S2N. 
    
    Args:
        mzml_file (str): path to input mzml
        param_dict (dict): dict with all node related parameters
    """

    validated_idents = prepare_idents_df(param_dict['quantification_evidences'])

    all_charges = validated_idents['Charge'].unique()
    isotope_lib = calc_isotope_envelopes(validated_idents['full_seq'], all_charges)
    rt_pickle_path = os.path.join(
        os.path.dirname(mzml_file),
        param_dict['rt_pickle_name']
    )
    lookup = pickle.load(
        open(
            os.path.join(
                os.path.dirname(mzml_file),
                param_dict['rt_pickle_name']
            ), 'rb'
        )
    )

    ####################
    impurity_data = StringIO(param_dict['impurity_matrix'])
    impurity_data = pd.read_csv(impurity_data)

    impurity_matrix = impurity_data.drop(columns='Mass-Tag').values
    impurity_matrix_inversed = np.linalg.inv(impurity_matrix.T)

    try:
        basename = os.path.basename(mzml_file).split('.', 1)[0]
        lookup = lookup[basename]
    except KeyError:
        raise KeyError(
            'mzML file {0} not found in RT lookup. Please update RT Lookup with the input mzML and rerun Node'.format(os.path.basename(mzml_file))
        )
    reader = pymzml.run.Reader(mzml_file)
    previous_ms2_idents = []

    tmt_data = {}
    S2I_data = {}
    P2T_data = {}

    cc_object = pyqms.chemical_composition.ChemicalComposition()

    for i, spec in enumerate(reader):
        if i % 500 == 0:
            print(f'Process spec {i}', end='\r')
        # if i == 3000:
        #     break
        if spec.ms_level == 1:
            ms1_spec = spec
            if len(previous_ms2_idents) > 0:
                for ident_info in previous_ms2_idents:
                    S2I = calculate_S2I(
                        ms1_spec,
                        ident_info['seq'],
                        ident_info['charge'],
                        ident_info['isolation_window_borders'],
                        isotope_lib,
                        cc_object
                    )
                    S2I_data[ident_info['ID']]['S2I_after'] = S2I
                    S2I_data[ident_info['ID']]['RT_after'] = ms1_spec.scan_time_in_minutes()
                    P2T = calculate_P2T(
                        ms1_spec,
                        ident_info['seq'],
                        ident_info['charge'],
                        ident_info['isolation_window_borders'],
                        isotope_lib,
                        cc_object
                    )
                    P2T_data[ident_info['ID']]['P2T_after'] = P2T
                    P2T_data[ident_info['ID']]['RT_after'] = ms1_spec.scan_time_in_minutes()
                previous_ms2_idents = []
        elif spec.ms_level == 2:
            ms2_spec = spec
            raw_channels = extract_TMT_signals(
                ms2_spec,
                param_dict['reporter_ion_mzs'],
                param_dict['reporter_ion_tolerance'],
                param_dict['reporter_ion_tolerance_unit']
            )
            # raw_channels = get_intensities(
            #     param_dict['reporter_ion_mzs'],
            #     ms2_spec.peaks('centroided'),
            #     tolerance_unit='da',
            #     tolerance_ppm=0.002
            # )
            # breakpoint()
            fixed_channels = fix_crosstalk(raw_channels, impurity_matrix_inversed)
            normalized_channels = fixed_channels / fixed_channels.sum()
            tmt_data[spec.ID] = fixed_channels

            ident_line = validated_idents[
                (validated_idents['Spectrum ID'] == spec.ID) &
                (validated_idents['Rank'] == 1)
            ]
            if len(ident_line) == 0:
                continue
            current_ident = ident_line[['Sequence', 'Modifications']]
            charge = ident_line['Charge'].iloc[0]
            try:
                full_seq = '#'.join(current_ident.iloc[0])
            except TypeError:
                full_seq = current_ident['Sequence']

            iso_mz = spec['isolation window target m/z']
            lower_border =  iso_mz - spec['isolation window lower offset']
            upper_border = iso_mz + spec['isolation window upper offset']
            ident_info = {
                'seq': full_seq,
                'charge': charge,
                'isolation_window_borders': (lower_border, upper_border),
                'RT': ms2_spec.scan_time_in_minutes(),
                'ID': ms2_spec.ID,
            }
            previous_ms2_idents.append(ident_info)
            S2I = calculate_S2I(
                ms1_spec,
                ident_info['seq'],
                ident_info['charge'],
                ident_info['isolation_window_borders'],
                isotope_lib,
                cc_object
            )
            if ms2_spec.ID not in S2I_data:
                S2I_data[ms2_spec.ID] = {
                    'S2I_before': S2I,
                    'S2I_after': None,
                    'RT_before': ms1_spec.scan_time_in_minutes(),
                    'RT_after': None,
                    'RT_MS2': ms2_spec.scan_time_in_minutes(),
                }
            P2T = calculate_P2T(
                ms1_spec,
                ident_info['seq'],
                ident_info['charge'],
                ident_info['isolation_window_borders'],
                isotope_lib,
                cc_object
            )
            if ms2_spec.ID not in P2T_data:
                P2T_data[ms2_spec.ID] = {
                    'P2T_before': P2T,
                    'P2T_after': None,
                    'RT_before': ms1_spec.scan_time_in_minutes(),
                    'RT_after': None,
                    'RT_MS2': ms2_spec.scan_time_in_minutes(),
                }
    interpol_S2I_data = {}
    interpol_P2T_data = {}
    for key, value in S2I_data.items():
        interpol_S2I = interpolate_qc(
            value['S2I_before'], value['S2I_after'],
            value['RT_before'], value['RT_after'],
            value['RT_MS2']
        )
        value.update({'S2I_interpol': interpol_S2I})
        interpol_S2I_data[key] = value
    for key, value in P2T_data.items():
        interpol_P2T = interpolate_qc(
            value['P2T_before'], value['P2T_after'],
            value['RT_before'], value['RT_after'],
            value['RT_MS2']
        )
        value.update({'P2T_interpol': interpol_P2T})
        interpol_P2T_data[key] = value
    # TODO Write output file
    with open(output_file, 'wt') as fout:
        fout.write('test')
