#!/usr/bin/env python3.4
import ursgal
from ursgal.chemical_composition import ChemicalComposition
import os
import pymzml
from pprint import pprint


class sirius_3_5(ursgal.UNode):
    """
    MetFrag v2.31 UNode.

    Get CL version at http://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.3.1-CL.jar
    """
    META_INFO = {
        'edit_version'                : 1.00,
        'name'                        : 'Sirius',
        'version'                     : 'v3.5',
        'release_date'                : '2017-01-01',
        'engine_type' : {
            'search_engine' : False,
        },
        'input_extensions'            : ['.mzml', '.mgf', '.txt', '.ms'], # simple mz i pairs in one row as input
        'input_multi_file'            : False,
        'output_extensions'           : ['.csv'],
        'compress_raw_search_results' : False,
        'create_own_folder'           : True,
        'in_development'              : True,
        'include_in_git'              : False,
        'utranslation_style'          : 'sirius_style_1',
        'engine' : {
            'platform_independent' : {
                'arc_independent' : {
                    'exe'            : 'sirius',
                    'url'            : '',
                    'zip_md5'        : '',
                    'additional_exe' : [],
                },
            },
        },
        'citation' : """
                     """
    }

    def __init__(self, *args, **kwargs):
        """Init paren object."""
        super(sirius_3_5, self).__init__(*args, **kwargs)

    def preflight(self):
        """Doc."""
        pprint(self.params['translations'])
        print()
        pprint(self.params)

        ws_dir = os.path.join(
            self.params['output_dir_path'],
            'sirius_workspace'
        )

        comp = ChemicalComposition('+{neutral_precursor_mol_formula}'.format(
            **self.params['translations']
        ))
        comp_mass = comp._mass()
        print(self.params['translations']['neutral_precursor_mol_formula'])
        print(comp_mass)
        print(comp)
        # exit(1)
        self.params['command_list'] = [
            self.exe,
            '-q',
            # '{allowed_elements_key}'.format(**self.params['translations']),
            # '{allowed_elements}'.format(**self.params['translations']),
            # '{neutral_precursor_mol_formula_key}'.format(**self.params['translations']),
            # '{neutral_precursor_mol_formula}'.format(**self.params['translations']),
            '{ionization_key}'.format(**self.params['translations']),
            '{ionization}'.format(**self.params['translations']),
            '{isotope_handling_key}'.format(**self.params['translations']),
            '{isotope_handling}'.format(**self.params['translations']),
            # '{noise_key}'.format(**self.params['translations']),
            # '{noise}'.format(**self.params['translations']),
            '{sirius_candidate_number_key}'.format(**self.params['translations']),
            '{sirius_candidate_number}'.format(**self.params['translations']),
            '{output_file_incl_path_key}'.format(**self.params['translations']),
            '{output_dir_path}'.format(output_dir_path=ws_dir),
            '{precursor_mass_tolerance_minus_key}'.format(**self.params['translations']),
            '{precursor_mass_tolerance_minus}'.format(**self.params['translations']),
            # '{neutral_precursor_mass_key}'.format(**self.params['translations']),
            # str(comp_mass),
            # '{instrument_key}'.format(**self.params['translations']),
            # '{instrument}'.format(**self.params['translations'])
        ]

        # # add flags if set
        if self.params['translations']['auto_charge'] is True:
            self.params['command_list'].append('-Z')
        if self.params['translations']['recalibrate_sirius'] is False:
            self.params['command_list'].append('--recalibrate')

        # add file handling stuff
        if self.params['input_file'].endswith('.mgf'):
            # only add mgf file with flag -2
            self.params['command_list'] += [
                '-2',
                os.path.join(
                    self.params['input_dir_path'],
                    self.params['input_file']
                )
            ]
        elif self.params['input_file'].endswith('.ms'):
            # no need for -1 or -2 flag, just append filename
            self.params['command_list'] += [
                self.params['input_file']
            ]
        # elif self.params['input_file'].endswith('.txt'):
        #     # enable checking if its ms1 or ms2 file,
        #     # then decide which flag to use
        #     self.params['command_list'] += [
        #         '{}'.format(**self.params['translations']),
        #         '{}'.format(**self.params['translations'])
        #     ]
        else:
            raise Exception(
                'Invalid input file (Ursgal should have checked in before)'
            )
        pprint(self.params['command_list'])

    def postflight(self):
        """DOC."""
        # convert output from xls to csv
        # unify output format
        pass

    def _write_ms_file(self, path):
        pass

