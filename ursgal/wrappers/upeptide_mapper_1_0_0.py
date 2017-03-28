#!/usr/bin/env python3.4
import ursgal
import importlib
import os
import sys
import pickle


class upeptide_mapper_1_0_0( ursgal.UNode ):
    """
        upeptide_mapper_1_0_0 UNode
    """
    META_INFO = {
        'engine_type' : {
            'search_engine' : False,
            'converter'     : True
        },
        'output_extension'  : '.csv',
        'output_suffix'     : 'pmap',
        'input_types'       : ['.csv'],
        'include_in_git'    : True,
        'in_development'    : True,
        'utranslation_style': 'upeptide_mapper_style_1',
        'engine': {
            'platform_independent' : {
                'arc_independent' : {
                    'exe' : 'upeptide_mapper_1_0_0.py',
                },
            },
        },
        'citation'          : 'Kremer, L. P. M., Leufken, J., '\
            'Oyunchimeg, P., Schulze, S. & Fufezan, C. (2016) '\
            'Ursgal, Universal Python Module Combining Common Bottom-Up '\
            'Proteomics Tools for Large-Scale Analysis. '\
            'J. Proteome res. 15, 788-794.',
    }

    def __init__(self, *args, **kwargs):
        super(upeptide_mapper_1_0_0, self).__init__(*args, **kwargs)

    def _execute( self ):
        '''
        Peptide from search engine csv file are mapped to the given database(s)

        '''
        print('[ -ENGINE- ] Executing conversion ..')
        self.time_point(tag = 'execution')
        upeptide_mapper_main = self.import_engine_as_python_function()
        if self.params['output_file'].lower().endswith('.csv') is False:
            raise ValueError('Trying to unify a non-csv file')

        output_file = os.path.join(
            self.params['output_dir_path'],
            self.params['output_file']
        )
        input_file  = os.path.join(
            self.params['input_dir_path'],
            self.params['input_file']
        )

        scan_rt_lookup_path = self.meta_unodes['ucontroller'].scan_rt_lookup_path

        # find the last search/denovo engine:
        last_engine = self.get_last_engine(
            history = self.stats['history'],
        )


        tmp_files = upeptide_mapper_main(
            input_file      = input_file,
            output_file     = output_file,
            params          = self.params,
            search_engine   = last_engine,
        )
        for tmp_file in tmp_files:
            self.created_tmp_files.append(tmp_file)

        self.print_execution_time(tag='execution')
        return output_file
