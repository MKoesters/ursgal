#!/usr/bin/env python3.4
import ursgal
import os
import pymzml
from pprint import pprint

class metfrag_2_3_1(ursgal.UNode):
    """
    MetFrag v2.31 UNode.

    Get CL version at http://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.3.1-CL.jar
    """
    META_INFO = {
        'edit_version'                : 1.00,
        'name'                        : 'metfrag',
        'version'                     : '2_3_1',
        'release_date'                : '2016-01-01',
        'engine_type' : {
            'search_engine' : False,
        },
        'input_extensions'            : ['.mgf', '.mzML', '.txt'], # simple mz i pairs in one row as input
        'input_multi_file'            : False,
        'output_extensions'           : ['.csv'],
        'compress_raw_search_results' : False,
        'create_own_folder'           : True,
        'in_development'              : True,
        'include_in_git'              : False,
        'utranslation_style'          : 'metfrag_style_1',
        'engine' : {
            'platform_independent' : {
                'arc_independent' : {
                    'exe'            : 'MetFrag2.3.1-CL.jar',
                    'url'            : '',
                    'zip_md5'        : '',
                    'additional_exe' : [],
                },
            },
        },
        'citation' : """
                     MetFrag relaunched: incorporating strategies beyond in silico fragmentation: 
                     C Ruttkies, E L Schymanski, S Wolf, J Hollender, S Neumann
                     Journal of Cheminformatics 2016 8:3
                     """
    }

    def __init__(self, *args, **kwargs):
        """Init paren object."""
        super(metfrag_2_3_1, self).__init__(*args, **kwargs)

    def preflight(self):
        """Doc."""
        input_file      = os.path.join(
            self.params['input_dir_path'],
            self.params['input_file']
        )
        output_file     = os.path.join(
            self.params['output_dir_path'],
            self.params['output_file']
        )
        if input_file.lower().endswith('.mzml'):
            peak_file = self._write_peak_list(input_file)
        else:
            peak_file = input_file
        param_file      = self._write_param_file(
            self.params,
            peak_file,
            output_file
        )
        self.params['command_list'] = [
            'java',
            '-Xmx{0}'.format(self.params['-xmx']),
            '-jar',
            self.exe,
            param_file
        ]

    def postflight(self):
        """DOC."""
        # convert output from xls to csv
        # unify output format
        pass

    def _write_peak_list(self, path):
        """DOC."""
        peak_file_path = os.path.join(
            self.params['output_dir_path'],
            'peak_list.txt'
        )
        with open(peak_file_path, 'wt') as fout:
            specReader = pymzml.run.Reader(path)
            for spec in specReader:
                # if spec.ms_level == 2:
                if spec['ms level'] == 2:
                    #for mz, i in spec.peaks('centroided'):
                    for mz, i in spec.centroidedPeaks:
                        fout.write('{mz}\t{i}\n'.format(mz=mz, i=i))
                    break
        return peak_file_path

    def _write_param_file(self, params, peak_file, output_path):
        param_output_path = os.path.join(
            os.path.dirname(
                self.exe
            ),
            'params.txt'
        )
        with open(param_output_path, 'wt') as fout:
            for key in self.params['translations']:
                if key.endswith('_key') and key[0] not in ['_', '-']:
                    tr_key = self.params['translations'][key]
                    str_key = key.replace('_key', '')
                    val    = self.params['translations']['_grouped_by_translated_key'][tr_key][str_key]
                    if isinstance(val, list):
                        val = ' '.join(val)
                    if key and val:
                        fout.write('{key:<40} = {val:<10}\n'.format(
                            key=tr_key,
                            val=val
                            )
                        )
            fout.write('{key:<40} = {val:<10}\n'.format(
                key='PeakListPath',
                val=peak_file
                )
            )
            fout.write('{key:<40} = {val:<10}\n'.format(
                key='SampleName',
                val=os.path.splitext(
                    os.path.basename(output_path)
                )[0]
                )
            )
            fout.write('{key:<40} = {val:<10}\n'.format(
                key='ResultsPath',
                val=os.path.dirname(output_path)
                )
            )
        return param_output_path
