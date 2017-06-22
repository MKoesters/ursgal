#!/usr/bin/env python3.4
import ursgal
import os
import csv
import pymzml


class metfrag_2_3_1(ursgal.UNode):
    """
    MetFrag v2.31 UNode.

    Get CL version at http://msbi.ipb-halle.de/~cruttkie/metfrag/MetFrag2.3.1-CL.jar
    """
    META_INFO = {
        'edit_version'                : 1.00,
        'name'                        : 'metfrag',
        'version'                     : 'v2.31',
        'release_date'                : '2016',
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
        from pprint import pprint
        #pprint(self.params)
        pprint(self.params['translations'])
        peak_input_file = self._parse_mzml(self.params['input_file_incl_path'])
        param_file      = self._write_param_file()
        self.params['command_list'] = [
            'java',
            '-Xmx{0}'.format(self.params['translations']['-xmx']),
            '-jar',
            self.exe,
            param_file
        ]


    def postflight(self):
        """DOC."""
        # convert output from xls to csv
        # unify output format
        pass

    def _parse_mzml(path):
        """DOC."""
        with open('PeakList.txt', 'wt') as fout:
            specReader = pymzml.run.Reader(path)
            for spec in specReader:
                if spec.ms_level == 2:
                    #for mz, i in spec.peaks('centroided'):
                    for mz, i in spec.centroidedPeaks:
                        fout.write('{mz}\t{i}'.format(mz=mz, i=i))
        return fout

    def _write_param_file():
        with open('MetFrag_param_file', 'wt') as fout:
            for key, value in self.params['translations']:
                fout.write('{key:<15} = {value:<10}'.format(key=key, value=value))
        return fout



