#!/usr/bin/env python3.4
import ursgal
import os
import rdkit


class midas_1_1_0(ursgal.UNode):
    """
    MIDAS v1.1 UNode.
    """
    META_INFO = {
        'edit_version'                : 1.00,
        'name'                        : 'MIDAS',
        'version'                     : 'v1.1',
        'release_date'                : '2014',
        'engine_type' : {
            'search_engine' : True,
        },
        'input_extensions'            : ['.mgf', '.mzML'], # simple mz i pairs in one row as input
        'input_multi_file'            : False,
        'output_extensions'           : ['.csv'],
        'compress_raw_search_results' : False,
        'create_own_folder'           : True,
        'in_development'              : True,
        'include_in_git'              : False,
        'utranslation_style'          : 'midas_style_1',
        'engine': {
            'platform_independent': {
                'arc_independent' : {
                    'exe'            : 'midas_1_1.py',
                    'url'            : '',
                    'zip_md5'        : '',
                    'additional_exe' : [],
                },
            },
        },
        'citation' : "",
    }

    def __init__(self, *args, **kwargs):
        """Init paren object."""
        super(midas_1_1_0, self).__init__(*args, **kwargs)

    def preflight(self):
        """Doc."""
        input_file = os.path.join(
            self.params['input_dir_path'],
            self.params['input_file_path']
        )
        config_file = os.path.join(
            os.path.dirname(rdkit.__file__)
        )
        output_file = os.path.join(
            self.params['output_dir_path'],
            self.params['output_file_path']
        )
        annotation_file = os.path.join(
            self.params['output_dir_path'],
            os.path.basename(input_file).strip('.csv') + 'annotation.txt'
        )
        config_file = self._create_conf_file(
            config_file,
            self.params['translations']
        )
        self.params['command_list'] += [
            'python3',
            self.exe,
            '-f',
            input_file,
            '-c',
            config_file,
            '-o',
            output_file,
            '-a',
            annotation_file
        ]

    def postflight(self):
        """DOC."""
        pass

    def _create_conf_file(self, out_path, translated_params):
        """Write config input file."""
        with open(out_path, 'wt') as fout:
            for para, key_name in translated_params:
                if para.endswith('_key'):
                    val = translated_params[para.strip('_key')]
                    fout.write('{0}={1]\n'.format(key_name, val))
        return out_path
