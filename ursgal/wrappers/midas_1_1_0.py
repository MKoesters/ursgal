#!/usr/bin/env python3.4
import ursgal
import os


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
        pass

    def postflight(self):
        """DOC."""
        # convert output from xls to csv
        # unify output format
        pass
