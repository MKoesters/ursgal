#!/usr/bin/env python3
"""
DOC
"""
import ursgal
import os
import sys

def main(path):
    """Run a simple test search using MetFrag 2.31."""
    uc = ursgal.UController()
    uc.params['MetFrag_database_type'] = 'Local SDF'
    uc.params['database']              = os.path.join(os.pardir, 'example_data', 'NucleosideDB.sdf.gz')
    print(path)
    path = uc.execute_unode(
        input_file=path,
        engine='metfrag_2_3_1',
    )
    print(path)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
