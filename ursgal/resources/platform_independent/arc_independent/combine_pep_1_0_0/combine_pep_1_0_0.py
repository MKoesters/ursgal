#!/usr/bin/python3.4

import csv
import os
import functools
import operator
import itertools
import collections


def sliding_window(iterable, window_size):
    '''Sliding window generator'''
    if window_size % 2 == 0:
        print('Warning! Window size must be uneven (to determine a '\
              'central value). Adjusted window size from {0} to {1}'\
              '.'.format(window_size, window_size+1))
        window_size += 1

    iterable = iter(iterable)
    half_window_size = int((window_size-1)/2)
    win = collections.deque(maxlen=window_size)

    try:

        # collect first bin (half the size of window_len)
        for __ in range(half_window_size):
            try:
                win.append(iterable.__next__())
            except StopIteration:
                continue

        # the first few values cannot get full bins since there are few values
        # to the left, so their bins are smaller and contain fewer values to the
        # left:
        iter_terminated = False
        for i in range(half_window_size):
            try:
                win.append(iterable.__next__())
            except StopIteration:
                # this would terminate the whole function...
                iter_terminated = True
                pass
            yield win, win[i]

        if iter_terminated:
            yield win, win[i+1]

        # windows from the center of the iterable, all have the full-length:
        for e in iterable:
            win.append(e)
            yield win, win[half_window_size]

        # the last few values also cannot get full bins since there are few values
        # to the right, so their bins are smaller and contain fewer values to the
        # right:
        for __ in range(half_window_size):
            win.popleft()
            yield win, win[half_window_size]

    except IndexError:
        raise StopIteration


class CombinedPEP(object):
    ''''''
    def __init__(self, *args):
        self.psm_dicts = {}
        self.score_dict = {}
        self.input_csv_fieldnames = set()

    @staticmethod
    def row_is_decoy(row):
        '''
        Check if a unified CSV row is a target or decoy PSM.
        Returns True if decoy and False if target.
        '''
        if row['Is decoy'].lower() == 'true':
            is_decoy = True
        elif row['Is decoy'].lower() == 'false':
            is_decoy = False
        else:
            raise Exception('Could not determine whether PSM is decoy or not.')
        if 'decoy_' in row.get('proteinacc_start_stop_pre_post_;', ''):
            is_decoy = True
        return is_decoy

    @staticmethod
    def naive_bayes(prob_list):
        '''
        Combines independent probabilities a, b, c
        like this:

                                a*b*c
        combined_prob = -------------------------
                        a*b*c + (1-a)*(1-b)*(1-c)

        For a straightforward explanation, see
        http://www.paulgraham.com/naivebayes.html
        '''
        # multiplying all probabilities: a*b*c
        multiplied_probs = \
            functools.reduce(operator.mul, prob_list, 1)
        # multiplying the opposite of all probabilities: (1-a)*(1-b)*(1-c)
        multiplied_opposite_probs = \
            functools.reduce(operator.mul, (1-p for p in prob_list), 1)
        return multiplied_probs / (multiplied_probs + multiplied_opposite_probs)

    def add_engine_result_csv(self, input_csv_path, input_engine_name):
        with open(input_csv_path, 'r', encoding='utf8') as csv_obj:
            reader = csv.DictReader(csv_obj)
            for fieldname in reader.fieldnames:
                self.input_csv_fieldnames.add(fieldname)
            self.psm_dicts[input_engine_name] = self.reader_to_psm_dict(
                reader, input_engine_name
            )
        print('Finished parsing {0}'.format(os.path.basename(input_csv_path)))

    def reader_to_psm_dict(self, reader, input_engine_name):
        '''
        Turns a CSV reader object into a dictionary.
        Dict key is a tuple of column values that can be specified by the user
        (usually a combinations of Seq/Mods/Spectrum).
        Dict value is the DictReader row.
        '''
        psm_dict = {}
        for row in reader:
            row_key = self.get_row_key(row, self.columns_for_grouping)
            psm_dict[row_key] = row
        return psm_dict

    @staticmethod
    def list_to_sorted_tuple(l):
        return tuple(sorted(set(l)))

    @staticmethod
    def get_row_key(row, columns_for_grouping):
        return tuple( row[c] for c in columns_for_grouping )

    @staticmethod
    def all_combinations(some_iterable):
        for i in range(len(some_iterable)):
            for combo in itertools.combinations(some_iterable, i+1):
                yield tuple(sorted(combo))

    def generate_psm_to_scores_dict(self, input_engines):
        for engine_combo in self.all_combinations(input_engines):
            engines_not_in_combo = {e for e in input_engines if e not in engine_combo}

            print("\nCombo:", engine_combo)
            print("not in Combo:", engines_not_in_combo)

            self.score_dict[engine_combo] = {}

            # get all shared PSMs:
            all_PSMs_of_combo_engines = set.intersection(
                *(set(self.psm_dicts[engine].keys()) for engine in engine_combo)
            )
            # remove all PSMs that are found by other engines:
            for other_eng in engines_not_in_combo:
                all_PSMs_of_combo_engines -= set(self.psm_dicts[other_eng].keys())
            print('# of PSMs in this combo:', len(all_PSMs_of_combo_engines))

            # For every PSM, use naive Bayes to calculate the combined PEP
            # ('Bayes PEP') among all engines and add it to self.score_dict:

            decoy_count_of_intersection = 0
            psm_count_of_intersection = 0
            for PSM_key in all_PSMs_of_combo_engines:
                engine_PEPs = []
                engine_decoy_bools = set()
                for eng in engine_combo:
                    d = self.psm_dicts[eng]
                    engine_PEPs.append(float(d[PSM_key]['PEP']))
                    engine_decoy_bools.add(self.row_is_decoy(d[PSM_key]))

                if len(engine_decoy_bools) != 1:
                    print('Engines don\'t agree whether PSM {0} is decoy or not'\
                          '!'.format(PSM_key))
                    PSM_is_decoy = True
                else:
                    PSM_is_decoy = list(engine_decoy_bools)[0]
                
                psm_count_of_intersection += 1
                if PSM_is_decoy:
                    decoy_count_of_intersection += 1
                    
                bayesPEP = self.naive_bayes(engine_PEPs)

                self.score_dict[engine_combo][PSM_key] = {
                    'Bayes PEP' : bayesPEP,
                    'Is decoy' : PSM_is_decoy,
                }

            print('Decoy proportion: {:.3%}%'.format(
                decoy_count_of_intersection / psm_count_of_intersection
            ))


            # Use the Bayes PEP to sort all PSMs that are shared by the current
            # combination of engines. Compute the PEP (similar to localized FDR)
            # for each PSM and store it in score_dict:

            sorted_score_dict_items = sorted(
                self.score_dict[engine_combo].items(),
                key = lambda kv: kv[1]['Bayes PEP']
            )

            sliding_window_rows = sliding_window(
                sorted_score_dict_items, self.window_size)

            for window_tup, central_tup in sliding_window_rows:
                n_false_positives = 2 * sum(
                    key_val[1]['Is decoy'] for key_val in window_tup
                )
                n_PSMs = len(window_tup)
                
                ##assertion, used for debugging! This makes sure that only PMSs
                ##from the correct combination of engines are used!
                #for w in window_tup:
                    #row_key = w[0]
                    #for eng in engine_combo:
                        #assert row_key in self.psm_dicts[eng], '{0} is not in {1}'.format(row_key, eng)

                intersection_PEP = n_false_positives / n_PSMs
                current_PSM_key = central_tup[0]
                self.score_dict[engine_combo][current_PSM_key]['combined PEP'] = \
                    intersection_PEP
        return

    def write_output_csv(self, output_csv_path):
        new_scores = ['combined PEP', 'Bayes PEP']
        fieldnames = list(self.input_csv_fieldnames) + new_scores + ['engines']
        with open(output_csv_path, 'w', encoding='utf8') as out_obj:
            writer = csv.DictWriter(out_obj, fieldnames=fieldnames)
            writer.writeheader()

            for engine_combo, combo_score_dict in self.score_dict.items():

                for psm_key, score_dict_val in combo_score_dict.items():

                    psm_rows_from_all_engines = []
                    for engine in engine_combo:
                        psm_rows_from_all_engines.append(
                            self.psm_dicts[engine][psm_key]
                        )
                    out_row = self.merge_rowdicts(
                        psm_rows_from_all_engines,
                        sep = self.join_sep,
                    )

                    # add column that lists all engines that found the PSM:
                    out_row['engines'] = self.join_sep.join(engine_combo)
                    # add columns with Bayes PEP and combined PEP:
                    for score_field in new_scores:
                        out_row[score_field] = score_dict_val[score_field]
                    writer.writerow(out_row)
        return

    def merge_rowdicts(self, rows_to_merge, sep=';'):
        '''
        Merges different DictReader rows (=dictionaries) to a single merged
        dictionary. If all values are in agreement, only that value is entered.
        If there are conflicting values, they are joined with a separator
        (i.e. ';').
        '''
        md = {}  # merged rowdict
        for row in rows_to_merge:
            for col_name, field_val in row.items():

                if col_name not in md:
                    # value is not there yet, so add it
                    md[col_name] = field_val
                    continue

                if md[col_name] == field_val:
                    # value is already there, no need to add it
                    continue

                # Different value(s) are already there! join them (i.e. with ';')
                old_vals = set(md[col_name].split(sep))
                old_vals.add(field_val)
                old_vals.discard('')
                md[col_name] = sep.join(sorted(old_vals))
        return md

def main(columns_for_grouping=None, input_csvs=None, output_csv=None,
         input_sep=None, output_sep=None, join_sep=None, pep_colname=None,
         input_engines=None, window_size=None):
    '''
    1. Set parsed attributes from command line as class attribute.
    2. Parse unified input CSV files (with PEPs) and buffer them as
       PSM-to-row dictionaries.
    3. calculate Bayes PEP for each PSM
    4. retrieve all possible engine combinations (like a VennDiagram)
    5. For each engine combination:
        a) retrieve list of shared PSMs
        b) sort PSMs by Bayes PEP
        c) loop over sorted PSMs and calculate combined PEP (sliding window)
    6. Return merged CSV output file with added Bayes PEP and combined PEP
       columns.
    '''
    c = CombinedPEP()

    # Set parsed attributes from command line as class attribute:
    c.columns_for_grouping = c.list_to_sorted_tuple(columns_for_grouping)
    c.input_sep = input_sep
    c.output_sep = output_sep
    c.join_sep = join_sep
    c.pep_colname = pep_colname
    c.window_size = window_size

    # Parse input CSV files and buffer them as PSM-to-row dicts:
    for input_csv, input_engine in zip(input_csvs, input_engines):
        c.add_engine_result_csv(input_csv, input_engine)

    # Calculate Bayes PEP and combined PEP for each PSM:
    c.generate_psm_to_scores_dict(input_engines)

    # Write merged CSV with Bayes PEP and combined PEP column:
    c.write_output_csv(output_csv)
    return

def parse_args(verbose=True):
    '''
    Parses command line arguments and returns them in a dict.
    Only used when executing this script from command line.
    '''
    import argparse
    # no need to import argparse when script is executed by importing
    # main function.

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input_csvs', type=str, nargs='+', required=True,
        help='Paths to unified input CSV files (2 or more)')
    parser.add_argument(
        '-c', '--columns_for_grouping', type=str, nargs='+', required=True,
        help='Column names by which the rows should be grouped')
    parser.add_argument(
        '-o', '--output_csv', type=str, help='Output CSV name')
    parser.add_argument(
        '-is', '--input_sep', type=str, default=',',
        help='Input file column delimiter character')
    parser.add_argument(
        '-os', '--output_sep', type=str, default=',',
        help='Output file column delimiter character')
    parser.add_argument(
        '-js', '--join_sep', type=str, default=';',
        help='Delimiter for multiple values in the same field')
    parser.add_argument(
        '-e', '--input_engines', type=str, required=True,
        help='The search engines of each input file (must be same order as input_csvs)')
    parser.add_argument(
        '-w', '--window_size', type=int, default=251,
        help='The size of the sliding window for PEP calculation.')

    args = vars(parser.parse_args())  # convert to dict
    if verbose:
        print('You are using the following settings:')
        for arg, val in sorted(args.items()):
            print('  {0: <13}:  {1}'.format(arg, val))
        print()
    return args


if __name__ == '__main__':
    command_line_args = parse_args()
    main(**command_line_args)