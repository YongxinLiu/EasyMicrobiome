#!/usr/bin/env python3

'''
    slime2 -- synthetic learning in microbial ecology
    Copyright (C) 2015  Scott W. Olesen

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


    swo@alum.mit.edu
'''

import argparse, hashlib, sys, time, random, itertools, pickle
import pandas as pd, numpy as np
import sklearn.ensemble

def hash_tag(strs, tag_length=6):
    '''create a short tag using md5'''

    m = hashlib.md5()
    m.update((''.join([s for s in strs])).encode('utf-8'))
    return m.hexdigest()[0: tag_length]

def parse_table_and_classes(table, klasses_fn, normalize=False, logt=None, keep=False):
    '''
    read an OTU table and a file that specifies the classes for each sample.

    the table has qiime format (rows=otus, columns=samples).
    the classes file can take one of two formats: one-column or two-column.

    keep : don't throw away columns that are not in the klasses file

    returns : (trimmed_table, classes)
        trimmed_table : pandas dataframe, with the samples without specified classes
            being dropped
        classes : list of classes in the same order as the indices in the dataframe
    '''

    # read in the classes files
    with open(klasses_fn) as f:
        lines = [l.rstrip() for l in f if not l.startswith("#")]

    # figure out how to parse the classes file
    n_fields = len(lines[0].split())
    if n_fields == 1:
        samples, klasses = parse_one_column_klasses(lines)
    elif n_fields == 2:
        samples, klasses = parse_two_column_klasses(lines)
    else:
        raise RuntimeError("got {} fields in classes file".format(n_fields))

    # read in the table
    raw_table = pd.read_table(table, index_col=0).transpose()

    # check to see that all the samples are columns
    # complain if the classes file gave a sample not that's not in the OTU table
    cols = list(raw_table.index.values)
    missing_cols = [s for s in samples if s not in cols]
    if len(missing_cols) > 0:
        raise RuntimeError("samples {} not a column in table, which has columns {}".format(missing_cols, cols))

    if keep:
        # figure out which rows/OTUs _would_ be removed if we were going to toss out samples
        tmp_table = raw_table.loc[list(samples)]
        trim_table = raw_table.loc[:, (tmp_table.sum(axis=0) != 0)]
    if not keep:
        # only keep samples in the OTU table if they have classes associated with them
        trim_table = raw_table.loc[list(samples)]

        # remove OTUs that have all-zero counts
        trim_table = trim_table.loc[:, (trim_table.sum(axis=0) != 0)]

    # if doing a log transformation, add the pseudocounts and proceed
    if logt is not None:
        trim_table += logt
        trim_table = np.log(trim_table)

    # if normalizing, do that
    if normalize:
        trim_table = trim_table.apply(lambda col: col.astype(float) / np.sum(col), axis=0)

    return trim_table, klasses

def parse_two_column_klasses(lines):
    '''
    two-column files have lines with two tab-separated fields: sample-tab-class.
    for example, lines would be like sick_guy1 tab sick, healthy_guy1 tab healthy, etc.
    '''

    samples, klasses = zip(*[line.split('\t') for line in lines])
    return samples, klasses

def parse_one_column_klasses(lines, comment="#"):
    '''
    one-column files have a header line with the name of the class, then the
    samples in that class, then a blank line before the next class. for example,
    lines would be like: sick, sick_guy1, sick_guy2, blank, healthy, healthy_guy1, etc.
    '''

    samples = []
    klasses = []

    read_klass = True # flag for asking if the next non-comment line is a class name
    for line in lines:
        if line.startswith(comment):
            # ignore comment lines
            continue
        elif line == "":
            # the next line after a blank is a class name
            read_klass = True
        elif read_klass:
            # the next lines after the class are samples
            klass = line
            read_klass = False
        else:
            sample = line
            samples.append(sample)
            klasses.append(klass)

    return samples, klasses


def create_classifier(constructor, otu_table, klasses, constructor_kwargs, sample_weights=None):
    '''
    initialize classifier from otu table and class file

    clf is the classifier constructor, e.g., sklearn.ensemble.RandomForestClassifier
    '''

    classifier = constructor(**constructor_kwargs)
    classifier.fit(table, klasses, sample_weights)

    # attach some extra data to the object for bookkeeping
    classifier.true_klasses = klasses
    classifier.predicted_klasses = classifier.predict(table)
    classifier.total_score = classifier.score(table, klasses)
    classifier.feature_names = list(table.columns.values)
    classifier.ordered_features = sorted(zip(classifier.feature_names, classifier.feature_importances_), key=lambda x: -x[1])

    return classifier


def tagged_name(fn, tag, suffix='txt'):
    '''format filenames for classifier output'''
    return "{}_{}.{}".format(tag, fn, suffix)

def categorize_classifications(targets, predictions):
    '''
    this is like an explicit confusion matrix. make a line for every sample.
    if a sample was X and was classified as X, just write "--". if it was X
    but classified as Y, write ">> X misclassified as Y".
    '''

    out = []
    for t, p in zip(targets, predictions):
        if t == p:
            out.append("--")
        else:
            out.append(">> {} misclassified as {}".format(t, p))

    return out

def save_results(classifier, tag):
    '''
    save the results of a classifier run. tag every output filename with a prefix
    so that they all appear next to each other in the directory.
    '''

    # pickle the whole classifier
    with open(tagged_name('classifier', tag, suffix='pkl'), 'wb') as f:
        pickle.dump(classifier, f, protocol=pickle.HIGHEST_PROTOCOL)

    # save the other information in text files
    with open(tagged_name('classes', tag), 'w') as f:
        f.write('\n'.join(classifier.predicted_klasses) + '\n')

    with open(tagged_name('featimp', tag), 'w') as f:
        cumul_imp = 0
        for of in classifier.ordered_features:
            cumul_imp += float(of[1])
            f.write("{}\t{}\t{:.3f}\n".format(of[0], of[1], cumul_imp))

    with open(tagged_name('scores', tag), 'w') as f:
        f.write("mean score: {}".format(classifier.total_score) + '\n')
        if hasattr(classifier, 'oob_score_'):
            f.write("oob score: {}".format(classifier.oob_score_) + '\n')

    with open(tagged_name('results', tag), 'w') as f:
        f.write('\n'.join(categorize_classifications(classifier.true_klasses, classifier.predicted_klasses)))

    with open(tagged_name('params', tag), 'w') as f:
        f.write('\n'.join(["{}: {}".format(*x) for x in classifier.get_params().items()]))

def int_or_none(x):
    '''
    take a string. if the string is 'none', return None object. if it's
    an integer, return that integer.
    '''

    assert(isinstance(x, str))
    if x.lower() == 'none':
        return None
    else:
        return int(x)

def int_float_str(x):
    '''
    take a string. if it's an integer, parse it that way. then try for a
    float. if that fails, just leave it as a string.
    '''

    assert(isinstance(x, str))
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return x

def assign_weights(weights_string, klasses):
    weight_vals = [float(x) for x in weights_string.split(",")]
    chunked_klasses = [g[0] for g in itertools.groupby(klasses)]
    weight_map = {k: w for w, k in zip(weight_vals, chunked_klasses)}

    weights = np.array([weight_map[k] for k in klasses])
    return weights


if __name__ == '__main__':
    p = argparse.ArgumentParser(description="slime2", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    g = p.add_argument_group('io')
    g.add_argument('otu_table')
    g.add_argument('classes', help='specifications of samples and their true classes')
    g.add_argument('--tag', '-t', default=None, help='tag for output data (default: use a hash tag)')
    g.add_argument('--shuffle', action='store_true', help='shuffle class labels?')
    g.add_argument('--weights', help='set of floats, comma separated')
    g.add_argument('--quiet', '-q', action='store_true', help='suppress stdout information?')

    g = p.add_argument_group('table preprocessing')
    g.add_argument('--normalize', action='store_true', help='normalize table of counts to rel. abunds.?')
    g.add_argument('--logtransform', default=None, type=int, help='log transform? if so, what pseudocount to add?')

    subparsers = p.add_subparsers(help="choose one classifier")

    sp = subparsers.add_parser("rf", help="random forest")
    sp.set_defaults(constructor=sklearn.ensemble.RandomForestClassifier)
    sp.set_defaults(constructor_kwarg_keys=['n_estimators', 'criterion', 'max_features', 'random_state', 'max_depth', 'oob_score', 'n_jobs', 'verbose'])
    sp.add_argument('--n_estimators', '-n', default=10, type=int, help='number of trees')
    sp.add_argument('--criterion', default='gini', choices=['gini', 'entropy'], help='function to measure quality of split')
    sp.add_argument('--max_features', '-f', type=int_float_str, default='auto')
    sp.add_argument('--random_state', '-r', type=int_or_none, default='none', help='random seed (none=random)')
    sp.add_argument('--max_depth', '-d', type=int_or_none, default='none', help='(none=no limit)')
    sp.add_argument('--no_oob_score', '-b', dest='oob_score', action='store_false')
    sp.add_argument('--n_jobs', '-j', type=int, default=1, help='-1=# of cores')
    sp.add_argument('--verbose', '-v', action='count', default=0, help='verbose output')

    sp = subparsers.add_parser("ab", help="adaboost")
    sp.set_defaults(constructor=sklearn.ensemble.AdaBoostClassifier)
    sp.set_defaults(constructor_kwarg_keys=['n_estimators', 'learning_rate', 'random_state'])
    sp.add_argument('--n_estimators', '-n', default=50, type=int, help='number of stumps')
    sp.add_argument('--learning_rate', '-a', default=1.0, type=float, help='shrink contribution of each classifier?')
    sp.add_argument('--random_state', '-r', type=int_or_none, default='none', help='random seed (none=random)')

    sp = subparsers.add_parser("load", help="load pickled classifier")
    sp.set_defaults(constructor="load")
    sp.add_argument('pickled_classifier', type=argparse.FileType('rb'))

    args = p.parse_args()

    # PREPROCESSING
    # generate the tag, unless supplied
    if args.tag is None:
        tag = hash_tag([open(args.otu_table).read(), open(args.classes).read()])
    else:
        tag = args.tag

    # save the command line
    with open(tagged_name('cmd', tag), 'w') as f:
        f.write(' '.join(sys.argv) + '\n')

    # CONSTRUCTION OF CLASSIFIER
    if args.constructor is None:
        raise RuntimeError('need to specify classifier on the command line')
    elif args.constructor == "load":
        table, klasses = parse_table_and_classes(args.otu_table, args.classes, normalize=args.normalize, logt=args.logtransform, keep=True)
        classifier = pickle.load(args.pickled_classifier)
        predicted_klasses = classifier.predict(table)
        prediction_probabilities = classifier.predict_proba(table)

        for sample, klass, prob in zip(table.index, predicted_klasses, prediction_probabilities):
            print(sample, klass, *prob, sep="\t")
    else:
        table, klasses = parse_table_and_classes(args.otu_table, args.classes, normalize=args.normalize, logt=args.logtransform, keep=False)

        # load the weights, if present
        if args.weights:
            sample_weights = assign_weights(args.weights, klasses)
        else:
            sample_weights = None

        # shuffle the classes, if requested
        if args.shuffle:
            random.shuffle(klasses)

        # extract the kwargs for the classifier's construction
        vargs = vars(args)
        constructor_kwargs = {k: vargs[k] for k in args.constructor_kwarg_keys}

        start_time = time.time()
        classifier = create_classifier(args.constructor, table, klasses, constructor_kwargs, sample_weights)

        if not args.shuffle:
            save_results(classifier, tag)
            print("saved results with tag {}".format(tag))

        end_time = time.time()

        if not args.quiet:
            print("walltime elapsed: {:.1f} seconds".format(end_time - start_time))

            if hasattr(classifier, 'oob_score_'):
                print("oob score: {:.5f}".format(classifier.oob_score_))

            print("top features:")
            for i in range(10):
                print("  {}\t{}".format(*classifier.ordered_features[i]))
