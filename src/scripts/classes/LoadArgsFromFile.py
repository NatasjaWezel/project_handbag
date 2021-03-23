# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script is part of the quantification pipeline of 3D experimental data of crystal structures that I wrote for my
# thesis in the Master Computational Science, University of Amsterdam, 2021.
#
# `LoadArgsFromFile` is a class that takes as input argparse action. It implements a function that checks if there is a
# file with input arguments, reads that if it exists, but does not overwrite other provided command line arguments.
#
# Inspired on
# https://stackoverflow.com/questions/27433316/how-to-get-argparse-to-read-arguments-from-a-file-with-an-option-rather-than-pre
#
# Author: Natasja Wezel
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import argparse


class LoadArgsFromFile(argparse.Action):
    """ Implements a call function to load arguments from a file, but does not overwrite command line arguments. """

    def __call__(self, parser, namespace, values, option_string=None):
        """ Call implements the functionality of reading command line arguments from a file. """

        # read all contents of the file
        with values as f:
            contents = f.readlines()

        # the contents are separated by an enter, so strip, and join with a space (to mimic command line arguments)
        contents = [x.strip() for x in contents]
        contents = " ".join(contents)

        # parse arguments in the file and store them in a blank namespace
        data = parser.parse_args(contents.split(), namespace=None)

        # loop over the namespaces in data
        for k, v in vars(data).items():

            # set arguments in the target namespace if they haven’t been set yet (with command line arguments)
            if getattr(namespace, k, None) is None:
                setattr(namespace, k, v)
