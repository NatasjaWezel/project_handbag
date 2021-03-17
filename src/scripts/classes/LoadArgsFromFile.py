# # # # # # # # # # you know what to do
# https://stackoverflow.com/questions/27433316/how-to-get-argparse-to-read-arguments-from-a-file-with-an-option-rather-than-pre
# # # #
import argparse


class LoadArgsFromFile(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        with values as f:
            contents = f.readlines()

        contents = [x.strip() for x in contents]

        contents = " ".join(contents)

        # parse arguments in the file and store them in a blank namespace
        data = parser.parse_args(contents.split(), namespace=None)

        for k, v in vars(data).items():
            # set arguments in the target namespace if they haven’t been set yet
            if getattr(namespace, k, None) is None:
                setattr(namespace, k, v)
