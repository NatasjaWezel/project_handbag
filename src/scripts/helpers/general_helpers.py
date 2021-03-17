import os
import csv


def check_if_label_exists(atom, fragment):
    """ Checks if label already exists in the fragment. Adds the letter 'a' to it if it does. """

    if atom.label in fragment.atoms.keys():
        atom.label += 'a'
        atom = check_if_label_exists(atom, fragment)

    return atom


def split(filehandler, delimiter=',', row_limit=10000,
          output_name_template='output_%s.csv', output_path='.'):
    """
    Splits a CSV file into multiple pieces.

    A quick bastardization of the Python CSV library.
    Arguments:
        `row_limit`: The number of rows you want in each output file
        `output_name_template`: A %s-style template for the numbered output files.
        `output_path`: Where to stick the output files
        `keep_headers`: Whether or not to print the headers in each output file.
    Example usage:

        >> from toolbox import csv_splitter;
        >> csv_splitter.split(csv.splitter(open('/home/ben/Desktop/lasd/2009-01-02 [00.00.00].csv', 'r')));

    """

    reader = csv.reader(filehandler, delimiter=delimiter)
    current_piece = 1
    current_out_path = os.path.join(output_name_template % current_piece)

    current_out_writer = csv.writer(open(current_out_path, 'w', newline=''))
    current_limit = row_limit

    for i, row in enumerate(reader):
        if i + 1 > current_limit:
            print(f"{current_piece}th file is done")

            current_piece += 1
            current_limit = row_limit * current_piece
            current_out_path = os.path.join(output_name_template % current_piece)

            current_out_writer = csv.writer(open(current_out_path, 'w', newline=''))

        current_out_writer.writerow(row)

    print(f"Splitted input into {current_piece} different files")
