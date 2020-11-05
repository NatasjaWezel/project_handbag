def check_if_label_exists(atom, fragment):
    """ Checks if label already exists in the fragment. Adds the letter 'a' to it if it does. """

    if atom.label in fragment.atoms.keys():
        atom.label += 'a'
        atom = check_if_label_exists(atom, fragment)

    return atom
