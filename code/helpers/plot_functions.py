import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_fragments(fragments):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    colors = ["red", "sandybrown", "gold", "chartreuse", "green", 
                "mediumturquoise", "dodgerblue", "darkblue", "slateblue",
                "mediumorchid", "fuchsia"]

    for i, fragment in enumerate(fragments):
        fragment.color = colors[i % len(colors)]

        # plot first atom of the fragment and label it
        atom = list(fragment.atoms.values())[0]
        ax.scatter(atom.x, atom.y, atom.z, color=fragment.color, label=fragment.from_entry + str(fragment.fragment_id))
        ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black') 
        
        ax = plot_atoms_bonds(ax=ax, fragment=fragment)

    ax.legend()
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    
    plt.show()

def plot_atoms_bonds(ax, fragment):
    """ This function plots the atoms and bonds of a fragment. """ 

    for atom in list(fragment.atoms.values())[1:]:
        ax.scatter(atom.x, atom.y, atom.z, color=fragment.color)
        ax.text(atom.x + .005, atom.y + .005 , atom.z + .005,  atom.label, size=8, zorder=1, color='black')                 
    
    for bond in fragment.bonds:
        # TODO: look into this cuz fragments.atoms.keys are supposed to be labels and bond[0] is supposed to be an atom
        if bond[0] in fragment.atoms.keys() and bond[1] in fragment.atoms.keys():
            x = [fragment.atoms[bond[0]].x, fragment.atoms[bond[1]].x]
            y = [fragment.atoms[bond[0]].y, fragment.atoms[bond[1]].y]
            z = [fragment.atoms[bond[0]].z, fragment.atoms[bond[1]].z]

            ax.plot(x, y, z, color=fragment.color)
    
    return ax