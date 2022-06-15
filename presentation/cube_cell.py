from visualize import *

# Create atom objects from sets of points
atoms = [Atom([0, 0, 0], 1), Atom([5, 0, 0], 1), Atom([-5, 0, 0], 1), Atom([0, 5, 0], 1), Atom([0, -5, 0], 1),
         Atom([0, 0, 5], 1), Atom([0, 0, -5], 1)]

# Create a system of from the atoms
selfsys = System()
selfsys.atoms = atoms

# Plot the atom objects
plot_atoms(selfsys.atoms, colors=['r', 'b'])
