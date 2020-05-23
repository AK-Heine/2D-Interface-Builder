from interfacebuilder.misc import *
import networkx as nx
import spglib
from collections import namedtuple
import ase.build


def recenter(atoms):
    """ Recenters atoms to be in the unit cell, with vacuum on both sides.

    The unit cell length c is always chosen such that it is larger than a and b.
    
    Returns:
        atoms : modified atoms object.

    Note:
        The ase.atoms.center() method is supposed to do that, but sometimes separates the layers. I didn't find a good way to circumvene that.

    """
    atoms = atoms.copy()
    atoms.wrap(pretty_translation=True)
    atoms.center(axis=(2))
    mp = atoms.get_center_of_mass(scaled=False)
    cp = (atoms.cell[0] + atoms.cell[1] + atoms.cell[2]) / 2
    pos = atoms.get_positions(wrap=False)
    pos[:, 2] += np.abs((mp - cp))[2]
    for z in range(pos.shape[0]):
        lz = atoms.cell.lengths()[2]
        if pos[z, 2] >= lz:
            pos[z, 2] -= lz
        if pos[z, 2] < 0:
            pos[z, 2] += lz
    atoms.set_positions(pos)
    newcell, newpos, newscal, numbers = (
        atoms.get_cell(),
        atoms.get_positions(wrap=False),
        atoms.get_scaled_positions(wrap=False),
        atoms.numbers,
    )
    z_pos = newpos[:, 2]
    span = np.max(z_pos) - np.min(z_pos)
    newcell[0, 2] = newcell[1, 2] = newcell[2, 0] = newcell[2, 1] = 0.0
    newcell[2, 2] = span + 100.0
    axes = [0, 1, 2]
    lengths = np.linalg.norm(newcell, axis=1)
    order = [x for x, y in sorted(zip(axes, lengths), key=lambda pair: pair[1])]
    while True:
        if (order == [0, 1, 2]) or (order == [1, 0, 2]):
            break
        newcell[2, 2] += 10.0
        lengths = np.linalg.norm(newcell, axis=1)
        order = [x for x, y in sorted(zip(axes, lengths), key=lambda pair: pair[1])]
    newpos = newscal @ newcell
    newpos[:, 2] = z_pos
    atoms = ase.Atoms(positions=newpos, numbers=numbers, cell=newcell, pbc=atoms.pbc)
    return atoms


def stack_atoms(atom1, atom2, weight=0.5, distance=4):
    """ Stacks two layered structures on top of each other.
    
    Args:
        atom1 (atoms): Lower layer.
        atom2 (atoms): Upper layer.
        weight (float, optional): Value between 0 and 1, defaults to 0.5. The unit cell of the reconstructed stack is :math:`B + w \cdot (T - B)`.
        distance (int, optional): Interlayer distance in Angström. Defaults to 4.
    
    Returns:
        atoms: Reconstructed stack.

    """

    bottom = atom1.copy()
    top = atom2.copy()
    c1 = np.linalg.norm(bottom.cell[2])
    c2 = np.linalg.norm(top.cell[2])
    cell1 = bottom.cell.copy()
    cell2 = top.cell.copy()
    cell1[2] /= c1
    cell2[2] /= c2
    cell = cell1 + weight * (cell2 - cell1)
    cell[2] /= np.linalg.norm(cell[2])
    cell1 = cell.copy()
    cell2 = cell.copy()
    cell1[2] *= c1
    cell2[2] *= c2

    bottom.set_cell(cell1, scale_atoms=True)
    top.set_cell(cell2, scale_atoms=True)

    zeroshift = np.min(bottom.get_positions()[:, 2])
    bottom.translate([0, 0, -zeroshift])
    zeroshift = np.min(top.get_positions()[:, 2])
    top.translate([0, 0, -zeroshift])
    bottom_thickness = np.max(bottom.get_positions()[:, 2]) - np.min(
        bottom.get_positions()[:, 2]
    )
    top.translate([0, 0, bottom_thickness])
    top.translate([0, 0, distance])
    bottom.extend(top)
    stack = recenter(bottom)
    return stack


class structure:
    """ A base class for structure analysis and manipulation relying on the ASE.
    
    Args:
        geometry (str): Path to structure file (.cif, .xyz ..) or atoms object.

    Attributes:
        atoms (Atoms): ASE atoms object.
        atom_indices (dict): Dictionary of atom index and label.
        species (dict): Dictionary of atom labels and counts.
        sg (spacegroup): Spglib spacegroup object.
        lattice (str): Description of Bravais lattice.
    """

    def __init__(self, geometry):
        if type(geometry) == ase.atoms.Atoms:
            self.atoms = geometry
        elif Path(geometry).is_file():
            try:
                self.atoms = ase.io.read(geometry)
            except:
                logging.error("Input structure not recognised.")
        assert type(self.atoms) == ase.atoms.Atoms, "Atoms not read correctly."

    def __repr__(self):
        return self.atoms.get_chemical_formula()

    def find_fragments(self, atoms=None):
        """ Finds unconnected structural fragments by constructing
        the first-neighbor topology matrix and the resulting graph
        of connected vortices. 

        Args:
            atoms (atoms): ASE atoms object.
        
        Note:
            Requires networkx library.
        
        Returns:
            dict: Dictionary with named tuples of indices and atoms, sorted by average z-value.

        """
        if atoms == None:
            atoms = self.atoms.copy()
        else:
            atoms = atoms.copy()
        nl = neighborlist.NeighborList(
            ase.neighborlist.natural_cutoffs(atoms),
            self_interaction=False,
            bothways=True,
        )
        nl.update(atoms)
        connectivity_matrix = nl.get_connectivity_matrix(sparse=False)

        con_tuples = {}  # connected first neighbors
        for row in range(connectivity_matrix.shape[0]):
            con_tuples[row] = []
            for col in range(connectivity_matrix.shape[1]):
                if connectivity_matrix[row, col] == 1:
                    con_tuples[row].append(col)

        pairs = []  # cleaning up the first neighbors
        for index in con_tuples.keys():
            for value in con_tuples[index]:
                if index > value:
                    pairs.append((index, value))
                elif index <= value:
                    pairs.append((value, index))
        pairs = set(pairs)

        graph = nx.from_edgelist(pairs)  # converting to a graph
        con_tuples = list(nx.connected_components(graph))

        fragments = {}
        i = 0
        for tup in con_tuples:
            fragment = namedtuple("fragment", ["indices", "atoms"])
            ats = ase.Atoms()
            indices = []
            for entry in tup:
                ats.append(atoms[entry])
                indices.append(entry)
            ats.cell = atoms.cell
            ats.pbc = atoms.pbc
            fragments[i] = fragment(indices, ats)
            i += 1
        fragments = {
            k: v
            for k, v in sorted(
                fragments.items(),
                key=lambda item: np.average(item[1][1].get_positions()[:, 2]),
            )
        }
        return fragments

    def standardize(self, to_primitive=True, symprec=1e-4):
        """ Wrapper of the spglib standardize() function with extra features.

        For 2D systems, the non-periodic axis is enforced as the z-axis.
        
        Args:
            to_primitive (bool): Reduces to primitive cell or not.
            symprec (float): Precision to determine new cell.

        Note:
            The combination of to_primitive=True and a larger value of symprec (1e-3) can be used to symmetrize a structure.
        """
        atoms = self.atoms.copy()
        pbc1 = self.find_nonperiodic_axes()
        lattice, positions, numbers = (
            atoms.get_cell(),
            atoms.get_scaled_positions(),
            atoms.numbers,
        )
        cell = (lattice, positions, numbers)
        newcell = spglib.standardize_cell(
            cell, to_primitive=to_primitive, no_idealize=False, symprec=symprec
        )
        if newcell == None:
            logging.error("Cell could not be standardized.")
            return None
        else:
            atoms = ase.Atoms(
                scaled_positions=newcell[1],
                numbers=newcell[2],
                cell=newcell[0],
                pbc=self.atoms.pbc,
            )
            pbc2 = self.find_nonperiodic_axes()
            if pbc1 != pbc2:
                old = [k for k, v in pbc1.items() if v]
                new = [k for k, v in pbc2.items() if v]
                assert len(old) == len(
                    new
                ), "Periodicity changed due to standardization."
                if len(new) == 2:
                    npbcax = list(set([0, 1, 2]) - set(new))[0]
                    atoms = ase.geometry.permute_axes(atoms, new + [npbcax])
                    assert self.is_2d(atoms), "Permutation to 2D not working."
            self.atoms = atoms

    def is_2d(self):
        """ Evaluates if given structure is qualitatively two-dimensional.

        Note:
            A 2D structure is considered 2D if only the z-axis is non-periodic.
        
        Returns:
            bool: 2D or not to 2D, that is the question.
        """
        atoms = self.atoms.copy()
        if list(atoms.pbc) == [True, True, False]:
            return True
        else:
            pbcax = self.find_nonperiodic_axes()
            if list(pbcax.values()) == [True, True, False]:
                return True
            else:
                return False

    def find_nonperiodic_axes(self):
        """ Evaluates if given structure is continuous along certain lattice directions.

        Returns:
            dict: Axis : Bool pairs.
        """

        atoms = self.atoms.copy()

        sc = ase.build.make_supercell(atoms, 2 * np.identity(3), wrap=True)
        fragments = self.find_fragments(sc)
        crit1 = True if len(fragments) > 1 else False
        pbc = dict(zip([0, 1, 2], [True, True, True]))
        if crit1:
            for axes in (0, 1, 2):
                spans = []
                for index, tup in fragments.items():
                    start = np.min(tup[1].get_positions()[:, axes])
                    cm = tup[1].get_center_of_mass()[axes]
                    spans.append([start, cm])
                spans = sorted(spans, key=lambda x: x[0])
                for j in range(len(spans) - 1):
                    nd = spans[j + 1][1] - spans[j][1]
                    if nd >= 30.0:
                        pbc[axes] = False
                        break
        return pbc
