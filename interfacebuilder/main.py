from interfacebuilder.misc import *
from interfacebuilder.structure import structure, recenter, stack_atoms
from interfacebuilder.plotting import interactive_plot
from interfacebuilder.backend import find_coincidence as _cxx_find_coincidence

import spglib


def _coincidence(a1, b1, m1, m2, n1, n2, theta):
    theta = theta / 180.0 * np.pi
    R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    Am = a1 @ np.array([m1, m2])
    Bn = R @ b1 @ np.array([n1, n2])
    return (m1, m2, n1, n2, theta, np.linalg.norm(Am - Bn))


def _standardize(collection, symprec=1e-4):
    atoms = collection.atoms
    cell = (atoms.get_cell()).tolist()
    pos = atoms.get_scaled_positions().tolist()
    numbers = atoms.get_atomic_numbers()

    cell, scaled_pos, numbers = spglib.standardize_cell(
        (cell, pos, numbers), to_primitive=True, symprec=symprec, no_idealize=False
    )

    atoms = ase.atoms.Atoms(
        scaled_positions=scaled_pos, numbers=numbers, cell=cell, pbc=True
    )
    axes = [0, 1, 2]
    lengths = atoms.cell.lengths()
    order = [x for x, y in sorted(zip(axes, lengths), key=lambda pair: pair[1])]
    if order != [0, 1, 2]:
        atoms = ase.geometry.permute_axes(atoms, order)
    collection._replace(atoms=atoms)
    return collection


class interface(interactive_plot):
    r"""Interface class to construct heterojunctions from two layered structures.

    Employs the algorithm outlined by Schwalbe-Koda (J. Phys. Chem. C 2016, 120, 20, 10895-10908) to find coincidence lattices of 2D crystals.

    Note:
        The scaling is :math:`\approx (2 \cdot N_{trans} + 1)^4 \cdot N_{angles}`. For the python backend, this is solved in by multiprocessing. The C++ backend is new and much faster.

    Example:
        >>> from interfacebuilder.main import interface
        >>> intf = interface(bottom, top)
        >>> intf.analyze_results(weight=0.5, distance=4)
        >>> intf.plot_results(jitter=0.05)
    
    Args:
        bottom (atoms): Input structure for the lower layer. Can be path to file or ase.atoms.Atoms object.
        top (atoms): Input structure for the upper layer. Can be path to file or ase.atoms.Atoms object.
    
    Keyword Args:
        N_translations (int): Number of translations or maximum super cell size to consider, defaults to 10.
        angle_stepsize (float): Stepsize between angles to rotate top layer. Defaults to 10.0.
        angle_limits (tuple): Lower and upper limit of angles in degrees to scan. Defaults to (0, 180).
        tolerance (float): Distance criterion to accept coincidence lattice points. Defaults to 0.1 (Angström).
        distance (float): Interlayer distance between the reconstructed stacks. Defaults to 4.0 (Angström).
        weight (float): Value between 0 and 1, defaults to 0.5. The unit cell of the reconstructed stack is :math:`B + w \cdot (T - B)`.
        prec (float): Precision to identify equivalent structures. Defaults to 1e-4.
        backend (str): Backend to perform calculations, either python or c++.
    
    Attributes:
        results (dict): Dictionary of angles in radians and list of lists containing (m1, m2, n1, n2, err).
        pairs (dict): Dictionary of angles in degrees and tuples of supercell matrices (M, N).
        solved (list): List of stacks as named tuples after analyzing the results. Contains (atoms, M, N, angle, stress).

    """

    def __init__(self, bottom, top, **kwargs):
        super().__init__()
        a1, b1 = self.__checks_and_setup(bottom, top)
        self.N_translations = kwargs.get("N_translations", 10)
        self.angle_stepsize = kwargs.get("angle_stepsize", 10.0)
        self.angle_limits = kwargs.get("angle_limits", (0, 180))
        self.tolerance = kwargs.get("tolerance", 0.1)
        self.distance = kwargs.get("distance", 4.0)
        self.weight = kwargs.get("weight", 0.5)
        self.prec = kwargs.get("prec", 1e-4)
        self.backend = kwargs.get("backend", "c++")

        self.results = self.find_coincidence(
            a1,
            b1,
            N=self.N_translations,
            stepsize=self.angle_stepsize,
            tolerance=self.tolerance,
            angle_limits=self.angle_limits,
            backend="c++",
        )
        if self.results != None:
            self.pairs = self.find_independent_pairs()

    def __checks_and_setup(self, bottom, top):
        if type(bottom) == ase.atoms.Atoms:
            self.bottom = structure(bottom.copy())
        elif Path(bottom).is_file():
            self.bottom = structure(Path(bottom))
        elif type(bottom) == structure:
            self.bottom = bottom
        else:
            logging.critical("Input for bottom structure not recognised.")
        if type(top) == ase.atoms.Atoms:
            self.top = structure(top.copy())
        elif Path(bottom).is_file():
            self.top = structure(Path(top))
        elif type(top) == structure:
            self.top = top
        else:
            logging.critical("Input for top structure not recognised.")
        assert self.bottom.is_2d(), "Bottom structure is not 2D!"
        assert self.top.is_2d(), "Top structure is not 2D!"
        self.bottom.standardize(to_primitive=True)
        self.top.standardize(to_primitive=True)
        a1 = self.bottom.atoms.cell[[0, 1], :2].T.copy()
        b1 = self.top.atoms.cell[[0, 1], :2].T.copy()
        return (a1, b1)

    def __py_find_coincidence(
        self,
        a1,
        b1,
        N=10,
        stepsize=1,
        tolerance=0.05,
        angle_limits=(0, 180.0),
        chunksize=10000,
        **kwargs
    ):
        import multiprocessing as mp

        chunksize = kwargs.get("chuncksize", 10000)
        start = time.time()
        cpus = mp.cpu_count()
        pool = mp.Pool(processes=cpus)
        nrange = range(-N, N + 1)
        angles = np.arange(angle_limits[0], angle_limits[1], stepsize)
        iterator = (
            (a1, b1, m1, m2, n1, n2, theta)
            for m1, m2, n1, n2, theta in itertools.product(
                nrange, nrange, nrange, nrange, angles
            )
        )
        ngrids = len(nrange) ** 4 * len(angles)
        logging.info("Initializing {} gridpoints ...".format(ngrids))
        logging.info("Parallelising on {} processors ...".format(cpus))
        res = pool.starmap_async(_coincidence, iterator, chunksize=chunksize)
        res = res.get()
        pool.close()
        pool.join()
        results = np.zeros((len(res), 6), dtype=np.float32)
        for i, j in enumerate(res):
            results[i] = j
        index = np.argwhere(
            (
                (
                    (
                        ((results[:, 0] == 0) & (results[:, 1] == 0))
                        | ((results[:, 2] == 0) & (results[:, 3] == 0))
                    )
                    | (results[:, 5] > tolerance)
                )
            )
            == False
        )
        results = results[index[:, 0]]
        return results

    def find_coincidence(self, a1, b1, **kwargs):
        """ Solves the system of linear equations given by :math:`\mathbf{A m = R B n}`.
    
        Args:
            a1 (ndarray): (2, 2) ndarray representing lower basis.
            b1 (ndarray): (2, 2) ndarray representing upper basis.            
            N (int, optional): Maximum supercell size. Defaults to 10.
            stepsize (int, optional): Angular step size. Defaults to 1 degree.
            tolerance (float, optional): Acceptance criterion. Defaults to 0.05.
            angle_limits (tuple, optional): Lower and upper limit of angular range. Defaults to (0, 180.0).
            backend (str, optional): Backend to perform calculations, either python or c++.
            chunksize (int): Distributed chunks in the asynchronous pool mapping if backend is python.
        
        Returns:
            dict: Pairs of angle : (m1, m2, n1, n2, err).

        """
        start = time.time()
        assert self.backend in ["python", "c++"], "Backend not accepted."
        if self.backend == "python":
            logging.info("Choosing python backend.")
            results = self.__py_find_coincidence(a1, b1, **kwargs)
        elif self.backend == "c++":
            logging.info("Choosing C++ backend.")
            a1t = list(a1[:, 0])
            a2t = list(a1[:, 1])
            b1t = list(b1[:, 0])
            b2t = list(b1[:, 1])
            angles = list(
                np.arange(
                    self.angle_limits[0], self.angle_limits[1], self.angle_stepsize
                )
                / 180
                * np.pi
            )
            N = self.N_translations
            tolerance = self.tolerance
            logging.info(
                "Initializing {} gridpoints ...".format(
                    ((2 * N + 1) ** 4) * len(angles)
                )
            )
            results = _cxx_find_coincidence(a1t, a2t, b1t, b2t, angles, N, tolerance)
            results = np.array(results, dtype=float)
            index = np.argwhere(
                (
                    (
                        ((results[:, 0] == 0) & (results[:, 1] == 0))
                        | ((results[:, 2] == 0) & (results[:, 3] == 0))
                    )
                )
                == False
            )
            results = results[index[:, 0]]
        end = time.time()
        logging.info(
            "Scanning angular and translational grid finished in {:.2f} seconds.".format(
                end - start
            )
        )
        if results.size > 2:
            logging.info("Found {:d} matching lattice points.".format(len(results)))
            d = {}
            for row in results:
                m1, m2, n1, n2, theta, diff = row
                if theta in d.keys():
                    d[theta].append([m1, m2, n1, n2, diff])
                elif theta not in d.keys():
                    d[theta] = [[m1, m2, n1, n2, diff]]
            return d
        else:
            logging.error("Found no matching lattice points.")
            return None

    def find_independent_pairs(self):
        r""" Finds linearly independent pairs of :math:`(\vec{m}_1, \vec{m}_2)` and :math:`(\vec{n}_1, \vec{n}_2)`.

        Constructs supercell matrices :math:`M=(\vec{m}_1^T, \vec{m}_2^T)` and :math:`N=(\vec{n}_1^T, \vec{n}_2^T)` with positive determinants.
        For the same angle, matrices are filtered such that the area of the unit cell :math:`(\vec{m}_1 \times \vec{m}_2)` is minimized and that the angle :math:`(\vec{m}_1 \angle \vec{m}_2)` is relatively close to 90°.

        
        Returns:
            dict: Pairs of angle : (M, N).

        """

        def _area(M):
            area = np.linalg.norm(np.cross(M[:, 0], M[:, 1]))
            return area

        def _angle(M):
            angle = np.arccos(
                np.dot(
                    M[:, 0] / np.linalg.norm(M[:, 0]), M[:, 1] / np.linalg.norm(M[:, 1])
                )
            )
            angle = angle / np.pi * 180
            return angle

        results = self.results
        d = {}
        for ang, indices in results.items():
            for i, pair1 in enumerate(indices):
                pairs = []
                theta = ang * 180.0 / np.pi
                m1, m2, n1, n2, stress = pair1
                for k, pair2 in enumerate(indices):
                    if k > i:
                        m3, m4, n3, n4, _ = pair2
                        M = np.abs(np.cross([m1, m2], [m3, m4])) < 1e-5
                        N = np.abs(np.cross([n1, n2], [n3, n4])) < 1e-5
                        dependent = M and N
                        if not dependent:
                            M = np.array([[m1, m2, 0], [m3, m4, 0], [0, 0, 1]])
                            N = np.array([[n1, n2, 0], [n3, n4, 0], [0, 0, 1]])
                            if (np.linalg.det(M) > 0) and (np.linalg.det(N) > 0):
                                pairs.append((M, N))
                if theta in d.keys():
                    d[theta] += pairs
                elif pairs != []:
                    d[theta] = pairs
        for theta, pairs in d.items():
            bestangle = [abs(_angle(j[0]) - 90) for j in pairs]
            bestangle = [abs(j - min(bestangle)) < 15 for j in bestangle]
            pairs = [j for k, j in enumerate(pairs) if bestangle[k]]
            smallest = min([_area(j[0]) for j in pairs])
            smallest = [abs(_area(j[0]) - smallest) < (1.5 * smallest) for j in pairs]
            pairs = [j for k, j in enumerate(pairs) if smallest[k]]
            d[theta] = pairs
        s = 0
        for k, v in d.items():
            s += len(v)
        logging.info(
            "Constructed {:d} linearly independent pairs of supercell matrices.".format(
                s
            )
        )
        return d

    def analyze_results(
        self, distance=3, weight=0.5, prec=1e-4,
    ):
        """ Analyzes results after grid point search.

        All structures are standardized via spglib to the primitive unit with enforced axis order :math:`a < b < c`.
        Structures are considered identical if they have the same number of atoms, similar unit cell area and stress on the lattice vectors.
        
        Args:
            distance (int, optional): Interlayer distance of reconstructed stacks. Defaults to 3.
            weight (float, optional): Value between 0 and 1, defaults to 0.5. The unit cell of the reconstructed stack is :math:`B + w \cdot (T - B)`.
            prec ([type], optional): Precision to identify equivalent structures. Defaults to 1e-4.
        
        """
        pairs = self.pairs
        self.weight = weight
        assert (self.weight <= 1.0) and (
            self.weight >= 0.0
        ), "Weight must be between 0 and 1."

        data = namedtuple("stack", ["atoms", "M", "N", "angle", "stress"])
        solved = []
        start = time.time()
        logging.info("Reconstructing heterostructures ...")
        for angle, indices in pairs.items():
            for row in indices:
                bottom = self.bottom.atoms.copy()
                top = self.top.atoms.copy()
                M, N = row
                try:
                    bottom = ase.build.make_supercell(bottom, M, wrap=False)
                    top.rotate(angle, v="z", rotate_cell=True)
                    top = ase.build.make_supercell(top, N, wrap=False)
                    bottom = recenter(bottom)
                    top = recenter(top)
                    stack = stack_atoms(
                        bottom, top, distance=distance, weight=self.weight
                    )
                    stress = self.stress_tensor(bottom.cell, top.cell, weight) * 100
                    # stress = np.linalg.norm(top.cell - bottom.cell)
                    solved.append(data(stack, M, N, angle, stress))
                except:
                    logging.error(
                        "ASE supercell generation didn't work for \n{} and \n{}".format(
                            M, N
                        )
                    )
        logging.info(
            "Standardizing representations and filtering unique structures ..."
        )
        self.__filter_unique_structures(solved, prec=prec)
        end = time.time()
        logging.info("Analysis finished in {:.2f} seconds.".format(end - start))

    def stress_tensor(self, cell1, cell2, weight):
        A = cell1.copy()
        B = cell2.copy()
        C = A + weight * (B - A)
        T1 = np.linalg.solve(A, C)
        T1 = T1[[0, 1], :2]
        e1, _ = np.linalg.eig(T1)
        stress1 = np.sum([abs(x - 1) ** 2 for x in e1])
        T2 = np.linalg.solve(B, C)
        T2 = T2[[0, 1], :2]
        e2, _ = np.linalg.eig(T2)
        stress2 = np.sum([abs(x - 1) ** 2 for x in e2])
        return np.sqrt(stress1 + stress2)

    def __spglib_standard(self, solved, prec):
        start = time.time()
        solved = [_standardize(j, prec) for j in solved]
        end = time.time()

        logging.info(
            "  Spglib standardization with enforced axes order finished after {:.2f} seconds ...".format(
                end - start
            )
        )
        return solved

    def __remove_doubles(self, solved):
        def _area(atoms):
            cell = atoms.cell
            area = np.linalg.norm(np.cross(cell[:, 0], cell[:, 1]))
            return area

        def _gamma(atoms):
            a = atoms.cell[:, 0]
            b = atoms.cell[:, 1]
            gamma = np.arccos(np.dot(a / np.linalg.norm(a), b / np.linalg.norm(b)))
            gamma *= 180 / np.pi
            return gamma

        start = time.time()
        data = np.around(
            [[_area(j.atoms), len(j.atoms), j.stress] for j in solved], decimals=3
        )
        u, indices = np.unique(data, return_index=True, axis=0)
        solved = [i for j, i in enumerate(solved) if j in indices]
        end = time.time()
        logging.info(
            "  Filtering structures finished in {:.2f} seconds ...".format(end - start)
        )
        return solved

    def __filter_unique_structures(self, solved, prec=1e-4):
        solved = self.__spglib_standard(solved, prec)
        solved = self.__remove_doubles(solved)
        logging.info("Found {:d} unique structures.".format(len(solved)))
        self.solved = solved
