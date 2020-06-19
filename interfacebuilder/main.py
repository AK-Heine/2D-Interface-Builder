from interfacebuilder.misc import *
from interfacebuilder.structure import structure, recenter, stack_atoms
from interfacebuilder.plotting import interactive_plot
from interfacebuilder.backend import backend_routine
from interfacebuilder.legacy_python_backend import python_backend

import spglib


class interface(interactive_plot, python_backend):
    r"""Interface class to construct heterojunctions from two layered structures.

    Employs the algorithm outlined by Schwalbe-Koda (J. Phys. Chem. C 2016, 120, 20, 10895-10908) to find coincidence lattices of 2D crystals.

    Note:
        The scaling is :math:`\approx (2 \cdot N_{trans} + 1)^4 \cdot N_{angles}`, making the C++ backend preferable. The python backend is deprecated.

    Example:
        >>> from interfacebuilder.main import interface
        >>> intf = interface("bottom.xyz", "top.xyz")
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
        mingamma (float): Minimum value of the gamma angle to filter out structures, defaults to 20.
        maxgamma (float): Maximum value of the gamma angle to filter out structures, defaults to 160.
        backend (str): Backend to perform calculations, either python or c++.
    
    Attributes:
        pairs (dict): Dictionary of angles in degrees and tuples of supercell matrices (M, N).
        solved (list): List of stacks as named tuples after analyzing the results. Contains (atoms, M, N, angle, stress_tot, stress_A, stress_B).

    """

    def __init__(self, bottom, top, **kwargs):
        super().__init__()
        a1, a2, b1, b2 = self.__checks_and_setup(bottom, top)
        self.N_translations = kwargs.get("N_translations", 10)
        self.angle_stepsize = kwargs.get("angle_stepsize", 10.0)
        self.angle_limits = kwargs.get("angle_limits", (0, 180))
        self.tolerance = kwargs.get("tolerance", 0.1)
        self.distance = kwargs.get("distance", 4.0)
        self.weight = kwargs.get("weight", 0.5)
        self.backend = kwargs.get("backend", "c++")
        self.mingamma = kwargs.get("mingamma", 20)
        self.maxgamma = kwargs.get("maxgamma", 160)
        angles = list(
            np.arange(self.angle_limits[0], self.angle_limits[1], self.angle_stepsize)
            / 180
            * np.pi
        )

        assert self.backend in ["c++", "python"], "Backend not accepted."
        if self.backend == "c++":
            logging.debug("Choosing C++ backend.")
            self.pairs = self.__cxx_backend(
                a1,
                a2,
                b1,
                b2,
                angles,
                N=self.N_translations,
                tolerance=self.tolerance,
                mingamma=self.mingamma,
                maxgamma=self.maxgamma,
            )

        if self.backend == "python":
            logging.debug("Choosing old python backend.")
            self.pairs = self.__py_backend_routine(
                a1, a2, b1, b2, angles, N=self.N_translations, tolerance=self.tolerance,
            )

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
        a1 = self.bottom.atoms.cell.T.copy()[:2, 0]
        a2 = self.bottom.atoms.cell.T.copy()[:2, 1]
        b1 = self.top.atoms.cell.T.copy()[:2, 0]
        b2 = self.top.atoms.cell.T.copy()[:2, 1]
        return (a1, a2, b1, b2)

    def __cxx_backend(self, a1, a2, b1, b2, angles, **kwargs):
        """ Python wrapper of the C++ backend. """
        N = self.N_translations
        tolerance = self.tolerance
        mingamma = self.mingamma
        maxgamma = self.maxgamma
        logging.info(
            "Initializing {} gridpoints ...".format(((2 * N + 1) ** 4) * len(angles))
        )
        pairs = backend_routine(
            a1, a2, b1, b2, angles, N, tolerance, mingamma, maxgamma
        )
        pairs = self.__reduce_pairs(pairs)
        if len(pairs) > 0:
            d = {}
            for row in pairs:
                theta, m1, m2, m3, m4, n1, n2, n3, n4 = row
                theta = round(theta * 180 / np.pi, 4)
                M = np.array([[m1, m2, 0], [m3, m4, 0], [0, 0, 1]])
                N = np.array([[n1, n2, 0], [n3, n4, 0], [0, 0, 1]])
                if theta in d.keys():
                    d[theta] += [(M, N)]
                elif theta not in d.keys():
                    d[theta] = [(M, N)]
            return d

    def __reduce_pairs(self, pairs):
        """ Reduces the number of linearly independent pairs.

        For all given pairs with the same area, a point-based evaluation filters out unwanted pairs.
        The point-based system prefers pairs high symmetry and mostly positive integers.

        Args:
            pairs (list): List of pairs theta, m1, m2, m3, m4, n1, n2, n3, n4.

        Returns:
            list: Reduced pairs.
        """
        from math import gcd
        from functools import reduce

        start = time.time()

        def _area(M):
            area = np.abs(np.linalg.det(M))
            return area

        def _is_symmetric(ms):
            if abs(ms[2] - ms[1]) < 1e-3:
                return 1
            else:
                return 0

        def _is_positive(ms):
            if sum(ms) > 0:
                return 1
            else:
                return 0

        def _is_simple(ms):
            if (ms[0] - ms[3]) < 1e-3:
                return 1
            else:
                return 0

        areagroup = []
        for row in pairs:
            theta, m1, m2, m3, m4, n1, n2, n3, n4 = row
            M = np.array([[m1, m2, 0], [m3, m4, 0], [0, 0, 1]])
            areagroup.append(_area(M))

        areagroup = set(areagroup)

        groups = []
        for ar in areagroup:
            subgroup = []
            for row in pairs:
                theta, m1, m2, m3, m4, n1, n2, n3, n4 = row
                M = np.array([[m1, m2, 0], [m3, m4, 0], [0, 0, 1]])
                if _area(M) == ar:
                    subgroup.append(row)
            groups.append(subgroup)

        newgroup = []
        for k, subgroup in enumerate(groups):
            if len(subgroup) > 1:
                symm1 = [_is_symmetric(f[1:5]) for f in subgroup]
                symm2 = [_is_symmetric(f[5:9]) for f in subgroup]
                pos1 = [_is_positive(f[1:5]) for f in subgroup]
                pos2 = [_is_positive(f[5:9]) for f in subgroup]
                simp1 = [_is_simple(f[1:5]) for f in subgroup]
                simp2 = [_is_simple(f[5:9]) for f in subgroup]
                it = zip(symm1, symm2, pos1, pos2, simp1, simp2)
                points = [sum(item) for item in it]
                order = sorted(zip(range(len(subgroup)), points), key=lambda x: x[1])[
                    ::-1
                ]
                indices = [n[0] for n in order if n[1] == order[0][1]]
                for i in indices:
                    newgroup.append(subgroup[i])
            else:
                newgroup.append(subgroup[0])

        anglegroup = list(set([j[0] for j in newgroup]))
        groups = []
        for ang in anglegroup:
            subgroup = []
            for row in newgroup:
                if row[0] == ang:
                    subgroup.append(row)
            groups.append(subgroup)

        newgroup = []
        for k, subgroup in enumerate(groups):
            if len(subgroup) > 1:
                areas = [f[1] * f[4] - f[2] * f[3] for f in subgroup]
                order = sorted(zip(range(len(subgroup)), areas), key=lambda x: x[1])
                indices = [n[0] for n in order if n[1] == order[0][1]]
                for i in indices:
                    newgroup.append(subgroup[i])
                # newgroup.append(subgroup[indices[0]])
            else:
                newgroup.append(subgroup[0])
        end = time.time()
        logging.info("Reduced to {} pairs.".format(len(newgroup)))
        return newgroup

    def analyze_results(self, distance=3, weight=0.5, prec=1e-4, optimize=True):
        r""" Analyzes results after grid point search.
        
        Given the target cell :math:`C = A + w \cdot (B - A)`, one can define a transformation :math:`T_A = C A^{-1}` and :math:`T_B = C B^{-1}`.
        These two-dimensional transformation matrices can be factorized by polar decomposition into a unitary rotation matrix :math:`U` and a strain matrix :math:`P` by :math:`T = U(\phi)P`.
        From :math:`P` one can define a stress matrix as :math:`\epsilon = P - \mathbb{I}` and a stress measure as 
        :math:`\bar{\epsilon} = 0.5 \cdot \sqrt{\epsilon_{xx}^2 + \epsilon_{yy}^2 + \epsilon_{xx} \cdot \epsilon_{yy} + \epsilon_{xy}^2}` .
        The sum of stress measures :math:`\bar{\epsilon}_A + \bar{\epsilon}_ B` can be minimized with respect to :math:`\phi` , yielding a new rotation matrix and a new strain matrix :math:`P' = U^{-1}(\phi_A) C A^{-1} `.
        The returned angles :math:`\phi_A` and :math:`\phi_B` can be used to rotate the individual cells before matching them to the lowest strain common unit cell, resulting in a modified 
        twist angle :math:`\theta' = \theta + \phi_B - \phi_A`.
        
        Args:
            distance (int, optional): Interlayer distance of reconstructed stacks. Defaults to 3.
            weight (float, optional): Value between 0 and 1, defaults to 0.5. The unit cell of the reconstructed stack is :math:`A + w \cdot (B - A)`.
            optimize (bool, optional): Turn on optimization of individual cell rotation to minimize strain. Defaults to true. 
        
        """
        pairs = self.pairs
        self.weight = weight
        assert (self.weight <= 1.0) and (
            self.weight >= 0.0
        ), "Weight must be between 0 and 1."

        data = namedtuple(
            "stack", ["atoms", "M", "N", "angle", "stress", "stress_A", "stress_B"]
        )
        solved = []
        start = time.time()
        logging.info("Building heterostructures ...")
        for angle, indices in pairs.items():
            for row in indices:
                bottom = self.bottom.atoms.copy()
                top = self.top.atoms.copy()
                M, N = row
                try:
                    bottom = ase.build.make_supercell(bottom, M, wrap=False)
                    top.rotate(angle, v="z", rotate_cell=True)
                    top = ase.build.make_supercell(top, N, wrap=False)
                    if optimize:
                        try:
                            newangle1, newangle2 = self.__optimize_cell_rotation(
                                bottom.cell, top.cell, weight
                            )
                            newangle1 *= 180 / np.pi
                            newangle2 *= 180 / np.pi
                            twist = angle + newangle2 - newangle1
                            bottom.rotate(newangle1, v="z", rotate_cell=True)
                            top.rotate(newangle2, v="z", rotate_cell=True)
                        except Exception as e:
                            logging.error(e)
                    bottom = recenter(bottom)
                    top = recenter(top)
                    stack = stack_atoms(
                        bottom, top, distance=distance, weight=self.weight
                    )
                    stress, stress_A, stress_B = self.__measure_stress(
                        bottom.cell, top.cell, weight
                    )
                    stress *= 100
                    twist = angle
                    solved.append(data(stack, M, N, twist, stress, stress_A, stress_B))
                except Exception as e:
                    logging.error(e)
        # self.__filter_unique_structures(solved, prec=prec)
        end = time.time()
        self.solved = solved
        logging.debug("Analysis finished in {:.2f} seconds.".format(end - start))
        logging.info("Analysis finished.")

    def __optimize_cell_rotation(self, cell1, cell2, weight):
        r""" Optimizes indivual cell rotation before matching the common unit cell.
       
        Args:
            cell1 (ndarray): Lower basis.
            cell2 (ndarray): Upper basis.
            weight (float): Weight factor for common unit cell.

        Returns:
            tuple: (:math:`\phi_A`, :math:`\phi_B`)
        """
        from scipy.linalg import polar
        from scipy.optimize import minimize

        def measure(P):
            eps = P - np.identity(2)
            meps = np.sqrt(
                (
                    eps[0, 0] ** 2
                    + eps[1, 1] ** 2
                    + eps[0, 0] * eps[1, 1]
                    + eps[1, 0] ** 2
                )
                / 4
            )
            return meps

        def find(angle, start, target):
            c = np.cos(angle)
            s = np.sin(angle)
            R = np.array([[c, -s], [s, c]])
            newP = np.linalg.inv(R) @ target @ np.linalg.inv(start)
            return measure(newP)

        def f(params, args):
            angle1, angle2 = params
            A, B, C = args
            m1 = find(angle1, A, C)
            m2 = find(angle2, B, C)
            return m1 + m2

        A = cell1.copy()[:2, :2]
        B = cell2.copy()[:2, :2]
        C = A + weight * (B - A)
        T1 = C @ np.linalg.inv(A)
        T2 = C @ np.linalg.inv(B)
        U1, P1 = polar(T1)  # this one goes counterclockwise
        U2, P2 = polar(T2)  # this one goes clockwise
        angle1 = np.arccos(U1[0, 0]) if U1[0, 0] < 0 else -np.arccos(U1[0, 0])
        angle2 = np.arccos(U2[0, 0]) if U2[0, 0] < 0 else -np.arccos(U2[0, 0])
        try:
            res = minimize(f, x0=[angle1, angle2], args=[A, B, C])
            newangle1, newangle2 = res.x
            return (newangle1, newangle2)
        except Exception as e:
            print(e)

    def __measure_stress(self, cell1, cell2, weight):
        """ Measures stress from the polar decomposition of the transformation :math:`UPA = C`.
        
        Args:
            cell1 (ndarray): Lower basis.
            cell2 (ndarray): Upper basis.
            weight (float): Weight factor for common unit cell.

        Returns:
            float: :math:`\bar{\epsilon}_A + \bar{\epsilon}_ B`.
        """
        from scipy.linalg import polar

        A = cell1.copy()[:2, :2]
        B = cell2.copy()[:2, :2]
        C = A + weight * (B - A)
        T1 = C @ np.linalg.inv(A)
        T2 = C @ np.linalg.inv(B)

        def measure(P):
            eps = P - np.identity(2)
            meps = np.sqrt(
                (
                    eps[0, 0] ** 2
                    + eps[1, 1] ** 2
                    + eps[0, 0] * eps[1, 1]
                    + eps[1, 0] ** 2
                )
                / 4
            )
            return meps

        U1, P1 = polar(T1)  # this one goes counterclockwise
        U2, P2 = polar(T2)  # this one goes clockwise
        # u is rotation, p is strain
        meps1 = measure(P1)
        meps2 = measure(P2)
        stress = meps1 + meps2
        return (stress, P1 - np.identity(2), P2 - np.identity(2))
