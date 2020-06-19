from interfacebuilder.misc import *
from interfacebuilder.structure import structure, recenter, stack_atoms
from interfacebuilder.plotting import interactive_plot

import spglib


def _coincidence(a1, b1, m1, m2, n1, n2, theta):
    R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    Am = a1 @ np.array([m1, m2])
    Bn = R @ b1 @ np.array([n1, n2])
    return (m1, m2, n1, n2, theta, np.linalg.norm(Am - Bn))


class python_backend:
    def py_find_coincidence(
        self,
        a_basis,
        b_basis,
        N=10,
        angles=[0],
        tolerance=0.05,
        chunksize=10000,
        **kwargs
    ):
        """ Solves the system of linear equations given by :math:`\mathbf{A m = R B n}`.
    
        Note:
            The python backend is deprecated and only remains for testing purposes.

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
        import multiprocessing as mp

        chunksize = kwargs.get("chuncksize", 10000)
        start = time.time()
        cpus = mp.cpu_count()
        pool = mp.Pool(processes=cpus)
        nrange = range(-N, N + 1)
        iterator = (
            (a_basis, b_basis, m1, m2, n1, n2, theta)
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
        return results

    def py_find_independent_pairs(self, results):
        r""" Finds linearly independent pairs of :math:`(\vec{m}_1, \vec{m}_2)` and :math:`(\vec{n}_1, \vec{n}_2)`.

        Constructs supercell matrices :math:`M=(\vec{m}_1^T, \vec{m}_2^T)` and :math:`N=(\vec{n}_1^T, \vec{n}_2^T)` with positive determinants.
        For the same angle, matrices are filtered such that the area of the unit cell :math:`(\vec{m}_1 \times \vec{m}_2)` is minimized and that the angle :math:`(\vec{m}_1 \angle \vec{m}_2)` is relatively close to 90°.

        
        Returns:
            dict: Pairs of angle : (M, N).

        """

        def _area(M):
            area = np.linalg.norm(np.cross(M[:, 0], M[:, 1]))
            return area

        d = {}
        for ang, indices in results.items():
            for i, pair1 in enumerate(indices):
                pairs = []
                theta = ang * 180.0 / np.pi
                m1, m2, n1, n2, stress = pair1
                for k, pair2 in enumerate(indices):
                    if k != i:
                        m3, m4, n3, n4, _ = pair2
                        M = np.array([[m1, m2, 0], [m3, m4, 0], [0, 0, 1]])
                        N = np.array([[n1, n2, 0], [n3, n4, 0], [0, 0, 1]])
                        if (np.linalg.det(M) > 0) and (np.linalg.det(N) > 0):
                            pairs.append((M, N))
                if theta in d.keys():
                    d[theta] += pairs
                elif pairs != []:
                    d[theta] = pairs
        # for theta, pairs in d.items():
        #     smallest = min([_area(j[0]) for j in pairs])
        #     smallest = [abs(_area(j[0]) - smallest) < (1.5 * smallest) for j in pairs]
        #     pairs = [j for k, j in enumerate(pairs) if smallest[k]]
        #     d[theta] = pairs
        s = 0
        for k, v in d.items():
            s += len(v)
        logging.info(
            "Constructed {:d} linearly independent pairs of supercell matrices.".format(
                s
            )
        )
        return d

    def py_backend_routine(self, a1, a2, b1, b2, angles, N, tolerance):
        a_basis = np.column_stack((a1, a2))
        b_basis = np.column_stack((b1, b2))
        results = self.py_find_coincidence(
            a_basis, b_basis, N=N, angles=angles, tolerance=tolerance, chuncksize=10000
        )
        pairs = self.py_find_independent_pairs(results)
        return pairs
