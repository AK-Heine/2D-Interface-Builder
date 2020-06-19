from interfacebuilder.misc import *
import spglib


class interactive_plot:
    def plot_results(self, jitter=0.05):
        """ Plots results interactively.

        Generates a matplotlib interface that allows to select the reconstructed stacks and save them to a file.
       
        Args:
            jitter (float, optional): Jitters data points to make picking easier. Defaults to 0.05.       

        """
        from matplotlib.widgets import Button
        import matplotlib.cm as cm
        from matplotlib.colors import ListedColormap

        def rand_jitter(arr, jitter):
            stdev = jitter * (max(arr) - min(arr)) + 0.01
            return arr + np.random.randn(len(arr)) * stdev

        data = np.array([[i.stress, len(i.atoms)] for i in self.solved], dtype=float)
        color = [i.angle for i in self.solved]
        norm = matplotlib.colors.Normalize(
            vmin=self.angle_limits[0], vmax=self.angle_limits[1], clip=True
        )
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "", ["darkgreen", "lightgreen", "lightblue", "royalblue"]
        )
        mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        color = [mapper.to_rgba(v) for v in color]

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=([6.4 * 3, 6.4]))
        ax[0].scatter(
            rand_jitter(data[:, 0], jitter),
            data[:, 1],
            color=color,
            alpha=0.75,
            picker=3.5,
            marker=".",
        )
        clb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax[0])
        clb.set_label(r"$\theta$ in °", rotation=270, labelpad=8)
        ax[0].set_ylim(np.min(data[:, 0]) - 0.01, np.max(data[:, 1]) + 0.01)
        ax[0].set_ylim(0, ax[0].get_ylim()[1] + 10)
        ax[0].set_xlabel(r"$\bar{\varepsilon}_A + \bar{\varepsilon}_B$ in %")
        ax[0].set_ylabel("Number of atoms")
        ax[0].set_title("Click a point to select structure.")
        ax[0].grid(
            axis="both", color="lightgray", linestyle="-", linewidth=1, alpha=0.2
        )
        ax[1].set_yticks([])
        ax[1].set_xticks([])
        ax[1].set_xlabel("")
        ax[1].set_ylabel("")
        ax[1].set_frame_on(False)
        ax[2].set_yticks([])
        ax[2].set_xticks([])
        ax[2].set_xlabel("")
        ax[2].set_ylabel("")
        ax[2].set_frame_on(False)
        axbutton = plt.axes([0.8, 0.05, 0.1, 0.05])
        axbutton2 = plt.axes([0.65, 0.05, 0.1, 0.05])
        fig.canvas.mpl_connect("pick_event", self.__onpick)

        def __save(stack):
            try:
                name = "{}_M{}{}{}{}_N{}{}{}{}_a{:.2f}.xyz".format(
                    self.current_stack.get_chemical_formula(), *self.current_scdata[1:]
                )
                stack.write(name, vec_cell=True)
                logging.info("Saved structure to {}".format(name))
            except:
                logging.error("You need to select a point first.")

        def __standardize(atoms):
            atoms = atoms.copy()
            cell = (atoms.get_cell()).tolist()
            pos = atoms.get_scaled_positions().tolist()
            numbers = atoms.get_atomic_numbers()
            cell, scaled_pos, numbers = spglib.standardize_cell(
                (cell, pos, numbers),
                to_primitive=True,
                symprec=1e-4,
                no_idealize=False,
            )
            atoms = ase.atoms.Atoms(
                scaled_positions=scaled_pos, numbers=numbers, cell=cell, pbc=True
            )
            axes = [0, 1, 2]
            lengths = atoms.cell.lengths()
            order = [x for x, y in sorted(zip(axes, lengths), key=lambda pair: pair[1])]
            if order != [0, 1, 2]:
                atoms = ase.geometry.permute_axes(atoms, order)
            self.current_stack = atoms
            self.__plot_stack(atoms, fig, ax[2], self.current_scdata)

        save = Button(axbutton, " Save this structure. ")
        save.on_clicked(lambda x: __save(self.current_stack))
        standard = Button(axbutton2, "spglib standardize")
        standard.on_clicked(lambda x: __standardize(self.current_stack))
        plt.show()

    def __onpick(self, event):
        point = event.artist
        mouseevent = event.mouseevent
        index = event.ind[0]
        fig = point.properties()["figure"]
        axes = fig.axes
        stack = self.solved[index].atoms.copy()
        M = self.solved[index].M
        N = self.solved[index].N
        angle = self.solved[index].angle
        m1, m2, m3, m4 = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
        n1, n2, n3, n4 = N[0, 0], N[0, 1], N[1, 0], N[1, 1]
        scdata = (
            int(len(stack)),
            int(m1),
            int(m2),
            int(m3),
            int(m4),
            int(n1),
            int(n2),
            int(n3),
            int(n4),
            float(angle),
            index,
        )
        self.current_scdata = scdata
        self.current_stack = stack
        self.__plot_stack(stack, fig, axes[2], scdata)
        basis1 = self.bottom.atoms.copy()
        basis2 = self.top.atoms.copy()
        basis2.rotate(angle, v="z", rotate_cell=True)
        self.__plot_lattice_points(
            fig, axes[1], basis1.cell, basis2.cell, scdata,
        )

    def __plot_stack(self, stack, fig, axes, scdata):
        from ase.visualize.plot import plot_atoms

        canvas = fig.canvas
        axes.clear()
        axes.set_yticks([])
        axes.set_xticks([])
        axes.set_xlabel("")
        axes.set_ylabel("")
        scdata = "#{:d}, {:d} atoms, twist angle of {:.2f}°".format(
            scdata[-1], len(stack), scdata[-2]
        )
        axes.set_title(scdata)
        plot_atoms(stack, axes, radii=0.3)
        axes.set_frame_on(False)
        canvas.draw()

    def __plot_lattice_points(self, fig, axes, basis1, basis2, scdata):
        import itertools
        from matplotlib import path, patches

        canvas = fig.canvas
        axes.clear()
        axes.set_yticks([])
        axes.set_xticks([])
        axes.set_xlabel("")
        axes.set_ylabel("")
        natoms, m1, m2, m3, m4, n1, n2, n3, n4, angle, _ = scdata
        sc1 = np.array([[m1, m2], [m3, m4]])
        sc2 = np.array([[n1, n2], [n3, n4]])
        N = self.N_translations + 15

        def plot_grid(fig, axes, basis, sc, color, fc, ec, alpha):
            from matplotlib import path, patches

            basis = basis[:2, :2].copy()
            a1 = basis[0, :].copy()
            a2 = basis[1, :].copy()
            p = itertools.product(range(-N + 1, N), range(-N + 1, N))
            points = np.array([n[0] * a1 + n[1] * a2 for n in p])
            axes.scatter(points[:, 0], points[:, 1], color=color, alpha=alpha, s=3)
            SC = sc @ basis
            path1 = [(0, 0), (SC[0, :]), (SC[0, :] + SC[1, :]), (SC[1, :]), (0, 0)]
            path1 = path.Path(path1)
            patch = patches.PathPatch(
                path1, facecolor=fc, edgecolor=ec, alpha=alpha, lw=2
            )
            axes.add_patch(patch)

        # first cell
        plot_grid(fig, axes, basis1, sc1, "darkred", "none", "darkred", 0.5)
        # second cell
        plot_grid(fig, axes, basis2, sc2, "darkblue", "none", "darkblue", 0.5)
        # supercell
        C = sc1 @ basis1[:2, :2] + self.weight * (
            sc2 @ basis2[:2, :2] - sc1 @ basis1[:2, :2]
        )
        it = itertools.product(range(-N + 1, N), range(-N + 1, N))
        a1 = C[0, :]
        a2 = C[1, :]
        points = np.array([n[0] * a1 + n[1] * a2 for n in it])
        path3 = [(0, 0), (C[0, :]), (C[0, :] + C[1, :]), (C[1, :]), (0, 0)]
        axes.scatter(
            points[:, 0],
            points[:, 1],
            facecolor="none",
            edgecolor="tab:green",
            s=20,
            linewidth=2,
        )
        p = path.Path(path3)
        patch = patches.PathPatch(
            p, facecolor="tab:purple", edgecolor="none", lw=2, alpha=0.5
        )
        axes.add_patch(patch)
        path3 = np.array(path3)
        xlim = (np.min(path3[:, 0]) - 4, np.max(path3[:, 0] + 4))
        ylim = (np.min(path3[:, 1]) - 4, np.max(path3[:, 1] + 4))
        axes.axis("equal")
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)
        axes.set_frame_on(False)

        scdata = """M = ({: 2d}, {: 2d}, {: 2d}, {: 2d})\nN = ({: 2d}, {: 2d}, {: 2d}, {: 2d})""".format(
            m1, m2, m3, m4, n1, n2, n3, n4
        )
        axes.set_title(scdata)

        canvas.draw()
