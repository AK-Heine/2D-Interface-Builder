from interfacebuilder.misc import *


class interactive_plot:
    def __init__(self):
        pass

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
            stdev = jitter * (max(arr) - min(arr))
            return arr + np.random.randn(len(arr)) * stdev

        data = np.array([[i.stress, len(i.atoms)] for i in self.solved], dtype=float)
        color = [i.angle for i in self.solved]
        norm = matplotlib.colors.Normalize(vmin=min(color), vmax=max(color), clip=True)
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "", ["darkgreen", "lightgreen", "lightblue", "royalblue"]
        )
        mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        color = [mapper.to_rgba(v) for v in color]

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=([6.4 * 3, 4.8 * 1.5]))
        ax[0].scatter(
            rand_jitter(data[:, 0], jitter),
            data[:, 1],
            color=color,
            alpha=0.75,
            picker=3.5,
            marker=".",
        )
        clb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax[0])
        clb.ax.set_title("angle in °")
        ax[0].set_ylim(np.min(data[:, 0]) - 0.01, np.max(data[:, 1]) + 0.01)
        ax[0].set_ylim(0, ax[0].get_ylim()[1] + 10)
        ax[0].set_xlabel("stress measure")
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

        save = Button(axbutton, " Save this structure. ")
        save.on_clicked(lambda x: __save(self.current_stack))
        plt.show()

    def __onpick(self, event):
        point = event.artist
        mouseevent = event.mouseevent
        index = event.ind[0]
        fig = point.properties()["figure"]
        axes = fig.axes
        stack = self.solved[index].atoms
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
        )
        self.current_scdata = scdata
        self.__plot_stack(stack, fig, axes[2], scdata)
        basis1 = self.bottom.atoms.copy()
        basis2 = self.top.atoms.copy()
        basis2.rotate(scdata[-1], v="z", rotate_cell=True)
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
        scdata = "{:d} atoms, twist angle of {:.2f}°".format(len(stack), scdata[-1])
        axes.set_title(scdata)
        self.current_stack = stack
        plot_atoms(stack, axes, radii=0.3)
        axes.set_frame_on(False)
        canvas.draw()

    def __plot_lattice_points(self, fig, axes, basis1, basis2, scdata, N=20):
        import itertools
        from matplotlib import path, patches

        canvas = fig.canvas
        axes.clear()
        axes.set_yticks([])
        axes.set_xticks([])
        axes.set_xlabel("")
        axes.set_ylabel("")

        sc1 = np.array(scdata[1:5]).reshape((2, 2))
        sc2 = np.array(scdata[5:9]).reshape((2, 2))

        # first cell
        A = basis1[[0, 1], :2]
        a1 = A[:, 0]
        a2 = A[:, 1]
        p = itertools.product(range(-N + 1, N), repeat=2)
        new = list()
        for n in p:
            point = n[0] * a1 + n[1] * a2
            new.append(point)
        Agrid = np.stack(new)
        axes.scatter(Agrid[:, 0], Agrid[:, 1], color="darkred", alpha=0.5, s=3)
        A = sc1 @ A
        path1 = [(0, 0), (A[:, 0]), (A[:, 0] + A[:, 1]), (A[:, 1]), (0, 0)]
        p = path.Path(path1)
        patch = patches.PathPatch(p, facecolor="red", lw=1, alpha=0.15)
        axes.add_patch(patch)

        # second cell
        B = basis2[[0, 1], :2]
        b1 = B[:, 0]
        b2 = B[:, 1]
        p = itertools.product(range(-N + 1, N), repeat=2)
        new = list()
        for n in p:
            point = n[0] * b1 + n[1] * b2
            new.append(point)
        Bgrid = np.stack(new)
        axes.scatter(Bgrid[:, 0], Bgrid[:, 1], color="darkblue", alpha=0.5, s=3)
        B = sc2 @ B
        path2 = [(0, 0), (B[:, 0]), (B[:, 0] + B[:, 1]), (B[:, 1]), (0, 0)]
        p = path.Path(path2)
        patch = patches.PathPatch(p, facecolor="darkblue", lw=1, alpha=0.15)
        axes.add_patch(patch)

        # supercell
        C = A + self.weight * (B - A)
        path3 = [(0, 0), (C[:, 0]), (C[:, 0] + C[:, 1]), (C[:, 1]), (0, 0)]
        p = path.Path(path3)
        patch = patches.PathPatch(
            p, facecolor="none", edgecolor="purple", lw=2, alpha=0.5
        )
        axes.add_patch(patch)
        path3 = np.array(path3)
        xlim = (np.min(path3[:, 0]) - 4, np.max(path3[:, 0] + 4))
        ylim = (np.min(path3[:, 1]) - 4, np.max(path3[:, 1] + 4))
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)
        axes.set_frame_on(False)
        natoms, m1, m2, m3, m4, n1, n2, n3, n4, angle = scdata
        scdata = """M = ({: 2d}, {: 2d}, {: 2d}, {: 2d})\nN = ({: 2d}, {: 2d}, {: 2d}, {: 2d})""".format(
            m1, m2, m3, m4, n1, n2, n3, n4
        )
        axes.set_title(scdata)
        canvas.draw()
