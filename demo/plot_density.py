"""
Plot the DTFE density field computed by the demo.

Usage (from the main directory of the repository):
    python demo/plot_density.py                         # uses demo/demo_output.den
    python demo/plot_density.py my_output.den --box 100 --out my_plot.png

Left panel : slab of thickness 2*D cells centred on the middle of the box along z,
             with D = ngrid/16 (can be changed with --halfwidth).
Right panel: projection along the full z axis.
Both panels show log10(rho/<rho>), i.e. the mean density along the line of sight in units of the mean density.

The output of DTFE ('.den' file, binary output type 101) is a raw array of ngrid^3 single-precision floats.
The grid size is deduced from the file size.

In a Jupyter notebook you can also run:  %run demo/plot_density.py
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt


def read_density(filename):
    """Read a DTFE binary density file and return it as an (n, n, n) array."""
    rho = np.fromfile(filename, dtype=np.float32)
    n = round(rho.size ** (1 / 3))
    if n**3 != rho.size:
        raise ValueError(f"'{filename}' has {rho.size} values, which is not a cubic grid.")
    return rho.reshape(n, n, n)


def main():
    parser = argparse.ArgumentParser(description="Plot a slab and a projection of a DTFE density field.")
    parser.add_argument("file", nargs="?", default="demo/demo_output.den", help="DTFE density file [demo/demo_output.den]")
    parser.add_argument("--box", type=float, default=50.0, help="box side length in Mpc [50]")
    parser.add_argument("--halfwidth", type=int, default=None, help="slab half-width D in cells [ngrid/16]")
    parser.add_argument("--out", default="demo/demo_density.png", help="save the figure to this file [demo/demo_density.png]")
    args, _ = parser.parse_known_args()   # 'known' so that it also works inside Jupyter

    rho = read_density(args.file)
    n = rho.shape[0]
    D = args.halfwidth if args.halfwidth is not None else max(n // 16, 1)
    delta = rho / rho.mean()
    print(f"grid: {n}^3   min/max of rho/<rho>: {delta.min():.3g} / {delta.max():.3g}   empty cells: {(rho <= 0).sum()}")

    i0, i1 = n // 2 - D, n // 2 + D
    dz = args.box / n
    slab = delta[:, :, i0:i1].mean(axis=2)
    proj = delta.mean(axis=2)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    panels = [(axes[0], slab, f"Slab {i0 * dz:.1f} < z < {i1 * dz:.1f} Mpc ({2 * D} cells)"),
              (axes[1], proj, "Projection along the full z axis")]
    for ax, img, title in panels:
        im = ax.imshow(np.log10(np.maximum(img, 1e-3)).T, origin="lower",
                       extent=[0, args.box, 0, args.box], cmap="inferno")
        ax.set_title(title)
        ax.set_xlabel("x [Mpc]")
        ax.set_ylabel("y [Mpc]")
        fig.colorbar(im, ax=ax, label=r"$\log_{10}(\rho/\langle\rho\rangle)$")
    fig.tight_layout()

    if args.out:
        fig.savefig(args.out, dpi=150)
        print(f"Figure saved to '{args.out}'")
    plt.show()


if __name__ == "__main__":
    main()
