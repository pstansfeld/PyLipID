==========================================================
PyLipID - A Python Package For Lipid Interactions Analysis
==========================================================

.. image:: https://github.com/pstansfeld/PyLipID/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/pstansfeld/PyLipID/actions/workflows/ci.yml
.. image:: https://img.shields.io/pypi/v/PyLipID
   :target: https://pypi.org/project/pylipid/
.. image:: https://img.shields.io/pypi/pyversions/PyLipID
   :target: https://pypi.org/project/pylipid/

.. image:: docs/static/pylipid_logo_smallsize.png
    :align: center


PyLipID is a Python package for analysing lipid interactions with membrane proteins from
Molecular Dynamics (MD) simulations. It supports both coarse-grained and all-atom force fields,
works with any trajectory format supported by `MDTraj <https://mdtraj.org>`_, and is designed
to run from Jupyter notebooks or Python scripts.

Full documentation and tutorials: `pylipid.readthedocs.io <https://pylipid.readthedocs.io>`_


Features
========

- **Dual-cutoff contact detection** — overcomes the 'rattling in cage' artefact common in
  coarse-grained simulations by defining contact start and end with separate distance thresholds.
- **Binding site detection** — identifies lipid binding sites by computing community structures
  in residue interaction networks using the Louvain algorithm.
- **Kinetics** — calculates lipid koff rates and residence times for both individual residues
  and binding sites.
- **Interaction metrics** — computes occupancy, lipid count, duration and residence time per
  residue and per binding site.
- **Bound pose analysis** — extracts and clusters representative lipid-bound poses for each
  binding site.
- **Surface area** — calculates solvent-accessible surface area contributions per residue and
  binding site.
- **Visualisation** — generates publication-ready figures and PyMOL scripts for structural
  visualisation.


Installation
============

Install from PyPI (recommended)::

    pip install pylipid

Install from source::

    git clone https://github.com/pstansfeld/PyLipID
    cd PyLipID
    pip install .

For development (editable install)::

    pip install -e .

**Requirements:** Python >= 3.9. All dependencies are installed automatically by pip.


Quick Start
===========

.. code-block:: python

    from pylipid.api import LipidInteraction

    # Point to your trajectory and topology files
    trajfile_list = ["run1/protein_lipids.xtc", "run2/protein_lipids.xtc"]
    topfile_list  = ["run1/protein_lipids.gro", "run2/protein_lipids.gro"]

    # Initialise — dual cutoffs in nm, lipid residue name as in your topology
    li = LipidInteraction(
        trajfile_list,
        topfile_list=topfile_list,
        cutoffs=[0.55, 0.8],   # [lower, upper] in nm
        lipid="CHOL",
        nprot=1,
    )

    # Collect contacts (the main calculation step)
    li.collect_residue_contacts()

    # Compute metrics
    li.compute_residue_duration()
    li.compute_residue_koff()
    li.compute_binding_nodes(threshold=4)
    li.compute_site_koff()

    # Plot and save
    li.plot(item="Residence Time")
    li.save_data(item="Dataset")
    li.save_pymol_script(pdb_file="receptor.pdb")

For a full walkthrough including coarse-grained and all-atom examples, see the
`tutorials on readthedocs <https://pylipid.readthedocs.io>`_.


What's New in v1.6.0
====================

This release, maintained by the `Stansfeld Lab <https://github.com/pstansfeld>`_,
fixes a number of bugs present in v1.5.14 and modernises the package for current
Python and dependency versions.

**Bug fixes**

- Fixed off-by-one error in binding site detection that silently dropped the last
  binding site on every run.
- Fixed ``calculate_scores`` returning ``None`` on error, causing a downstream
  ``TypeError`` in pose analysis.
- Fixed correlation coefficient matrix computation (``sparse_corrcoef``) returning a
  ``numpy.matrix`` object incompatible with NetworkX 3.x and newer NumPy (root cause
  of the ``IndexError`` reported in issue #33).
- Fixed ``analyze_bound_poses`` referencing a non-existent column name and using
  incorrect DataFrame indexing.
- Fixed ``compute_residue_lipidcount`` missing its ``return`` statement for single-residue
  queries.
- Fixed ``write_PDB`` crashing with ``ValueError: invalid literal for int()`` when chain
  IDs are absent in GROMACS ``.gro`` topologies (issues #29, #26).
- Fixed hardcoded figure number in ``plot_koff`` causing curves to overwrite each other
  when called for multiple residues.
- Fixed deprecated ``plt.get_cmap()`` call removed in Matplotlib 3.9.
- Fixed ``LogNorm`` crash in ``plot_corrcoef`` when the correlation matrix contains zeros.
- Fixed ``np.float128`` (Linux-only) replaced with portable ``np.longdouble``.
- Fixed typo in output folder name ``Reisidue_koffs`` -> ``Residue_koffs``.

**Performance**

- ``collect_residue_contacts`` now issues a single ``md.compute_distances`` call per
  residue rather than one per lipid molecule, giving a significant speedup on systems
  with large lipid counts. A nested progress bar now shows per-residue progress within
  each trajectory so the calculation no longer appears frozen.

**Dependency changes**

- Removed ``python-louvain`` dependency; binding site detection now uses the Louvain
  implementation built into NetworkX (>= 2.6). Results are now reproducible across
  runs (fixed random seed).
- All dependencies now carry minimum version bounds.
- Python < 3.9 (end-of-life) is no longer supported.

**Packaging**

- Migrated from ``setup.py`` to ``pyproject.toml``.
- CI migrated from Travis to GitHub Actions (Python 3.9--3.12, Ubuntu and macOS).


Citation
========

If you use PyLipID in scientific research, please cite:

.. image:: https://img.shields.io/badge/DOI-10.1021/acs.jctc.1c00708-blue
   :target: https://doi.org/10.1021/acs.jctc.1c00708

::

    @article{song_pylipid_2022,
        author  = {Song, Wanling and Corey, Robin A. and Ansell, T. Bertie and
                   Cassidy, C. Keith and Horrell, Michael R. and Duncan, Anna L.
                   and Stansfeld, Phillip J. and Sansom, Mark S.P.},
        title   = {PyLipID: A Python package for analysis of protein-lipid
                   interactions from MD simulations},
        journal = {J. Chem. Theory Comput.},
        year    = {2022},
        doi     = {10.1021/acs.jctc.1c00708},
        url     = {https://doi.org/10.1021/acs.jctc.1c00708},
    }


Licence
=======

MIT licence. See `LICENSE.txt <LICENSE.txt>`_ for details.
