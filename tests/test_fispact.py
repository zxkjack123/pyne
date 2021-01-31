"""fispact tests """
import os
import nose 
from nose.tools import assert_equal, assert_true, assert_almost_equal, assert_raises
import numpy as np

import warnings
from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)

from pyne.mesh import HAVE_PYMOAB
from pyne.mesh import Mesh, StatMesh, MeshError
from pyne import fispact
from pyne.material import Material

thisdir = os.path.dirname(__file__)
fispactii_path="fispii.out"

fo=fispact.read_fis_out(fispactii_path)

def test_read_fis_out():
    """test simple single values have been read """
    assert_equal(fo.file_name, fispactii_path)
    assert_equal(fo.version, "FISPACT-II")
    assert_true(fo.isFisII)
    assert_equal(fo.tot_irrad_time, 8.640000E+06)
    assert_equal(fo.tot_fluence, 8.640000E+16)
    assert_equal(fo.ave_flux, 1.000000E+10)
    assert_equal(fo.cpu_time, 1.9417)
    assert_equal(fo.num_irrad_step, 1)

def test_read_time_step():
    """test reading time steps for fispact-II"""
    ts1 = fo.timestep_data[0]
    ts2 = fo.timestep_data[4]
    ts3 = fo.timestep_data[-1]
    assert_equal(len(fo.timestep_data), 11)
    assert_equal(ts1.step_length, 8.6400E+06)
    assert_equal(ts2.step_length, 3.6000E+03)
    assert_equal(ts3.step_length, 3.1558E+10)
    assert_equal(ts1.alpha_act, 4.803006E+02)
    assert_equal(ts1.beta_act, 2.050816E+11)
    assert_equal(ts1.gamma_act, 4.525181E+09)
    assert_equal(ts1.total_act, 2.09607E+11)
    assert_equal(ts1.total_act_no_trit, 2.09607E+11)
    assert_equal(ts1.alpha_heat, 1.38577E-13)
    assert_equal(ts1.beta_heat, 1.77732E-05)
    assert_equal(ts1.gamma_heat, 3.65898E-05)
    assert_equal(ts1.total_heat, 5.43630E-05)
    assert_equal(ts1.total_heat_no_trit, 5.43630E-05)
    assert_equal(ts1.num_nuclides, 170)
    assert_equal(ts1.num_fissions, "-2.01548-191")
    assert_equal(ts1.neutron_flux, 1.00000E+10)
    assert_equal(ts1.initial_mass, 1.0)
    assert_equal(ts1.total_mass, 1.0)
    assert_equal(ts1.density, 7.93)
    assert_equal(ts1.actinide_burn, 0)

    assert_equal(ts3.appm_h1, 6.1126E-03)
    assert_equal(ts3.appm_h2, 4.4613E-10)
    assert_equal(ts3.appm_h3, 1.8511E-20)
    assert_equal(ts3.appm_he3, 1.2419E-08)
    assert_equal(ts3.appm_he4, 6.8193E-03)


def test_read_spectra():
    """test read of spectra data each time step for fispact-II """
    ts1 = fo.timestep_data[0]
    ts2 = fo.timestep_data[5]
    ts3 = fo.timestep_data[-1]
    assert_equal(ts3.gspec[0], 1.88660E+03)
    assert_equal(ts3.gspec[-1], 0.00000E+00)
    assert_equal(len(ts1.gspec), 24)
    assert_equal(len(ts2.gspec), 24)
    assert_equal(len(ts3.gspec), 24)


def test_read_dominant():
    """test read of dominant data each time step for fispact-II"""
    ts1 = fo.timestep_data[0]

    assert_equal(len(ts1.dom_data[0]), 96)
    assert_equal(len(ts1.dom_data[0]), len(ts1.dom_data[1]))
    assert_equal(len(ts1.dom_data[1]), len(ts1.dom_data[2]))

    assert_equal(ts1.dom_data[0][0], "Mn56")
    assert_equal(float(ts1.dom_data[1][0]), 1.2883E+11)
    assert_equal(float(ts1.dom_data[2][0]), 61.46E+00)
    assert_equal(ts1.dom_data[3][0], "Mn56")
    assert_equal(float(ts1.dom_data[4][0]), 5.2249E-05)
    assert_equal(float(ts1.dom_data[5][0]), 96.11E+00)
    assert_equal(ts1.dom_data[6][0], "Mn56")
    assert_equal(float(ts1.dom_data[7][0]), 2.9243E-04)
    assert_equal(float(ts1.dom_data[8][0]), 84.68E+00)
    assert_equal(ts1.dom_data[9][0], "Mn56")
    assert_equal(float(ts1.dom_data[10][0]), 3.5299E-05)
    assert_equal(float(ts1.dom_data[11][0]), 96.47E+00)
    assert_equal(ts1.dom_data[12][0], "Mn56")
    assert_equal(float(ts1.dom_data[13][0]), 1.6949E-05)
    assert_equal(float(ts1.dom_data[14][0]), 95.36E+00)


def test_read_composition():
    """test read of composition data each time step for fispact-II"""
    ts1 = fo.timestep_data[0]
    ts3 = fo.timestep_data[-1]

    assert_equal(len(ts1.composition[0]), 36)
    assert_equal(len(ts3.composition[0]), 36)

def test_read_inv():
    """test read of inventory data for each time step for fispact-II """
    ts1 = fo.timestep_data[0]
    ts3 = fo.timestep_data[-1]

    assert_equal(len(ts1.inventory), ts1.num_nuclides)
    assert_equal(len(ts3.inventory), ts3.num_nuclides)
    assert_equal(ts1.inventory[0,0], "H1")
    assert_equal(float(ts1.inventory[0,1]), 6.67568E+16)
    assert_equal(float(ts1.inventory[0,2]), 1.117E-07)
    assert_equal(float(ts1.inventory[0,3]), 0)
    assert_equal(float(ts1.inventory[0,4]), 0)
    assert_equal(float(ts1.inventory[0,5]), 0)
    assert_equal(float(ts1.inventory[0,6]), 0)
    assert_equal(float(ts1.inventory[0,7]), 0)

    assert_equal(ts1.inventory[-1,0], "Os189")
    assert_equal(float(ts1.inventory[-1,1]), 2.53243E+04)
    assert_equal(float(ts1.inventory[-1,2]), 7.946E-18 )
    assert_equal(float(ts1.inventory[-1,3]), 0)
    assert_equal(float(ts1.inventory[-1,4]), 0)
    assert_equal(float(ts1.inventory[-1,5]), 0)
    assert_equal(float(ts1.inventory[-1,6]), 0)
    assert_equal(float(ts1.inventory[-1,7]), 0)
 
def test_write_fluxin_single():
    """This function tests the mesh_to_fispact_fluxin function for a single energy
    group case.
    """

    if not HAVE_PYMOAB:
        raise SkipTest

    output_dir = os.path.join(thisdir, "files_test_fispact")
    forward_fluxin = os.path.join(thisdir, "files_test_fispact",
                                  "fluxin_single_forward.txt")
    flux_mesh = Mesh(structured=True,
                     structured_coords=[[0, 1, 2], [0, 1, 2], [0, 1]])
    tag_flux = flux_mesh.tag(name="flux", size=1, dtype=float)
    flux_data = [1, 2, 3, 4]
    ves = flux_mesh.structured_iterate_hex("xyz")
    for i, ve in enumerate(ves):
        flux_mesh.flux[i] = flux_data[i]

    # test forward writting
    fispact.mesh_to_fispact_fluxin(flux_mesh, flux_tag="flux",
                                   fispact_files_dir=output_dir, reverse=False)

    # test flux file of the first ve
    output = os.path.join(output_dir, ''.join(["ve0.flx"]))
    with open(output) as f:
        written = f.readlines()
    with open(forward_fluxin) as f:
        expected = f.readlines()
    assert_equal(written, expected)

    # remove generated output
    for i in range(4):
        output = os.path.join(output_dir, ''.join(["ve", str(i), ".flx"]))
        if os.path.isfile(output):
            os.remove(output)

def test_write_fluxin_multiple():
    """This function tests the mesh_to_fispact_fluxin function for a multiple
    energy group case.
    """

    if not HAVE_PYMOAB:
        raise SkipTest

    output_dir = os.path.join(thisdir, "files_test_fispact")
    forward_fluxin = os.path.join(thisdir, "files_test_fispact",
                                  "fluxin_multiple_forward.txt")
    reverse_fluxin = os.path.join(thisdir, "files_test_fispact",
                                  "fluxin_multiple_reverse.txt")

    flux_mesh = Mesh(structured=True,
                     structured_coords=[[0, 1, 2], [0, 1], [0, 1]])
    flux_data = [[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14]]
    flux_mesh.tag("flux", flux_data, 'nat_mesh', size=7, dtype=float)

    # test forward writting of first ve
    fispact.mesh_to_fispact_fluxin(flux_mesh, "flux", output_dir, False)
    output = os.path.join(output_dir, ''.join(["ve0.flx"]))
    with open(output) as f:
        written = f.readlines()

    with open(forward_fluxin) as f:
        expected = f.readlines()
    assert_equal(written, expected)

    # remove generated output
    for i in range(4):
        output = os.path.join(output_dir, ''.join(["ve", str(i), ".flx"]))
        if os.path.isfile(output):
            os.remove(output)

    # test reverse writting of first ve
    fispact.mesh_to_fispact_fluxin(flux_mesh, "flux", output_dir, True)
    output = os.path.join(output_dir, ''.join(["ve0.flx"]))
    with open(output) as f:
        written = f.readlines()
    with open(reverse_fluxin) as f:
        expected = f.readlines()
    assert_equal(written, expected)

    # remove generated output
    for i in range(4):
        output = os.path.join(output_dir, ''.join(["ve", str(i), ".flx"]))
        if os.path.isfile(output):
            os.remove(output)


def test_write_fispact_input_single_ve():
    """This function tests the write_fispact_input function for test.
    """
    h2o = {10010000: 0.11191487328808077, 80160000: 0.8880851267119192}
    mat = Material(h2o, density=1.0)
    output_dir = os.path.join(thisdir, "files_test_fispact")
    filename = os.path.join(output_dir, "test")
    fispact.write_fispact_input_single_ve(filename=filename, material=mat,
        decay_times=['1 s'], total_flux=1e-2)
    os.remove(filename+".i")

def test_write_fispact_input():
    """This function tests the write_fispact_input function for test.
    """
    flux_mesh = Mesh(structured=True,
                     structured_coords=[[0, 1, 2], [0, 1, 2], [0, 1]])
    flux_data = [[1, 2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14],
                 [15, 16, 17, 18, 19, 20, 21],
                 [22, 23, 24, 25, 26, 27, 28]]
    flux_mesh.tag("flux", flux_data, 'nat_mesh', size=7, dtype=float)
    flux_mesh.tag("n_flux_total", np.sum(flux_data), 'nat_mesh', size=1, dtype=float)

    cell_fracs = np.zeros(6, dtype=[('idx', np.int64),
                                    ('cell', np.int64),
                                    ('vol_frac', np.float64),
                                    ('rel_error', np.float64)])
    cell_fracs[:] = [(0, 11, 1, 0.0), (1, 11, 0.5, 0.0), (1, 12, 0.5, 0.0),
                     (2, 11, 0.5, 0.0), (2, 13, 0.5, 0.0), (3, 13, 1, 0.0)]
    flux_mesh.tag_cell_fracs(cell_fracs)
    cell_mats = {11: Material({'H': 1.0}, density=1.0),
                 12: Material({'He': 1.0}, density=1.0),
                 13: Material({}, density=0.0, metadata={'name': 'void'})}

    # write fispact input files
    fispact_files_dir = os.path.join(thisdir, "files_test_fispact")
    fispact.write_fispact_input(flux_mesh, cell_fracs, cell_mats,
        fispact_files_dir=fispact_files_dir, decay_times=['1 s'])

