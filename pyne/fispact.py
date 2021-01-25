"""
Module for parsing FISPACT output data
FISPACT and FISPACT-II are bateman equation solvers for transmutation
and fission product yield calculations.  FISPACT-II is developed and maintained
by the UKAEA.
it supports not only neutron irradiaiton but gamma, proton, deuteron and triton
irradiation. it has support for self shielding and can read endf format nuclear
data

this module has methods for parsing the fispact output file,
extracting data, processing the data
"""
import os
import numpy as np
from pyne.mesh import Mesh, HAVE_PYMOAB
from pyne.material import Material, from_atom_frac
from pyne.nucname import name

try:
    basestring
except NameError:
    basestring = str

try:
    import pypact as pp
    HAVE_PYPACT = True
except ImportError:
    HAVE_PYPACT = False
    warn("The pypact optional dependency could not be imported. "
         "The functions related to FISPACT input file writor may be incomplete.", ImportWarning)

if HAVE_PYMOAB:
    from pyne.mesh import mesh_iterate
else:
    warn("The PyMOAB optional dependency could not be imported. "
         "Some aspects of the mesh module may be incomplete.", ImportWarning)

class FispactOutput():
    """ fispact output data"""
    def __init__(self):
        """ define data"""
        self.file_name = ""
        self.sumdat = []
        self.timestep_data = []
        self.num_cool_step = 0   # number of steps after zero keyword
        self.num_irrad_step = 0  # number of steps with flux > 0
        self.version = ""
        self.isFisII = False
        self.cpu_time = 0.0
        self.tot_irrad_time = 0.0
        self.tot_fluence = 0.0
        self.ave_flux = 0.0
        self.time_days = []


class FispactTimeStep():
    """ data for an individual time step can be heating or cooling """
    def __init__(self):
        self.step_num = 1
        self.step_length = 0
        self.flux_amp = 0
        self.is_cooling = False
        self.num_nuclides = 0
        self.alpha_act = 0     # bq
        self.beta_act = 0      # bq
        self.gamma_act = 0     # bq
        self.total_act = 0     # bq
        self.total_act_no_trit = 0   # bq
        self.alpha_heat = 0    # kw
        self.beta_heat = 0     # kw
        self.gamma_heat = 0    # kw
        self.total_heat = 0    # kw
        self.total_heat_no_trit = 0   # kw
        self.num_fissions = 0
        self.neutron_flux = 0   # n/cm**2/s
        self.initial_mass = 0   # kg
        self.total_mass = 0     # kg
        self.density = 0        # g/cc
        self.actinide_burn = 0   # %
        self.appm_h1 = 0
        self.appm_h2 = 0
        self.appm_h3 = 0
        self.appm_he3 = 0
        self.appm_he4 = 0

        self.dom_data = []
        self.inventory = []
        self.gspec = []
        self.composition = []


def read_fis_out(path):
    """ parse a fispact output file
        returns fo, a fispact output object
    """

    fo = FispactOutput()
    fo.file_name = path
    with open(path) as f:
        lines = f.read().splitlines()

    fo.version = check_fisp_version(lines)
    fo.isFisII = isFisII(lines)
    fo.sumdat = read_summary_data(lines)

    fo.ave_flux = read_parameter(lines, "Mean flux")
    fo.tot_irrad_time = read_parameter(lines, "Total irradiation time")
    fo.tot_fluence = read_parameter(lines, "Total fluence")
    fo.num_irrad_step = read_parameter(lines, "Number of on-times")

    if isFisII:
        search_string = "fispact run time"
    else:
        search_string = "CPU Time used for case"
    fo.cpu_time = read_parameter(lines, search_string)

    # find where each time step starts
    time_step_inds = []
    for line in lines:
        if len(line) > 0:
            if line[0:7] == "1 * * *":
                time_step_inds.append(lines.index(line))

    # parse all time steps except setup step
    i = 1
    while i < len(time_step_inds)-1:
        data = lines[time_step_inds[i]:time_step_inds[i+1]]
        fo.timestep_data.append(read_time_step(data, i))
        i = i + 1
    # final timestep
    data = lines[time_step_inds[-1]:]
    fo.timestep_data.append(read_time_step(data, i))

    return fo


def read_time_step(lines, i):
    """ reads a particular time step """
    ts = FispactTimeStep()
    ts.step_num = i + 1

    ts.step_length = float(lines[0][50:60])

    ind = find_ind(lines, "TOTAL NUMBER OF NUCLIDES PRINTED IN INVENTORY")
    ts.num_nuclides = int(lines[ind][50:])

    ind = find_ind(lines, "ALPHA BECQUERELS")
    ts.alpha_act = float(lines[ind][22:34])
    ts.beta_act = float(lines[ind][54:66])
    ts.gamma_act = float(lines[ind][87:99])

    ind = find_ind(lines, "TOTAL ACTIVITY FOR ALL MATERIALS ")
    ts.total_act = float(lines[ind][40:51])

    ind = find_ind(lines, "TOTAL ACTIVITY EXCLUDING TRITIUM ")
    ts.total_act_no_trit = float(lines[ind][40:51])

    ind = find_ind(lines, "TOTAL ALPHA HEAT")
    ts.alpha_heat = float(lines[ind][40:51])
    ts.beta_heat = float(lines[ind + 1][40:51])
    ts.gamma_heat = float(lines[ind + 2][40:51])
    ts.total_heat = float(lines[ind + 2][90:101])
    ts.initial_mass = float(lines[ind + 3][40:51])
    ts.total_heat_no_trit = float(lines[ind + 3][90:101])
    ts.total_mass = float(lines[ind + 4][40:51])
    ts.neutron_flux = float(lines[ind + 5][40:51])

    ts.num_fissions = lines[ind + 6][39:51]
    # added check for E as if <=1E-100 the E is dropped
    if "E" in ts.num_fissions:
        ts.num_fissions = float(ts.num_fissions)

    ts.actinide_burn = float(lines[ind + 6][90:101])

    ind = find_ind(lines, "DENSITY")
    ts.density = float(lines[ind][78:86])

    if ts.total_act > 0.0:
        # added check for E as if <=1E-100 the E is dropped
        ind = find_ind(lines, "APPM OF He  4 ")
        ts.appm_he4 = lines[ind][23:33]
        if "E" in ts.appm_he4:
            ts.appm_he4= float(ts.appm_he4)        
        ts.appm_he3 = lines[ind+1][23:33]
        if "E" in ts.appm_he3:
            ts.appm_he3 = float(ts.appm_he3)
        ts.appm_h3 = lines[ind+2][23:33]
        if "E" in ts.appm_h3:
            ts.appm_h3 = float(ts.appm_h3)
        ts.appm_h2 = lines[ind+3][23:33]
        if "E" in ts.appm_h2:
            ts.appm_h2 = float(ts.appm_h2)
        ts.appm_h1 = lines[ind+4][23:33]
        if "E" in ts.appm_h1:
            ts.appm_h1 = float(ts.appm_h1)
        ind = 1

        ts.dom_data = parse_dominant(lines)
        ts.composition = parse_composition(lines)
        ts.gspec = parse_spectra(lines)
    ts.inventory = parse_inventory(lines)

    return ts


def check_fisp_version(data):
    """ Checks which version of fispact was used to produced data

        requires a list with each element being a line from the fispact output file

        returns a string of the version name
    """

    sub = "FISPACT VERSION 07.0/0"
    data = data[:50]
    if next((s for s in data if sub in s), None):
        v = "FISP07"
    else:
        v = "FISPACT-II"
    return v


def isFisII(data):
    """boolean check if file is fispact-ii output """
    v = check_fisp_version(data)
    if v == "FISPACT-II":
        return True
    else:
        return False


def read_summary_data(data):
    """ Processes the summary block at the end of the file"""

    if isFisII(data):
        cool_str = " -----Irradiation Phase-----"
    else:
        cool_str = "  COOLING STEPS"

    start_ind = data.index(cool_str)
    end_ind = [i for i, line in enumerate(data) if "0 Mass" in line]
    sum_lines = data[start_ind+1:end_ind[0]]
    sum_data = []
    time_yrs = []
    act = []
    act_un = []
    dr = []
    dr_un = []
    heat = []
    heat_un = []
    ing = []
    ing_un = []
    inhal = []
    inhal_un = []
    trit = []
    to = 0

    for l in sum_lines:
        if isFisII(data):
            if l[1] == "-":
                to = time_yrs[-1]
            else:
                time_yrs.append(float(l[24:32]) + to)
                act.append(l[35:43])
                dr.append(l[58:66])
                heat.append(l[81:89])
                ing.append(l[104:112])
                inhal.append(l[127:135])
                trit.append(l[150:158])

        else:
            time_yrs.append(l[20:28])
            act.append(l[31:39])
            dr.append(l[54:62])
            heat.append(l[77:85])
            ing.append(l[100:108])
            inhal.append(l[123:131])
            trit.append(l[146:154])

    sum_data.append(time_yrs)
    sum_data.append(act)
    sum_data.append(dr)
    sum_data.append(heat)
    sum_data.append(ing)
    sum_data.append(inhal)
    sum_data.append(trit)
    sum_data.append(act_un)
    sum_data.append(dr_un)
    sum_data.append(heat_un)
    sum_data.append(ing_un)
    sum_data.append(inhal_un)

    return sum_data


def parse_dominant(data):
    """parse dominant nuclides section and return a list of lists """
    p1_ind = find_ind(data, "DOMINANT NUCLIDES")
    data = data[p1_ind:]
    d1_ind = find_ind(data, "(Bq) ")
    d2_ind = find_ind(data, "GAMMA HEAT")
    topset = data[d1_ind+2:d2_ind-1]
    topset = np.array(topset)
    lowerset = data[d2_ind+3:]

    act_nuc = []
    act = []
    act_percent = []
    heat_nuc = []
    heat = []
    heat_percent = []
    dr_nuc = []
    dr = []
    dr_percent = []
    gheat_nuc = []
    gheat = []
    gheat_percent = []
    bheat_nuc = []
    bheat = []
    bheat_percent = []

    for l in topset:
        act_nuc.append(l[7:13].replace(" ", ""))
        act.append(l[15:25])
        act_percent.append(l[27:36])
        heat_nuc.append(l[38:44].replace(" ", ""))
        heat.append(l[46:56])
        heat_percent.append(l[58:67])
        dr_nuc.append(l[69:75].replace(" ", ""))
        dr.append(l[77:87])
        dr_percent.append(l[89:98])

    for l in lowerset:
        gheat_nuc.append(l[7:13].replace(" ", ""))
        gheat.append(l[15:25])
        gheat_percent.append(l[27:36])
        bheat_nuc.append(l[38:44].replace(" ", ""))
        bheat.append(l[46:56])
        bheat_percent.append(l[58:67])

    dom_data = []
    dom_data.append(act_nuc)
    dom_data.append(act)
    dom_data.append(act_percent)
    dom_data.append(heat_nuc)
    dom_data.append(heat)
    dom_data.append(heat_percent)
    dom_data.append(dr_nuc)
    dom_data.append(dr)
    dom_data.append(dr_percent)
    dom_data.append(gheat_nuc)
    dom_data.append(gheat)
    dom_data.append(gheat_percent)
    dom_data.append(bheat_nuc)
    dom_data.append(bheat)
    dom_data.append(bheat_percent)

    return dom_data


def parse_composition(data):
    """ parse compostions section
        returns a list of 2 lists, one with name of element,
        one with the number of atoms
    """
    p1 = find_ind(data, "COMPOSITION  OF  MATERIAL  BY  ELEMENT")
    p2 = find_ind(data, "GAMMA SPECTRUM AND ENERGIES/SECOND")
    data = data[p1+5:p2-3]
    ele_list = []
    atoms = []

    for l in data:
        ele_list.append(l[12:14])
        atoms.append(float(l[20:30]))

    composition = []
    composition.append(ele_list)
    composition.append(atoms)

    return composition


def parse_spectra(data):
    """ reads gamma spectra data for each timestep
        returns list of length 24 corresponding to 24 gamma energy groups
        data is in gamma/s/cc
    """
    p1 = find_ind(data, "GAMMA SPECTRUM AND ENERGIES/SECOND")
    data = data[p1+7:p1+31]
    spectra = []
    for l in data:
        spectra.append(float(l[130:141]))
    return spectra


def parse_inventory(data):
    """ parse inventory data
        returns a list of lists with all data from the inventory
        section of the times step in order:
        nuclide name,
        # of atoms,
        mass in grams,
        activity in bq,
        beta energy in kw
        alpha energy in kw
        gamma energy in kw
        dose rate in Sv/hr
    """
    inv = []
    p2 = find_ind(data, "0  TOTAL NUMBER OF NUCLIDES PRINTED IN INVENTORY")
    data = data[4:p2]
    for nuc in data:
        nuc_data = [nuc[2:8].replace(" ", ""), float(nuc[14:25]), 
                    float(nuc[28:37]), float(nuc[40:49]), 
                    float(nuc[52:61]), float(nuc[64:72]),
                    float(nuc[75:84]), float(nuc[87:96])]
        inv.append(nuc_data)

    return np.array(inv)


def find_ind(data, sub):
    """ finds index in data whic contains sub string """
    for i, s in enumerate(data):
        if sub in s:
            ind = i
    return ind


def read_parameter(data, sub):
    """ finds and cleans integral values in each timestep"""
    ind = find_ind(data, sub)
    line = data[ind]
    line = line.split("=")
    line = line[1].strip()
    line = line.split(" ")
    param = float(line[0])
    return param

################################
# functions using pypact module
################################
def mesh_to_fispact_fluxin(flux_mesh, flux_tag, fluxin_dir="./", reverse=False,
                           sub_voxel=False, cell_fracs=None, cell_mats=None):
    """This function creates FISPACT-II flux files from fluxes tagged on a PyNE
    Mesh object. Fluxes are printed in the order of the flux_mesh.__iter__().

    Parameters
    ----------
    flux_mesh : PyNE Mesh object
        Contains the mesh with fluxes tagged on each volume element.
    flux_tag : string
        The name of the tag of the flux mesh. Flux values for different energy
        groups are assumed to be represented as vector tags.
    fluxin_dir : string
        The directory of the FISPACT-II flux files.
    reverse : bool
        If true, fluxes will be printed in the reverse order as they appear in
        the flux vector tagged on the mesh.
    sub_voxel: bool, optional
        If true, sub-voxel r2s work flow will be sued. Flux of a voxel will
        be duplicated c times. Where c is the cell numbers of that voxel.
    cell_fracs : structured array, optional
        The output from dagmc.discretize_geom(). A sorted, one dimensional
        array, each entry containing the following fields:

            :idx: int
                The volume element index.
            :cell: int
                The geometry cell number.
            :vol_frac: float
                The volume fraction of the cell withing the mesh ve.
            :rel_error: float
                The relative error associated with the volume fraction.

        The array must be sorted with respect to both idx and cell, with
        cell changing fastest.
    cell_mats : dict, optional
        Maps geometry cell numbers to PyNE Material objects.

        The cell_fracs and cell_mats are used only when sub_voxel=True.
        If sub_voxel=False, neither cell_fracs nor cell_mats will be used.

    """
    tag_flux = flux_mesh.get_tag(flux_tag)

    # find number of e_groups
    e_groups = tag_flux[list(mesh_iterate(flux_mesh.mesh))[0]]
    e_groups = np.atleast_1d(e_groups)
    num_e_groups = len(e_groups)

    # Establish for loop bounds based on if forward or backward printing
    # is requested
    if not reverse:
        start = 0
        stop = num_e_groups
        direction = 1
    else:
        start = num_e_groups - 1
        stop = -1
        direction = -1

    if not sub_voxel:
        for i, mat, ve in flux_mesh:
            output = u""
            # print flux data to file
            output = _output_flux(ve, tag_flux, output, start, stop, direction)
            filename = os.path.join(fluxin_dir, ''.join(["ve", str(i), ".flx"]))
            with open(filename, "w") as f:
                f.write(output)
    else:
        ves = list(flux_mesh.iter_ve())
        for row in cell_fracs:
            if len(cell_mats[row['cell']].comp) != 0:
                output = _output_flux(ves[row['idx']], tag_flux, output, start,
                                      stop, direction)

def _output_flux(ve, tag_flux, output, start, stop, direction):
    """
    This function is used to get neutron flux for fluxin for fispact

    Parameters
    ----------
    ve : entity, a mesh sub-voxel
    tag_flux : array, neutron flux of the sub-voxel
    output : string
    start : int
    stop : int
    direction: int
    """

    count = 0
    flux_data = np.atleast_1d(tag_flux[ve])
    for i in range(start, stop, direction):
        output += u"{:.6E} ".format(flux_data[i])
        # fluxin formatting: create a new line
        # after every 6th entry
        count += 1
        if count % 6 == 0:
            output += u"\n"

    output += u"\n"
    output += u"{:.6E}\n\n".format(1.0)
    return output


def write_fispact_input_single_ve(filename, material):
    """
    """
    id = pp.InputData(name=filename)
    
    # control setup
    id.overwriteExisting()
    id.enableJSON()
    id.approxGammaSpectrum()
    id.readXSData(709)
    id.readDecayData()
    id.enableHalflifeInOutput()
    id.enableHazardsInOutput()
    id.setProjectile(pp.PROJECTILE_NEUTRON)
    id.enableSystemMonitor()
    id.readGammaGroup()
    id.enableInitialInventoryInOutput()
    id.setLogLevel(pp.LOG_SEVERITY_ERROR)
    
    # thresholds
    id.setXSThreshold(1e-12)
    id.setAtomsThreshold(1e5)
    
    ## set target
    #id.setDensity(19.5)
    #id.setMass(1.0)
    #id.addElement('Ti', percentage=80.0)
    #id.addElement('Fe', percentage=14.8)
    #id.addElement('Cr', percentage=5.2)
    id.setDensity(material.density)
    id.setFuel()
    atom_dens = material.to_atom_dens()
    for nuc, comp in material.comp.items():
        # add isotpoes in unit of atoms/kg
        id.addIsotope(name(nuc), atom_dens[nuc]/material.density*1000.0) 
    
    # irradiate and cooling times
    id.addIrradiation(300.0, 1.1e15)
    id.addCooling(10.0)
    id.addCooling(100.0)
    id.addCooling(1000.0)
    id.addCooling(10000.0)
    id.addCooling(100000.0)
    
    # validate data
    id.validate()
    
    #print(pp.to_string(id))
    
    # write to file
    pp.to_file(id, '{}.i'.format(id.name))
