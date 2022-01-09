#!/usr/bin/env python
#xsimspe - X-ray simulated spectra using XmiMsim

progname = "xsimspe"
version = 0.1

try:
    import gi
except ImportError:
    USEAPI = False
else:
    USEAPI = True
    gi.require_version("XmiMsim", "1.0")
    from gi.repository import XmiMsim as xms
import os
import numpy as np
import subprocess
import time
import xml.etree.ElementTree as et
from xml.dom import minidom

try:
    from mendeleev import element as m_elem
except ImportError:
    print('mendeleev library is required\nInstall it with pip:\n  pip install mendeleev')

if os.name == "posix":
    dtd_file = "/usr/share/xmimsim/xmimsim-1.0.dtd" if USEAPI else "http://www.xmi.UGent.be/xml/xmimsim-1.0.dtd"
    user = os.getenv('USER')
    if user == 'rosario':
        root_dir = "/home/rosario/progetti/xmimsim"
        work_dir = os.path.join(root_dir, "scripts/xsimspe")
    else:
        root_dir = f"/home/{user}/rosario_sim"
        work_dir = os.path.join(root_dir, "xsimspe")
else:
    dtd_file = "http://www.xmi.UGent.be/xml/xmimsim-1.0.dtd"
    root_dir = "C:\\Users\\XRAYLab\\rosario_sim"
    work_dir = os.path.join(root_dir, "xsimspe")

if not os.path.exists(root_dir):
    os.makedirs(root_dir, exist_ok = True)

if not os.path.exists(work_dir):
    os.makedirs(work_dir, exist_ok = True)

data_dir = os.path.join(root_dir, "outdata")
if not os.path.exists(data_dir):
    os.makedirs(data_dir, exist_ok = True)
    
inputs_dir = os.path.join(data_dir, "input_files")
if not os.path.exists(inputs_dir):
    os.makedirs(inputs_dir, exist_ok = True)
    
NPHOTONSINTERVAL = 100
NPHOTONSLINES = 1000
WFRACNUM = 3
MAXNUMCPU = os.cpu_count() if os.cpu_count() <= 4 else os.cpu_count() -2
SINGLECORE = False
DRYRUN = False

#_____________________________________________________________________
def read_elements(file_to_read : str) -> list :
    if not isinstance(file_to_read, str):
        raise TypeError(f'String expected: {file_to_read}')
    elements = list()
    with open(file_to_read) as f2r:
        for _ in f2r:
            _ = _.strip()
            print(f"loading  {_:12} data")
            layer_elements = {}
            for _elem in _.split(sep = ','):
                _elem = m_elem(_elem.strip())
                layer_elements[_elem.symbol] = _elem
            elements.append(layer_elements)
    return elements

#_____________________________________________________________________
def layer_density(elements_densities, weight_fractions):
    return np.average(elements_densities, weights = weight_fractions)
    
#_____________________________________________________________________
def gen_out_file_name(symbols, weights, thickness):
    out_fname = f'{progname}_'
    for items in zip(symbols,np.round(weights, 2)):
        for item in items:
            out_fname += str(item)
        out_fname += "_"
    out_fname += f"{thickness}_hydrocerussite"
    return out_fname

#_____________________________________________________________________
def gen_input_file_from_API(layer_elements, w_fraction, thickness = None, dry_run = False):
    input_template_fname = os.path.join(work_dir, "input_template.xmsi")
    input_template = xms.Input.read_from_xml_file(input_template_fname)
    new_input = xms.Input.init_empty()

    # COMPOSITION
    e_symbol = []
    e_density = []
    atomic_num = []
    for k, v in layer_elements.items():
        e_symbol.append(k)
        e_density.append(v.density)
        atomic_num.append(v.atomic_number)
        
    n_layers = input_template.composition.n_layers
    ref_layer_index = input_template.composition.reference_layer -1
    layers = [input_template.composition.get_layer(index) for index in range(n_layers)]
    if not thickness:
        thickness = layers[ref_layer_index].thickness
    new_ref_layer = xms.Layer.new(atomic_num,
                                  w_fraction,
                                  layer_density(e_density,w_fraction),
                                  thickness)
    layers[ref_layer_index] = new_ref_layer
    composition = xms.Composition.new(layers, reference_layer = ref_layer_index + 1)
    new_input.set_composition(composition)

    # GENERAL
    out_file_name = gen_out_file_name(e_symbol,w_fraction,layers[ref_layer_index].thickness)
    output_file = os.path.join(data_dir, out_file_name + ".xmso")
    general = input_template.general.copy()
    general.n_photons_interval = NPHOTONSINTERVAL
    general.n_photons_lines = NPHOTONSLINES
    general.outputfile = output_file
    new_input.set_general(general)

    # GEOMETRY
    new_input.set_geometry(input_template.geometry.copy())

    # EXCITATION
    new_input.set_excitation(input_template.excitation.copy())

    # ABSORBERS
    new_input.set_absorbers(input_template.absorbers.copy())

    # DETECTOR
    new_input.set_detector(input_template.detector.copy())
    if not dry_run:
        new_input.write_to_xml_file(os.path.join(inputs_dir, out_file_name + ".xmsi"))
    return os.path.join(inputs_dir, out_file_name + ".xmsi")

#______________________________________________________________________
def gen_input_file(layer_elements, w_fraction, thickness = None, dry_run = False):
    lt = time.localtime()
    timestamp = f'{lt.tm_year}-{lt.tm_mon}-{lt.tm_mday} {lt.tm_hour}:{lt.tm_min}:{lt.tm_sec} (CET)'
    if not thickness:
        thickness = 0.002
    input_template_fname = os.path.join(work_dir, "input_template.xmsi.in")
    e_symbol = []
    e_density = []
    atomic_num = []
    for k, v in layer_elements.items():
        e_symbol.append(k)
        e_density.append(v.density)
        atomic_num.append(v.atomic_number)
    reference_layer = et.Element("layer")
    for a, w in zip(atomic_num, w_fraction):
        element = et.Element("element")
        reference_layer.append(element)
        atnum = et.SubElement(element, "atomic_number")
        atnum.text = str(a)
        wfrac = et.SubElement(element, "weight_fraction")
        wfrac.text = str(w)
    density = et.SubElement(reference_layer, "density")
    density.text = str(layer_density(e_density, w_fraction))
    thness = et.SubElement(reference_layer, "thickness")
    thness.text = str(thickness)
    
    #xmlstr = minidom.parseString(et.tostring(reference_layer)).toprettyxml(indent = "  ")
    #xmlstr = xmltxt[xmltxt.find('\n') + 1 : -1]
    xmlstr = et.tostring(reference_layer).decode('utf-8')
    out_file_name = gen_out_file_name(e_symbol,w_fraction,thickness)
    output_file = os.path.join(data_dir, out_file_name + ".xmso")
    input_template_str = ""
    with open(input_template_fname) as it:
        for line in it:
            line = line.replace('\n', '')
            if "@outputfile@" in line:
                line = line.replace("@outputfile@", output_file)
            elif "@n_photons_interval@" in line:
                line = line.replace("@n_photons_interval@", str(NPHOTONSINTERVAL))
            elif "@n_photons_line@" in line:
                line = line.replace("@n_photons_line@", str(NPHOTONSLINES))
            elif "@reference_layer@" in line:
                line = line.replace("@reference_layer@", xmlstr)
            elif "@dtd_file@" in line:
                line = line.replace("@dtd_file@", dtd_file)
            elif "@timestamp@" in line:
                line = line.replace("@timestamp@", timestamp)
            input_template_str += line.strip()
    input_template_str = minidom.parseString(input_template_str).toprettyxml(indent = "  ")        
    if not dry_run:
        with open(os.path.join(inputs_dir, out_file_name + ".xmsi"), "w") as inout:
            inout.write(input_template_str)
    return os.path.join(inputs_dir, out_file_name + ".xmsi")

#_______________________________________________________________________
def singlecore_processing():
    """One input file for core"""
    if os.name == "posix":
        command = ["xmimsim", "--disable-gpu", "--set-threads", "1"]
    else:
        command = ["C:\\Program Files\\XMI-MSIM 64-bit\\Bin\\xmimsim-cli.exe",
                   "--disable-gpu",
                   "--set-threads", "1"]
    processes = set()
    # w_fraction loop
    print(f"using {MAXNUMCPU} cores")
    for proc_num, weights in enumerate(w_fraction):
        #print(gen_input_file(layer_elements, weights, dry_run = True))
        if USEAPI:
            Ifile = gen_input_file_from_API(layer_elements, weights)
        else:
            Ifile = gen_input_file(layer_elements, weights)
        #command_string = f"echo proc: {proc_num:3} {Ifile}"
        if os.name == 'posix':
            processes.add(subprocess.Popen(command + [Ifile]))
        else:
            processes.add(subprocess.Popen(command + [Ifile], shell = True))
        print(f"processing {Ifile}")
        while len(processes) >= MAXNUMCPU:
            time.sleep(0.5)
            processes.difference_update([
                p for p in processes if p.poll() is not None])
    for p in processes:
        p.wait()

#_____________________________________________________________________
def multicore_processing(disable_gpu = True):
    """Multicore for each input file at a time"""
    if os.name == 'posix':
        command = ["xmimsim"]
    else:
        command = ["C:\\Program Files\\XMI-MSIM 64-bit\\Bin\\xmimsim-cli.exe"]
    if disable_gpu:
        command  += ["--set-threads", str(MAXNUMCPU), "--disable-gpu"]
    else:
        command += ["--enable-gpu"]
    #weight fractions loop
    for weights in w_fraction:
        if USEAPI:
            Ifile = gen_input_file_from_API(layer_elements, weights)
        else:
            Ifile = gen_input_file(layer_elements, weights, dry_run = True)
        print(f"processing {Ifile} with {MAXNUMCPU} core")
        p = subprocess.Popen(command + [Ifile])
        p.wait()
    
if __name__ == "__main__":
    elements_file = os.path.join(work_dir, "elements.txt")
    input_elements = read_elements(elements_file)
    # to do - loop on input_elements
    layer_elements = input_elements[1]
    num_elements = len(layer_elements.items())
    wfrac = [np.linspace(0,1,WFRACNUM) for _ in range(num_elements)]
    w_fraction = np.array(np.meshgrid(*wfrac))
    w_fraction = w_fraction.T.reshape(-1,num_elements)[1:]
    #normalize w_fraction
    w_fraction_sum = w_fraction.sum(axis = 1)
    for c in range(w_fraction.shape[1]):
        w_fraction[:,c] = w_fraction[:,c]/w_fraction_sum
    w_fraction = np.unique(w_fraction, axis = 0)
    USEAPI = False
    start = time.time()
    if SINGLECORE:
        singlecore_processing()
    else:
        multicore_processing()
    stop = time.time()
    print(f'time enlapsed: {stop - start}')

