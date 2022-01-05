#!/usr/bin/env python
#xsimspe - X-ray simulated spectra using XmiMsim

progname = "xsimspe"
version = 0.1

import gi
gi.require_version("XmiMsim", "1.0")
from gi.repository import XmiMsim as xms
import os
import numpy as np
import subprocess
import time

try:
    from mendeleev import element as m_elem
except ImportError:
    print('mendeleev library is required\nInstall it with pip:\n  pip install mendeleev')

if os.name == "posix":
    root_dir = "/home/rosario/progetti/xmimsim"
else:
    root_dir = "C:/Users/utente/"
if not os.path.isdir(root_dir):
    raise ValueError(f'{root_dir} is not a directory or not exist')
    
work_dir = os.path.join(root_dir, "scripts/xsimspe")
if not os.path.exists(work_dir):
    os.makedirs(work_dir, exist_ok = True)

data_dir = os.path.join(root_dir, "outdata")
if not os.path.exists(data_dir):
    os.makedirs(data_dir, exist_ok = True)
    
inputs_dir = os.path.join(data_dir, "input_files")
if not os.path.exists(inputs_dir):
    os.makedirs(inputs_dir, exist_ok = True)
    
NPHOTONSINTERVAL = 5
NPHOTONSLINES = 5
WFRACNUM = 6
MAXNUMCPU = os.cpu_count() - 2

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
    out_fname = f'{progname}-{version}_'
    for items in zip(symbols,np.round(weights, 2)):
        for item in items:
            out_fname += str(item)
        out_fname += "_"
    out_fname += f"{thickness}"
    return out_fname

#_____________________________________________________________________
def gen_input_file(layer_elements, w_fraction, thickness = None, dry_run = False):
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

if __name__ == "__main__":
    if os.name == "posix":
        command = ["/bin/bash", "-c"]
    else:
        command = ["C:/Program Files/XMI-SIM 64-bit/Bin/xmimsim.exe", "--disable_gpu"]
    elements_file = os.path.join(work_dir, "elements.txt")
    input_elements = read_elements(elements_file)
    # to do - loop on input_elements
    layer_elements = input_elements[0]
    num_elements = len(layer_elements.items())
    # w_fraction loop
    wfrac = [np.linspace(0,1,WFRACNUM) for _ in range(num_elements)]
    w_fraction = np.array(np.meshgrid(*wfrac))
    w_fraction = w_fraction.T.reshape(-1,num_elements)[1:]
    processes = set()
    print(f"using {MAXNUMCPU} cores")
    for proc_num, weights in enumerate(w_fraction):
        _file = gen_input_file(layer_elements, weights, dry_run = True)
        #command_string = f"echo proc: {proc_num:3} {_file}"
        command_string = _file
        processes.add(subprocess.Popen(command + [command_string]))
        while len(processes) >= MAXNUMCPU:
            #time.sleep(0.5)
            processes.difference_update([
                p for p in processes if p.poll() is not None])
        
        #print(gen_input_file(layer_elements, weights, dry_run = True))

