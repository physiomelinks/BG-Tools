"""
libcuflynx_run.py MODULE_DIR TYPE SUBTYPE ORIGINAL.cellml WORK_DIR

Run by tests/test_module_export.py in a Python with libcuflynx (circulatory_autogen's src/ on PYTHONPATH) and
Myokit: builds a one-instance model ("bg") of the exported module version from its default parameterisation,
with an empty host parameters file, runs it and bg2cellml's own CellML from 1.5 times their initial states,
and prints {"states", "difference", "moved"} as JSON: the largest difference of a state relative to its
largest value, and the largest relative change of a state over the run.
"""
import contextlib
import io
import json
import os
import sys

import numpy as np

module_dir, module_type, subtype, original, work = sys.argv[1:6]

from libcuflynx.scripts.script_generate_with_new_architecture import generate_with_new_architecture
import myokit
import myokit.formats.cellml as mc

resources = os.path.join(work, 'resources')
os.makedirs(resources, exist_ok=True)
prefix = 'bgcheck'
with open(os.path.join(resources, f'{prefix}_module_array.json'), 'w') as f:
    json.dump([{'name': 'bg', 'module_type': module_type, 'module_subtype': subtype}], f)
with open(os.path.join(resources, f'{prefix}_parameters.json'), 'w') as f:
    json.dump({'schema_version': 1}, f)
config = {
    'file_prefix': prefix, 'input_param_file': f'{prefix}_parameters.json', 'model_type': 'cellml',
    'solver': 'CVODE_myokit', 'resources_dir': resources,
    'generated_models_dir': os.path.join(work, 'generated_models'),
    'module_library_dirs': [module_dir], 'use_builtin_modules': False, 'DEBUG': False,
}
log = io.StringIO()
with contextlib.redirect_stdout(log):
    ok = generate_with_new_architecture(False, config)
if not ok:
    sys.exit(f'generation failed:\n{log.getvalue()[-3000:]}')
generated = os.path.join(work, 'generated_models', prefix, f'{prefix}_flat.cellml')

def run(path, rename=lambda name: name):
    model = mc.CellMLImporter().model(path)
    states = {rename(v.name()): v.qname() for v in model.states()}
    simulation = myokit.Simulation(model)
    simulation.set_state([1.5 * x for x in simulation.state()])
    scale = max(abs(x) for x in simulation.state()) or 1
    simulation.set_tolerance(1e-12 * scale, 1e-10)
    result = simulation.run(10.0, log=list(states.values()), log_interval=0.1)
    return {name: np.asarray(result[q]) for name, q in states.items()}

a = run(original)
b = run(generated, rename=lambda name: name[:-len('_bg')] if name.endswith('_bg') else name)
if set(a) - set(b):
    sys.exit(f'states missing from the libcuflynx model: {sorted(set(a) - set(b))} (it has {sorted(b)})')
difference = max(float(np.max(np.abs(b[s] - a[s])) / (np.max(np.abs(a[s])) or 1)) for s in a)
moved = max(float(np.max(np.abs(a[s] - a[s][0])) / (np.max(np.abs(a[s])) or 1)) for s in a)
print(json.dumps({'states': len(a), 'difference': difference, 'moved': moved}))
