# Generate CellML from BG-RDF models

> [!WARNING]
> These tools are actively being developed and nothing should be considered stable until version `1.0` has been released.

## Installation

### As a package

```sh
uv add https://github.com/CellDL/BG-Tools/releases/download/v0.9.2/bg2cellml-0.9.2-py3-none-any.whl
```

### Development

```sh
$ git clone https://github.com/CellDL/BG-Tools
$ cd BG-Tools
$ git submodule update --init --recursive
$ uv sync --all-groups
```

## Usage

See `cli/bgtool/__main__.py` for details.

### A module for libcuflynx and PhLynx

```sh
$ PYTHONPATH=cli python -m bgtool --module glycogen/bondgraph_v01 --output OUT_DIR glycogen.celldl.svg
```

writes a module version as a module library (such as circulatory-autogen-modules) holds one, in the formats
of libcuflynx's `schemas/module_config.schema.json` (PhLynx's keys) and `schemas/parameters.schema.json`:
`<type>_<subtype>_modules.cellml`, `_units.cellml`, `_modules_config.json` (with `licence`, `creator` and
`required_citations`), `_references.bib`, `_tests.yaml`, `_verification_config.json` and
`parameterisations/default/default_parameters.json` (nested JSON, `schema_version` 1, each value with its units,
data_reference and `sourced: "no"`), every value of the bond graph a parameter there. `--bib FILE` and
`--citation KEY` (repeatable) give the references and required citations (by default a placeholder citing
libcuflynx), `--licence` and `--creator` the config's. `--variable-mapping` also writes
`<type>_<subtype>_variable_mapping.csv`.

Ports: each chemostat species and compartment (all its copies) is a general port,
`chemostat_<species>_<location>`, of its potential (a boundary condition) and the module's uptake of the species. A
signal input port block (`bgf:SignalInputPort`) is an entrance port of its `bgf:portType` (by default
`<name>_port`) whose variable, the block's output, is a boundary condition; a signal output port block
(`bgf:SignalOutputPort`) an exit port (`multi_port`). A flow or potential source with a `bgf:portType` is a port of
its flow and potential, the flow into the module at an entrance and the pressure that of the node the source feeds.
Any other element or block with a `bgf:portType` is a port too: `bgf:portVariables` names the variables it gives
(`"q"` a store's amount, `"v"`, `"u"`, a parameter; by default a store's amount, a block's output, an element's flow),
`bgf:portDirection` makes it an `entrance`, `exit` or `general` port, `bgf:multiPort` sets `multi_port` (`True`,
`sum`), and the elements of one port type and direction with a `bgf:portOrder` make one port of several variables in
that order. A port variable with a value is a boundary condition, a state stays a state. Domain constants no equation
uses (RT in a hydraulic model) are left out. See `bg2cellml/module_export.py`.

`--reference VERSION_DIR` (often `OUTPUT_DIR` itself) regenerates an existing version keeping its form: each
variable's units named as there (where the definitions agree), its parameters' `data_reference`, `sourced` and other
keys (with the bond graph's values; a different value is reported), its config's other keys (e.g.
`required_citations_uncertain`), licence, creators, citations and references, its `tests.yaml` and
`verification_config.json`, and its order of variables, parameters and ports.
`tests/test_module_export.py` builds each example module with libcuflynx and runs it against bg2cellml's own
CellML when `LIBCUFLYNX_SRC` (circulatory_autogen's `src/`) and `LIBCUFLYNX_PYTHON` (a Python with its
dependencies) are set.

`bgtool` checks the CellML by running it with libOpenCOR first; where libOpenCOR's compiler cannot run (it reports
"the given code could not be compiled"), `--save-errors` writes the files anyway.

### Signals

Signal bonds (`bgf:SignalBond`) carry a variable of an element, a junction or a signal block into a signal input of
a template (a signal block's input, a Michaelis-Menten reaction's enzyme amount, a Kedem-Katchalsky flux's osmotic
pressures, ...) or scale one of an element's parameters (`bgf:modulates`, `bgf:modulationExponent`): in its
constitutive relation `p` becomes `p s^n`; a flow or potential source's value `v` becomes `v0 s^n`, `v0` a new
parameter. Signal blocks (BG-RDF's `templates/signal.ttl`) each give one variable, named as the block. Their
variables' units may be their signals' (`bgf:unitsOf`, or none for an input or for a parameter, which then takes its
value's), found once the signals are connected; a value in other compatible units is converted. A signal of the
wrong units into a fixed-unit input, a modulation by a signal with units, and an input with no signal bond are
issues. A signal input port's or a constant's value in other units of the same kind takes the units of the input or
parameter it sets. An element or block may give a domain constant the model's value (`bgf:parameterValue [
bgf:varName "RT" ; ... ]`). A model with `bgf:unitSystem "SI"` has cubic metres for the templates' litres and Pa for
kPa; values in other units of the same kind are converted to their variables'.

Where a bond graph's causality differs from how its equations are written (a resistance between a node a flow source
feeds and a compliance gives the node's pressure), the equations are solved for the variables they must give
(`bg2cellml/causality.py`). Names are written as they are (`alpha`, not a Greek letter). `tests/models/` has an
example of each element: `vessel_inertance.ttl`, `angII_vessel.ttl`,
`baroreceptor.ttl`, `renin_secretion.ttl`, `aldosterone_transport.ttl`, `kedem_katchalsky.ttl`
(`tests/test_signals.py`), `kidney_organ.ttl` (the Na/K-ATPase, arithmetic blocks, ports of several variables) and
`flow_fed_node.ttl` (causality, SI units) (`tests/test_organ_modules.py`).

### A variable mapping for libcuflynx's methods LaTeX

```sh
$ PYTHONPATH=cli python -m bgtool --variable-mapping --output OUT_DIR glycogen.celldl.svg
```

also writes `OUT_DIR/<name>_variable_mapping.csv` beside the CellML, and `--module` writes
`<type>_<subtype>_variable_mapping.csv` with the module. It is libcuflynx's variable mapping
(`cuflynx-methods-latex --mapping`): `variable_name` (`component/variable`), `latex`, `kind`, `component` and
`units`, with one more column, `bondgraph`, each variable's bond graph name (`q_G6P_c`). A variable keeping its
bond graph name is written in bond graph notation (`q^{\mathrm{G6P}}_{c}`); one the model names its own way
(bgf:NamingRule, bgf:VariableName) by libcuflynx's default rule for its new name. See
`bg2cellml/variable_mapping.py`.

### Reaction and species annotations

A model's elements may carry identifiers.org terms (`https://identifiers.org/<prefix>:<id>`), as the CellDL
editor's Annotations panel writes them. A reaction: `bqbiol:isVersionOf` its Rhea reaction (`rhea:17825`) and an
SBO term (`SBO:0000176` biochemical, `SBO:0000655` transport reaction), `bqbiol:hasPart` each Rhea reaction a
lumped reaction combines, and `semsim:hasMediatorParticipant` its enzyme or transporter (`uniprot:P19367`). A
species element (a store, or a chemostat with `bgf:hasSpecies`): `bqbiol:is` its ChEBI term and `bqbiol:isPartOf`
its GO compartment, the same on every element of that species and location. A reaction's participants are the
species elements bonded to it, directly or through 0- and 1-junctions: the bond into the reaction makes them
substrates, out of it products, with the bonds' `bgf:bondCount` as coefficient. `bg2cellml/annotations.py`'s
`model_annotations(model)` reads all this as plain data.

```sh
$ PYTHONPATH=cli python -m bgtool --annotate --output OUT_DIR glycogen.celldl.svg
```

writes `OUT_DIR/<name>.ttl` beside the CellML (and `--omex` puts it in the archive as its OMEX metadata). Each
variable is `bqmodel:isDerivedFrom` its bond graph element. Reactions and species that are annotated, or have
annotated participants, also get [OMEX Metadata 1.2](https://doi.org/10.1515/jib-2021-0020) composite
annotations, with `local:` the annotation file's own namespace and `opb:` `https://identifiers.org/opb:`:
the reaction's flux variable `bqbiol:isVersionOf opb:OPB_00592` (chemical molar flow rate) and
`bqbiol:isPropertyOf local:process_<element>`; the process `bqbiol:isVersionOf` its Rhea and SBO terms (or
`bqbiol:hasPart` each Rhea reaction), with `semsim:hasSourceParticipant local:source_<element>_<k>` (substrates)
and `semsim:hasSinkParticipant local:sink_<element>_<k>` (products), each with `semsim:hasMultiplier` (an
`xsd:double`) and `semsim:hasPhysicalEntityReference local:entity_<species>_<location>`, and a mediator
`local:mediator_<element>` referring to its UniProt term; the entity `bqbiol:is` its ChEBI term and
`bqbiol:isPartOf` its GO term; and a store's amount `bqbiol:isVersionOf opb:OPB_00425` (molar amount of
chemical), `bqbiol:isPropertyOf` its entity. Variables are named by their `cmeta:id`, which naming rules don't
change.

## Example models:

* [An RC electrical circuit](https://github.com/CellDL/BG-RDF/blob/main/examples/example_RC.ttl) 
* [A simple biochemical reaction](https://github.com/CellDL/BG-RDF/blob/main/examples/example_B1.ttl)
* [All examples](https://github.com/CellDL/BG-RDF/blob/main/examples)

### Examples in a Jupyter notebook:

1. Clone the repository and its sub-modules as above.
2. Install with `uv sync --group notebook`
3. Start JupyterLab with `uv run jupyter lab`
4. Open `./notebooks/bg2cellml.ipynb` in Jupyter.
