# Medchem and reaction-informatics integrations

## Lilly Medchem Rules

Upstream: https://github.com/IanAWatson/Lilly-Medchem-Rules

Version checked for this integration: **2.1.0** (2026-06-24).

The upstream implementation is Apache-2.0 licensed and contains roughly 300
rejection/demerit rules. Its current Python driver orchestrates LillyMol
executables and query files; it is not an RDKit package.

Chem-Predict therefore exposes two backends:

### RDKit-native incremental port

```python
from chem_predict.medchem import LillyNativeRules

result = LillyNativeRules().assess("CCN(N=O)CC")
print(result.passed)
print(result.total_demerits)

for hit in result.hits:
    print(hit.rule_id, hit.severity, hit.demerit_total)
```

The native port currently includes:

- default lower and hard upper heavy-atom rejection thresholds;
- nitrosamine hard rejection;
- acid-halide hard rejection;
- aldehyde/thioaldehyde hard rejection;
- ester demerit (35 per match);
- nitro demerit (60 per match).

This is **not yet the complete Lilly ruleset**. Every ported rule records its
upstream query name, upstream version and port status. Rules translated from
non-SMARTS Lilly query syntax are explicitly distinguishable from direct
SMARTS ports.

The design goal is to expand this table incrementally and validate each port
against the upstream test molecules and/or official LillyMol output.

### Official backend adapter

If an upstream Lilly-Medchem-Rules checkout and its required LillyMol
executables are available locally:

```python
from chem_predict.medchem import LillyOfficialBackend

backend = LillyOfficialBackend("/path/to/Lilly-Medchem-Rules")
result = backend.run_smiles("CCN(N=O)CC")
print(result.returncode)
print(result.stdout)
```

This keeps exact-upstream reproducibility available while the RDKit-native
rules are ported.

## SynKit

Upstream: https://github.com/TieuLongPhan/SynKit

Release/API checked: **SynKit 1.6.2**, released 2026-08-18.

SynKit is kept optional because its full dependency set is substantially
larger than the Chem-Predict core.

Install:

```bash
python -m pip install -e ".[synkit]"
```

The first adapter uses the documented/current tuple ITS interface:

```python
from chem_predict.integrations import reaction_to_its

reaction = "[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]"
result = reaction_to_its(reaction)

for change in result.changed_bonds:
    print(change)
```

Internally the adapter calls:

```python
synkit.IO.rsmi_to_its(
    reaction_smiles,
    format="tuple",
    ...
)
```

and extracts changed bonds whose paired `order=(reactant, product)` differs.

This boundary is intentionally narrow. Future uses can add adapters for:

- reaction-centre extraction;
- graph-transformation rule extraction/application;
- reaction-network analysis;
- supplied mechanism verification;
- Lewis-labelled/electron-flow graphs.

Chem-Predict should not expose SynKit internals throughout the rest of the
package. If SynKit changes API, only this adapter should require modification.
