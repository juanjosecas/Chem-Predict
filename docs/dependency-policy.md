# Dependency verification policy

Chem-Predict must not rely on remembered third-party APIs when adding or
modifying integration code.

For each external library used by a new feature:

1. check the current official documentation or maintained upstream API
   reference;
2. record the checked version/date and the specific APIs used;
3. add executable tests around the integration boundary;
4. avoid deprecated aliases when the current supported API is available;
5. prefer a narrow adapter boundary so future library changes do not require
   refactoring domain logic.

## Current verification log

### RDKit

Checked: 2026-09-23  
Documentation release checked: RDKit 2026.03.6  
Package minimum currently supported by Chem-Predict: RDKit >= 2025.09.4

APIs used by the initial chemistry and nitrosamine modules include:

- \`Chem.MolFromSmiles\`
- \`Chem.MolToSmiles\`
- \`Chem.MolFromSmarts\`
- \`Chem.SanitizeMol\`
- \`Atom.GetTotalNumHs\`
- \`Mol.GetSubstructMatches\`
- \`Mol.GetRingInfo\` / \`RingInfo.AtomRings\`
- \`rdChemReactions.ReactionFromSmarts\`
- \`ChemicalReaction.RunReactants\`

RDKit reaction products are explicitly sanitized before Chem-Predict accepts
them as generated products.

Whenever an integration is extended, re-check the upstream documentation even
if the same library already appears in this file.


### RDKit visualization

Checked: 2026-09-23
Documentation release checked: RDKit 2026.03.6

Visualization APIs used:

- `Draw.MolToSVG`
- `rdMolDraw2D.MolDraw2DSVG`
- `rdMolDraw2D.MolDraw2DCairo`
- `MolDraw2D.DrawReaction`
- `rdChemReactions.ReactionFromSmarts(..., useSmiles=...)`

Reaction-center detection uses stable RDKit atom-map and bond APIs rather than
depending on undocumented drawing internals.
