from __future__ import annotations

from rdkit import Chem

from chem_predict.chemistry import canonicalize_smiles, mol_from_smiles
from chem_predict.medchem.models import (
    MedchemAssessment,
    MedchemHit,
    MedchemRule,
    MedchemSeverity,
)


LILLY_RULESET_VERSION = "2.1.0"
LILLY_SOURCE_URL = "https://github.com/IanAWatson/Lilly-Medchem-Rules"

# Initial RDKit-native port. This is deliberately a small, auditable subset.
# Each entry points back to the upstream query name. Additional rules should be
# added only with a regression fixture against the official implementation.
LILLY_NATIVE_RULES: tuple[MedchemRule, ...] = (
    MedchemRule(
        id="lilly_nitrosamine",
        name="Nitrosamine",
        smarts="[#7]-[N]=O",
        severity=MedchemSeverity.REJECT,
        source_rule="queries/nitrosamine.qry",
        source_version=LILLY_RULESET_VERSION,
        port_status="direct_smarts_port",
    ),
    MedchemRule(
        id="lilly_acid_halide",
        name="Acid halide",
        smarts="[#6X3](=[#8])-[#9,#17,#35,#53]",
        severity=MedchemSeverity.REJECT,
        source_rule="queries/acid_halide.qry",
        source_version=LILLY_RULESET_VERSION,
        port_status="rdkit_translation",
    ),
    MedchemRule(
        id="lilly_aldehyde",
        name="Aldehyde / thioaldehyde",
        smarts="[#6H1](=[#8,#16])",
        severity=MedchemSeverity.REJECT,
        source_rule="queries/aldehyde.qry",
        source_version=LILLY_RULESET_VERSION,
        port_status="rdkit_translation",
    ),
    MedchemRule(
        id="lilly_ester",
        name="Ester",
        smarts="[O;D1]=[C;D3;R0]-[O;D2;!$([O]-c)]",
        severity=MedchemSeverity.DEMERIT,
        demerit=35,
        source_rule="queries/ester.qry",
        source_version=LILLY_RULESET_VERSION,
        port_status="rdkit_translation",
        notes="Translated from LillyMol query syntax; excludes aryl esters via the upstream Environment_no_Match condition.",
    ),
    MedchemRule(
        id="lilly_nitro",
        name="Nitro group",
        smarts="[$([N+](=O)[O-]),$([N](=O)=O)]",
        severity=MedchemSeverity.DEMERIT,
        demerit=60,
        source_rule="queries/nitro.qry",
        source_version=LILLY_RULESET_VERSION,
        port_status="rdkit_approximation",
        notes="Upstream query contains extra topology/environment constraints; validate against official LillyMol output.",
    ),
)


class LillyNativeRules:
    """RDKit-native, incrementally ported Lilly Medchem Rules.

    Not yet equivalent to the complete upstream ~300-rule implementation.
    """

    def __init__(
        self,
        *,
        demerit_cutoff: int = 100,
        min_heavy_atoms: int = 7,
        hard_max_heavy_atoms: int = 40,
        rules: tuple[MedchemRule, ...] = LILLY_NATIVE_RULES,
    ) -> None:
        self.demerit_cutoff = demerit_cutoff
        self.min_heavy_atoms = min_heavy_atoms
        self.hard_max_heavy_atoms = hard_max_heavy_atoms
        self.rules = rules

    def assess(self, smiles: str) -> MedchemAssessment:
        mol = mol_from_smiles(smiles)
        canonical = canonicalize_smiles(smiles)
        warnings: list[str] = ["partial_lilly_rdkit_port_not_full_ruleset"]
        hits: list[MedchemHit] = []
        rejected = False

        heavy_atoms = mol.GetNumHeavyAtoms()
        if heavy_atoms < self.min_heavy_atoms:
            rejected = True
            hits.append(
                MedchemHit(
                    rule_id="lilly_too_few_atoms",
                    rule_name=f"Fewer than {self.min_heavy_atoms} heavy atoms",
                    severity=MedchemSeverity.REJECT,
                    count=1,
                    demerit_each=0,
                    demerit_total=0,
                )
            )
        if heavy_atoms > self.hard_max_heavy_atoms:
            rejected = True
            hits.append(
                MedchemHit(
                    rule_id="lilly_too_many_atoms_hard",
                    rule_name=f"More than {self.hard_max_heavy_atoms} heavy atoms",
                    severity=MedchemSeverity.REJECT,
                    count=1,
                    demerit_each=0,
                    demerit_total=0,
                )
            )

        total_demerits = 0
        for rule in self.rules:
            query = Chem.MolFromSmarts(rule.smarts)
            if query is None:
                raise ValueError(f"Invalid internal SMARTS for {rule.id}: {rule.smarts}")
            matches = mol.GetSubstructMatches(query, uniquify=True)
            if not matches:
                continue
            count = len(matches)
            subtotal = rule.demerit * count
            total_demerits += subtotal
            if rule.severity == MedchemSeverity.REJECT:
                rejected = True
            hits.append(
                MedchemHit(
                    rule_id=rule.id,
                    rule_name=rule.name,
                    severity=rule.severity,
                    count=count,
                    demerit_each=rule.demerit,
                    demerit_total=subtotal,
                    atom_matches=tuple(tuple(int(i) for i in m) for m in matches),
                )
            )
            if rule.port_status == "rdkit_approximation":
                warnings.append(f"{rule.id}_is_approximate")

        if total_demerits >= self.demerit_cutoff:
            rejected = True

        return MedchemAssessment(
            smiles=canonical,
            passed=not rejected,
            rejected=rejected,
            total_demerits=total_demerits,
            demerit_cutoff=self.demerit_cutoff,
            hits=tuple(hits),
            warnings=tuple(dict.fromkeys(warnings)),
            metadata={
                "ruleset": "Lilly Medchem Rules",
                "ruleset_version": LILLY_RULESET_VERSION,
                "source": LILLY_SOURCE_URL,
                "heavy_atoms": heavy_atoms,
            },
        )
