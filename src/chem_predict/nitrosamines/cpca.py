from __future__ import annotations

from rdkit import Chem

from chem_predict.chemistry import mol_from_smiles
from chem_predict.nitrosamines.models import (
    CPCAAssessment,
    CPCAFeature,
    CPCASiteResult,
    NitrosamineSite,
)
from chem_predict.nitrosamines.sites import find_nitrosamine_sites


FDA_CPCA_SOURCE = "FDA Recommended Acceptable Intake Limits for NDSRIs (August 2023)"
FDA_CPCA_SOURCE_VERSION = "FDA web cross-check 2026-08-06"
FDA_CPCA_URL = "https://www.fda.gov/media/170794/download"

CATEGORY_AI_NG_PER_DAY: dict[int, float] = {
    1: 26.5,
    2: 100.0,
    3: 400.0,
    4: 1500.0,
    5: 1500.0,
}

_CARBOXYLIC_ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
if _CARBOXYLIC_ACID is None:  # pragma: no cover
    raise RuntimeError("Could not compile internal carboxylic-acid SMARTS")


def _feature(feature_id: str, label: str, score: int) -> CPCAFeature:
    return CPCAFeature(id=feature_id, label=label, score=score)


def _has_tertiary_alpha_carbon(mol: Chem.Mol, site: NitrosamineSite) -> bool:
    for idx in site.alpha_carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        carbon_neighbors = sum(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())
        if carbon_neighbors == 3:
            return True
    return False


def _methylene_alpha_is_ethyl(mol: Chem.Mol, site: NitrosamineSite) -> bool:
    """FDA Table A footnote for the 0,2 alpha-H pattern."""

    central_n = site.nitrosamine_n_index
    for idx in site.alpha_carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetTotalNumHs() != 2:
            continue
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() == central_n or neighbor.GetAtomicNum() != 6:
                continue
            if neighbor.GetHybridization() == Chem.HybridizationType.SP3 and neighbor.GetTotalNumHs() == 3:
                return True
    return False


def _alpha_hydrogen_score(mol: Chem.Mol, site: NitrosamineSite) -> int | None:
    pair = tuple(sorted(site.alpha_hydrogens))
    mapping = {
        (0, 2): 3,
        (0, 3): 2,
        (1, 2): 3,
        (1, 3): 3,
        (2, 2): 1,
        (2, 3): 1,
    }
    score = mapping.get(pair)
    if pair == (0, 2) and score is not None and _methylene_alpha_is_ethyl(mol, site):
        return 2
    return score


def _ring_feature(mol: Chem.Mol, site: NitrosamineSite) -> CPCAFeature | None:
    n_idx = site.nitrosamine_n_index
    rings = sorted((ring for ring in mol.GetRingInfo().AtomRings() if n_idx in ring), key=len)
    if not rings:
        return None

    ring = rings[0]
    size = len(ring)
    atomic_numbers = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in ring]

    if size == 5 and sorted(atomic_numbers) == [6, 6, 6, 6, 7]:
        return _feature("pyrrolidine_ring", "N-nitroso group in a pyrrolidine ring", 3)
    if size == 6 and 16 in atomic_numbers:
        return _feature(
            "six_membered_sulfur_ring",
            "N-nitroso group in a 6-membered ring containing sulfur",
            3,
        )
    if size == 6 and sorted(atomic_numbers) == [6, 6, 6, 6, 7, 8]:
        return _feature("morpholine_ring", "N-nitroso group in a morpholine ring", 1)
    if size in {5, 6}:
        return _feature("five_or_six_membered_ring", "N-nitroso group in a 5- or 6-membered ring", 2)
    if size == 7:
        return _feature("seven_membered_ring", "N-nitroso group in a 7-membered ring", 1)
    return None


def _path_qualifies_for_long_chain(mol: Chem.Mol, path: tuple[int, ...]) -> bool:
    if len(path) < 5:
        return False
    first_five = set(path[:5])
    for ring in mol.GetRingInfo().AtomRings():
        if len(first_five.intersection(ring)) > 4:
            return False
    return True


def _has_five_heavy_atom_path(mol: Chem.Mol, start_idx: int, blocked: set[int]) -> bool:
    stack: list[tuple[int, tuple[int, ...]]] = [(start_idx, (start_idx,))]
    while stack:
        current, path = stack.pop()
        if len(path) >= 5 and _path_qualifies_for_long_chain(mol, path):
            return True
        if len(path) >= 5:
            continue
        atom = mol.GetAtomWithIdx(current)
        for neighbor in atom.GetNeighbors():
            idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 1 or idx in blocked or idx in path:
                continue
            stack.append((idx, path + (idx,)))
    return False


def _long_chain_feature(mol: Chem.Mol, site: NitrosamineSite) -> CPCAFeature | None:
    central_n = mol.GetAtomWithIdx(site.nitrosamine_n_index)
    if central_n.IsInRing():
        return None

    blocked = {site.nitrosamine_n_index, site.nitroso_n_index, site.oxygen_index}
    if all(_has_five_heavy_atom_path(mol, alpha, blocked) for alpha in site.alpha_carbon_indices):
        return _feature(
            "long_chain_both_sides",
            "At least 5 consecutive non-hydrogen atoms on both sides of an acyclic N-nitroso group",
            1,
        )
    return None


def _is_carboxylic_acid_carbon(mol: Chem.Mol, carbon_idx: int) -> bool:
    atom = mol.GetAtomWithIdx(carbon_idx)
    if atom.GetAtomicNum() != 6:
        return False
    oxygen_double = False
    hydroxyl_oxygen = False
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() != 8:
            continue
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            oxygen_double = True
        elif bond.GetBondType() == Chem.BondType.SINGLE and other.GetTotalNumHs() >= 1:
            hydroxyl_oxygen = True
    return oxygen_double and hydroxyl_oxygen


def _carbonyl_like_ewg(mol: Chem.Mol, neighbor_idx: int, alpha_idx: int) -> bool:
    atom = mol.GetAtomWithIdx(neighbor_idx)
    if atom.GetAtomicNum() != 6 or _is_carboxylic_acid_carbon(mol, neighbor_idx):
        return False

    has_double_hetero = False
    non_alpha_heavy: list[Chem.Atom] = []
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetIdx() == alpha_idx:
            continue
        if other.GetAtomicNum() > 1:
            non_alpha_heavy.append(other)
        if bond.GetBondType() == Chem.BondType.DOUBLE and other.GetAtomicNum() in {7, 8, 16}:
            has_double_hetero = True
    if not has_double_hetero:
        return False

    # FDA notes conflicting data for ketones; only retain carbonyl-like groups
    # carrying an additional heteroatom substituent.
    return any(a.GetAtomicNum() not in {1, 6} for a in non_alpha_heavy)


def _alpha_has_ewg(mol: Chem.Mol, site: NitrosamineSite, alpha_idx: int) -> bool:
    alpha = mol.GetAtomWithIdx(alpha_idx)
    central_n = site.nitrosamine_n_index

    for neighbor in alpha.GetNeighbors():
        idx = neighbor.GetIdx()
        if idx == central_n:
            continue

        z = neighbor.GetAtomicNum()
        if z in {9, 17, 35, 53}:  # halogens
            return True
        if z in {7, 8, 16}:  # directly attached N/O/S substituents
            return True
        if neighbor.GetIsAromatic():
            continue  # aryl is an activating feature counted separately
        if z != 6:
            continue

        if sum(n.GetAtomicNum() == 9 for n in neighbor.GetNeighbors()) >= 3:
            return True

        for bond in neighbor.GetBonds():
            other = bond.GetOtherAtom(neighbor)
            if bond.GetBondType() == Chem.BondType.TRIPLE and other.GetAtomicNum() == 7:
                return True

        if _carbonyl_like_ewg(mol, idx, alpha_idx):
            return True

        # Kept finite and explicit. This catalogue must be validated/expanded
        # independently rather than silently treating every unsaturation as EWG.
        for bond in neighbor.GetBonds():
            other = bond.GetOtherAtom(neighbor)
            if other.GetIdx() == alpha_idx or other.GetAtomicNum() != 6:
                continue
            if bond.GetBondType() in {Chem.BondType.DOUBLE, Chem.BondType.TRIPLE}:
                return True

    return False


def _ewg_feature(mol: Chem.Mol, site: NitrosamineSite) -> CPCAFeature | None:
    sides = sum(_alpha_has_ewg(mol, site, idx) for idx in site.alpha_carbon_indices)
    if sides == 1:
        return _feature(
            "alpha_ewg_one_side",
            "Electron-withdrawing group bonded to alpha-carbon on one side",
            1,
        )
    if sides >= 2:
        return _feature(
            "alpha_ewg_both_sides",
            "Electron-withdrawing groups bonded to alpha-carbons on both sides",
            2,
        )
    return None


def _alpha_has_beta_hydroxyl(mol: Chem.Mol, site: NitrosamineSite, alpha_idx: int) -> bool:
    alpha = mol.GetAtomWithIdx(alpha_idx)
    for beta in alpha.GetNeighbors():
        if beta.GetIdx() == site.nitrosamine_n_index or beta.GetAtomicNum() != 6:
            continue
        if beta.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        for neighbor in beta.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() >= 1:
                return True
    return False


def _beta_hydroxyl_feature(mol: Chem.Mol, site: NitrosamineSite) -> CPCAFeature | None:
    sides = sum(_alpha_has_beta_hydroxyl(mol, site, idx) for idx in site.alpha_carbon_indices)
    if sides == 1:
        return _feature(
            "beta_hydroxyl_one_side",
            "Hydroxyl group bonded to beta-carbon on one side",
            1,
        )
    if sides >= 2:
        return _feature(
            "beta_hydroxyl_both_sides",
            "Hydroxyl groups bonded to beta-carbons on both sides",
            2,
        )
    return None


def _aryl_alpha_feature(mol: Chem.Mol, site: NitrosamineSite) -> CPCAFeature | None:
    for alpha_idx in site.alpha_carbon_indices:
        alpha = mol.GetAtomWithIdx(alpha_idx)
        for neighbor in alpha.GetNeighbors():
            if neighbor.GetIdx() == site.nitrosamine_n_index:
                continue
            if neighbor.GetIsAromatic():
                return _feature("aryl_on_alpha", "Aryl group bonded to alpha-carbon", -1)
    return None


def _beta_methyl_feature(mol: Chem.Mol, site: NitrosamineSite) -> CPCAFeature | None:
    for alpha_idx in site.alpha_carbon_indices:
        alpha = mol.GetAtomWithIdx(alpha_idx)
        for beta in alpha.GetNeighbors():
            if beta.GetIdx() == site.nitrosamine_n_index or beta.GetAtomicNum() != 6:
                continue
            for neighbor in beta.GetNeighbors():
                if neighbor.GetIdx() == alpha_idx or neighbor.GetAtomicNum() != 6:
                    continue
                if (
                    neighbor.GetHybridization() == Chem.HybridizationType.SP3
                    and neighbor.GetTotalNumHs() == 3
                ):
                    return _feature("methyl_on_beta", "Methyl group bonded to beta-carbon", -1)
    return None


def _score_to_category(score: int) -> int:
    if score >= 4:
        return 4
    if score == 3:
        return 3
    if score == 2:
        return 2
    return 1


def _category_five(site: NitrosamineSite, reason: str) -> CPCASiteResult:
    return CPCASiteResult(
        site=site,
        potency_score=None,
        alpha_hydrogen_score=None,
        deactivating_features=(),
        activating_features=(),
        category=5,
        ai_ng_per_day=CATEGORY_AI_NG_PER_DAY[5],
        rationale=(reason,),
    )


def _assess_site(mol: Chem.Mol, site: NitrosamineSite) -> CPCASiteResult:
    if not site.cpca_applicable:
        return CPCASiteResult(
            site=site,
            potency_score=None,
            alpha_hydrogen_score=None,
            deactivating_features=(),
            activating_features=(),
            category=None,
            ai_ng_per_day=None,
            rationale=("FDA CPCA scope exclusion applies",),
            warnings=site.exclusion_reasons,
        )

    alpha_h = site.alpha_hydrogens
    if max(alpha_h) == 0:
        return _category_five(site, "No alpha-hydrogens are present")
    if max(alpha_h) <= 1:
        return _category_five(site, "Neither side has more than one alpha-hydrogen")
    if _has_tertiary_alpha_carbon(mol, site):
        return _category_five(site, "A tertiary alpha-carbon is present")

    alpha_score = _alpha_hydrogen_score(mol, site)
    if alpha_score is None:
        return CPCASiteResult(
            site=site,
            potency_score=None,
            alpha_hydrogen_score=None,
            deactivating_features=(),
            activating_features=(),
            category=None,
            ai_ng_per_day=None,
            rationale=(
                f"Alpha-hydrogen pattern {tuple(sorted(alpha_h))} is not defined in FDA Table A",
            ),
            warnings=("cpca_alpha_h_pattern_not_defined",),
        )

    deactivating: list[CPCAFeature] = []
    activating: list[CPCAFeature] = []

    if mol.HasSubstructMatch(_CARBOXYLIC_ACID):
        deactivating.append(
            _feature(
                "carboxylic_acid",
                "Carboxylic acid group anywhere in the molecule",
                3,
            )
        )

    for candidate in (
        _ring_feature(mol, site),
        _long_chain_feature(mol, site),
        _ewg_feature(mol, site),
        _beta_hydroxyl_feature(mol, site),
    ):
        if candidate is not None:
            deactivating.append(candidate)

    for candidate in (
        _aryl_alpha_feature(mol, site),
        _beta_methyl_feature(mol, site),
    ):
        if candidate is not None:
            activating.append(candidate)

    score = (
        alpha_score
        + sum(feature.score for feature in deactivating)
        + sum(feature.score for feature in activating)
    )
    category = _score_to_category(score)

    return CPCASiteResult(
        site=site,
        potency_score=score,
        alpha_hydrogen_score=alpha_score,
        deactivating_features=tuple(deactivating),
        activating_features=tuple(activating),
        category=category,
        ai_ng_per_day=CATEGORY_AI_NG_PER_DAY[category],
        rationale=(
            f"alpha-H score={alpha_score}",
            f"deactivating score={sum(feature.score for feature in deactivating)}",
            f"activating score={sum(feature.score for feature in activating)}",
        ),
        warnings=("ewg_detection_uses_a_finite_open_pattern_catalog",),
    )


def assess_cpca(smiles: str) -> CPCAAssessment:
    """Apply FDA CPCA structural rules to each detected N-nitrosamine site.

    This is an auditable structural screening implementation, not a regulatory
    determination. Current FDA tables and compound-specific/read-across limits
    take precedence over a calculated CPCA value.
    """

    mol = mol_from_smiles(smiles)
    sites = find_nitrosamine_sites(smiles)
    results = tuple(_assess_site(mol, site) for site in sites)
    warnings: list[str] = []

    if not results:
        warnings.append("no_r2n_n_eq_o_site_detected")
        return CPCAAssessment(
            smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            sites=(),
            overall_category=None,
            overall_ai_ng_per_day=None,
            warnings=tuple(warnings),
            source=FDA_CPCA_SOURCE,
            source_version=FDA_CPCA_SOURCE_VERSION,
        )

    if len(results) > 2:
        warnings.append(
            "fda_recommends_contacting_the_agency_for_more_than_two_n_nitroso_groups"
        )

    categorized = [
        result
        for result in results
        if result.category is not None and result.ai_ng_per_day is not None
    ]
    if len(categorized) != len(results) or len(results) > 2:
        overall_category = None
        overall_ai = None
        warnings.append("overall_cpca_not_assigned")
    else:
        most_stringent = min(
            categorized,
            key=lambda result: (result.ai_ng_per_day, result.category or 99),
        )
        overall_category = most_stringent.category
        overall_ai = most_stringent.ai_ng_per_day

    warnings.append("verify_current_fda_table_before_regulatory_use")
    return CPCAAssessment(
        smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        sites=results,
        overall_category=overall_category,
        overall_ai_ng_per_day=overall_ai,
        warnings=tuple(dict.fromkeys(warnings)),
        source=FDA_CPCA_SOURCE,
        source_version=FDA_CPCA_SOURCE_VERSION,
    )
