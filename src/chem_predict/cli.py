from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from chem_predict.chemistry import apply_reaction, canonicalize_smiles
from chem_predict.core import Conditions, StressType
from chem_predict.degradation import DegradationEngine
from chem_predict.rules import RuleRegistry
from chem_predict.applicability import SimilarityDomain
from chem_predict.nitrosamines import NitrosationContext, assess_nitrosation_context, enumerate_secondary_amine_nitrosation_products


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="chem-predict")
    subparsers = parser.add_subparsers(dest="command", required=True)

    normalize = subparsers.add_parser("normalize", help="Canonicalize a SMILES string")
    normalize.add_argument("smiles")

    apply = subparsers.add_parser("apply", help="Apply reaction SMARTS")
    apply.add_argument("reaction_smarts")
    apply.add_argument("reactants", nargs="+")
    apply.add_argument("--max-products", type=int, default=1000)

    predict = subparsers.add_parser("predict", help="Run rules from a JSON registry")
    predict.add_argument("rules_json")
    predict.add_argument("smiles")
    predict.add_argument("--stress", action="append", choices=[item.value for item in StressType])
    predict.add_argument("--ph", type=float)
    predict.add_argument("--temperature-c", type=float)
    predict.add_argument("--duration-h", type=float)
    predict.add_argument("--oxygen", action=argparse.BooleanOptionalAction, default=None)
    predict.add_argument("--light", action=argparse.BooleanOptionalAction, default=None)
    predict.add_argument("--max-products-per-rule", type=int, default=1000)

    mixture = subparsers.add_parser("mixture", help="Enumerate a bounded multi-step reaction network")
    mixture.add_argument("rules_json")
    mixture.add_argument("components", nargs="+", help="SMILES for each available species (quote each one)")
    mixture.add_argument("--depth", type=int, default=2)
    mixture.add_argument("--max-species", type=int, default=500)
    mixture.add_argument("--max-steps", type=int, default=1000)
    mixture.add_argument("--max-combinations-per-rule", type=int, default=1000)
    mixture.add_argument("--max-products-per-rule", type=int, default=100)
    mixture.add_argument("--ph", type=float)
    mixture.add_argument("--temperature-c", type=float)
    mixture.add_argument("--duration-h", type=float)
    mixture.add_argument("--stress", action="append", choices=[item.value for item in StressType])
    mixture.add_argument("--oxygen", action=argparse.BooleanOptionalAction, default=None)
    mixture.add_argument("--light", action=argparse.BooleanOptionalAction, default=None)
    mixture.add_argument("--collect-errors", action="store_true")
    mixture.add_argument("--output", help="Write JSON to this file instead of stdout")

    domain = subparsers.add_parser("domain", help="Check structural coverage against local reference SMILES")
    domain.add_argument("references", help="UTF-8 file, one SMILES or reaction per line")
    domain.add_argument("query")
    domain.add_argument("--mode", choices=["molecule", "reaction"], default="molecule")
    domain.add_argument("--threshold", type=float, default=0.5)
    domain.add_argument("--neighbors", type=int, default=5)

    nitrosation = subparsers.add_parser("nitrosation", help="Screen precursor/context and enumerate direct products")
    nitrosation.add_argument("smiles")
    nitrosation.add_argument("--ph", type=float)
    nitrosation.add_argument("--nitrite", action=argparse.BooleanOptionalAction, default=None)

    return parser


def _main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "mixture":
        conditions = Conditions(
            stresses=frozenset(StressType(value) for value in (args.stress or [])),
            ph=args.ph, temperature_c=args.temperature_c, duration_h=args.duration_h,
            oxygen=args.oxygen, light=args.light)
        network = DegradationEngine(RuleRegistry.from_json(args.rules_json)).predict_mixture(
            args.components, conditions, max_depth=args.depth, max_species=args.max_species,
            max_steps=args.max_steps, max_combinations_per_rule=args.max_combinations_per_rule,
            max_products_per_rule=args.max_products_per_rule,
            on_error="collect" if args.collect_errors else "raise")
        payload = json.dumps(network.to_dict(), indent=2, ensure_ascii=False) + "\n"
        if args.output:
            Path(args.output).write_text(payload, encoding="utf-8")
        else:
            print(payload, end="")
        return 0

    if args.command == "domain":
        references = [s.strip() for s in Path(args.references).read_text(encoding="utf-8").splitlines()
                      if s.strip() and not s.lstrip().startswith("#")]
        result = SimilarityDomain(references, mode=args.mode, threshold=args.threshold).assess(
            args.query, k=args.neighbors)
        print(json.dumps(asdict(result), indent=2))
        return 0

    if args.command == "nitrosation":
        result = assess_nitrosation_context(args.smiles, NitrosationContext(
            nitrite_present=args.nitrite, ph=args.ph))
        print(json.dumps({"assessment": asdict(result),
                          "direct_products": enumerate_secondary_amine_nitrosation_products(args.smiles)}, indent=2))
        return 0

    if args.command == "normalize":
        print(canonicalize_smiles(args.smiles))
        return 0

    if args.command == "apply":
        for outcome in apply_reaction(
            args.reaction_smarts,
            args.reactants,
            max_products=args.max_products,
        ):
            print(".".join(outcome))
        return 0

    if args.command == "predict":
        registry = RuleRegistry.from_json(args.rules_json)
        conditions = Conditions(
            stresses=frozenset(StressType(value) for value in (args.stress or [])),
            ph=args.ph,
            temperature_c=args.temperature_c,
            duration_h=args.duration_h,
            oxygen=args.oxygen,
            light=args.light,
        )
        predictions = DegradationEngine(registry).predict(
            args.smiles,
            conditions,
            max_products_per_rule=args.max_products_per_rule,
        )
        for prediction in predictions:
            products = ".".join(prediction.products)
            print(f"{prediction.rule_id}\t{prediction.priority}\t{products}")
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


def main(argv: list[str] | None = None) -> int:
    try:
        return _main(argv)
    except (ValueError, TypeError, OSError) as exc:
        import sys
        print(f"chem-predict: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
