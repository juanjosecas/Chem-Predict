from __future__ import annotations

import argparse

from chem_predict.chemistry import apply_reaction, canonicalize_smiles
from chem_predict.core import Conditions, StressType
from chem_predict.degradation import DegradationEngine
from chem_predict.rules import RuleRegistry


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
    predict.add_argument("--oxygen", action=argparse.BooleanOptionalAction, default=None)
    predict.add_argument("--light", action=argparse.BooleanOptionalAction, default=None)
    predict.add_argument("--max-products-per-rule", type=int, default=1000)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

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


if __name__ == "__main__":
    raise SystemExit(main())
