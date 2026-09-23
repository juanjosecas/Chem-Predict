from chem_predict.core import Conditions, StressType
from chem_predict.degradation import DegradationEngine
from chem_predict.rules import RuleRegistry

registry = RuleRegistry.from_json("examples/rules_demo.json")
engine = DegradationEngine(registry)

predictions = engine.predict(
    "CC=O",
    Conditions(stresses=frozenset({StressType.REDUCTION})),
)

for prediction in predictions:
    print(prediction.rule_id, prediction.products)
