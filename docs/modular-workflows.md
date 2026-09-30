# Mezclas, condiciones y cálculos modulares

## Qué se incorporó

| Capacidad | API | Interpretación |
|---|---|---|
| Mezclas y rutas multietapa | `engine.predict_mixture()` | Productos compatibles con reglas suministradas |
| Comparación de escenarios | `compare_scenarios()` | Mismo conjunto de reglas y límites, distintas condiciones |
| Redes de reacción | `ReactionNetwork.to_graph_dict()` | Grafo bipartito con correaccionantes y coproductos |
| Exportación | `to_dict()`, `to_rows()` | JSON reproducible o registros para pandas/CSV |
| Cobertura estructural | `SimilarityDomain` | Vecinos Morgan/Tanimoto; no confianza calibrada |
| Rendimiento con modelo externo | `YieldPredictor` | Predicción con escala y procedencia explícitas |
| Especiación de nitrosación | `nitrosation_speciation()` | Fracciones de equilibrio con pKa suministrados |
| Presupuesto de nitrito | `nitrite_limited_bound()` | Cota estequiométrica de mononitrosación |
| Exposición térmica | `first_order_exposure()` | Pérdida de precursor con parámetros cinéticos suministrados |

## Mezclas y redes

```python
from chem_predict import Conditions, DegradationEngine, Rule, RuleRegistry

# Regla de demostración de la API, no regla farmacéutica curada.
registry = RuleRegistry([
    Rule(
        id="demo_substitution", name="Sustitución de demostración",
        reaction_smarts="[C:1][Br:2].[O;H1;-1:3]>>[C:1][O;+0:3].[Br-:2]",
        source="Demostración",
    )
])
engine = DegradationEngine(registry)
network = engine.predict_mixture(
    {"sustrato": "CCBr", "reactivo": "[OH-]"},
    Conditions(ph=8, temperature_c=25, duration_h=2),
    max_depth=2,
    max_species=500,
    max_steps=1000,
    max_combinations_per_rule=1000,
    max_products_per_rule=100,
)
print(network.to_rows())
print(network.truncated, network.warnings)

import json
from pathlib import Path
Path("network.json").write_text(json.dumps(network.to_dict(), indent=2), encoding="utf-8")
```

Cada componente es una especie disponible; un SMILES con puntos sigue siendo
una entrada molecular con fragmentos, no una lista implícita de reactivos.
Para reglas con varios reactivos, suministrar cada especie en su propia entrada.
Los nombres originales se conservan aunque dos entradas representen la misma estructura.
Las reglas se filtran por condiciones y cada combinación se ejecuta una sola vez.
Los productos quedan disponibles a partir de la siguiente ronda. Los ciclos se
conservan como pasos, pero no generan una expansión infinita.

La red representa disponibilidad cualitativa. No consume reactivos, calcula
concentraciones ni aplica leyes de velocidad. Una especie puede ocupar varios
slots de una regla. Todos los coproductos pueden participar en pasos posteriores.
La prioridad ordena la ejecución y la selección si se alcanza un límite; no es
una probabilidad. En búsquedas truncadas ese orden introduce un sesgo de selección.

`max_depth` define el horizonte solicitado. `truncated` indica que se alcanzó
otro límite de recursos; los límites y condiciones quedan guardados en el JSON.
El límite por regla cuenta combinaciones a lo largo de toda la búsqueda.
El límite de productos restringe la generación cruda de RDKit, antes de
deduplicar: una molécula simétrica puede agotarlo aunque tenga pocos productos únicos.
`on_error="collect"` registra fallos de sanitización por combinación y sigue con
las demás reglas. El modo predeterminado es `"raise"`; errores de reglas o
entradas inválidas siempre deben corregirse.

El método antiguo `predict()` mantiene su formato de salida y ahora omite las
reglas con aridad distinta de uno. Utilizar `predict_mixture()` para ellas.

```python
from chem_predict.degradation import compare_scenarios

results = compare_scenarios(registry, ["CCBr", "[OH-]"], {
    "acido": Conditions(ph=3, temperature_c=25),
    "neutro": Conditions(ph=7, temperature_c=25),
    "caliente": Conditions(ph=7, temperature_c=60),
}, max_depth=2)
```

Los escenarios solo difieren si las reglas tienen ventanas aplicables. Cambiar
la temperatura no altera por sí solo la velocidad ni la prioridad de una regla.
`ConditionWindow` también acepta `duration_min_h` y `duration_max_h`.
Un valor desconocido no satisface una restricción que requiera ese dato.

## Dominio estructural

```python
from chem_predict.applicability import SimilarityDomain

domain = SimilarityDomain(["CC=O>>CCO", "CCC=O>>CCCO"],
                          mode="reaction", threshold=0.6)
assessment = domain.assess("CCCC=O>>CCCCO", k=2)
print(assessment.in_domain, assessment.nearest_similarity, assessment.neighbors)
```

En modo reacción se concatenan fingerprints Morgan de reactivos, agentes y
productos. Se preserva el rol de cada lado y se ignoran los números de mapeo;
no se normalizan sales o tautómeros. La multiplicidad de moléculas idénticas
no queda representada en estos fingerprints binarios. En modo `molecule`
solo se representa una estructura. El umbral predeterminado 0.5 es una decisión
de búsqueda, no un límite validado. Ajustarlo con datos de validación independientes
y analizar error frente a cobertura. No incluir el conjunto de prueba como referencia.
No se representan pH, temperatura, escala o tiempo de reacción.

## Modelos de rendimiento

```python
from chem_predict.yields import YieldPredictor

# `model` ya cargado en su entorno compatible, con predict(list[str])
# que devuelve (predicciones, salidas_crudas), como rxnfp/SimpleTransformers.
# mean/std deben ser los usados al entrenar ese checkpoint, en unidades %.
predictor = YieldPredictor.from_rxn_yields_model(
    model,
    model_id="nombre/checkpoint-version",
    output_scale="standardized",
    training_mean=training_mean,
    training_std=training_std,
    domain=domain,
)
predictions = predictor.predict(["CC=O>>CCO"])
```

También se puede pasar una función `predict_batch(list[str]) -> list[float]`
al constructor. Las escalas disponibles son `percent`, `fraction` y
`standardized`. El adaptador no cambia el texto, tokeniza ni instala modelos;
la preparación del input debe coincidir con el entrenamiento. Rendimientos
fuera de 0–100 se conservan con una advertencia, evitando ocultarlos con clipping.
Una predicción de rendimiento sintético no equivale a conversión de degradación,
probabilidad de formación ni proporción en una mezcla.

## Nitrosación y exposición térmica

```python
from chem_predict.nitrosamines import (
    excipient_nitrite_mass, nitrite_limited_bound, nitrosation_speciation,
)

nitrite_mg = excipient_nitrite_mass(excipient_mass_mg=100, nitrite_ppm=10)
bound = nitrite_limited_bound("CCNCC", amine_mass_mg=100, nitrite_mass_mg=nitrite_mg)
print(bound.product_bounds_mg)

# Valores ilustrativos; sustituir por pKa adecuados al sistema.
speciation = nitrosation_speciation(3.2, amine_pka=10, nitrous_acid_pka=3.2)
print(speciation.nitrous_acid_fraction, speciation.free_amine_fraction)
```

La masa de nitrito se expresa como ion NO2−, no como NaNO2. La masa del precursor
debe corresponder a la estructura suministrada, no a una sal o mezcla distinta.
Las cotas de diferentes productos posicionales comparten un mismo presupuesto:
no deben sumarse. El cálculo solo considera productos de mononitrosación directa
de aminas secundarias neutras no amídicas; las nitrosamidas y rutas de aminas
terciarias requieren modelos separados. No se neutralizan sales automáticamente.

El screening conserva centros amídicos anotados pero ya no los presenta como
precursores de la vía de aminas ni enumera sus productos en esa función.
Un nivel positivo de nitrito en excipientes respalda la presencia de fuente;
si contradice `nitrite_present=False`, se rechaza la entrada.

```python
from chem_predict.degradation import ThermalStage, first_order_exposure

# Parámetros ilustrativos, no calibrados para un compuesto real.
loss = first_order_exposure(
    [ThermalStage(25, 24), ThermalStage(40, 2)],
    k_reference_per_h=0.001, reference_temperature_c=25,
    activation_energy_kj_mol=50,
)
print(loss.remaining_fraction, loss.converted_fraction)
```

Cinética de primer orden y Arrhenius, sin cambios de mecanismo, pH o medio.
La conversión calculada es desaparición del precursor. Para obtener rendimiento
de productos hacen falta cinética de ramificación y balance de especies.

## CLI

```bash
chem-predict mixture rules.json 'CCBr' '[OH-]' --depth 2 --ph 8 --output network.json
chem-predict domain references.txt 'CC=O>>CCO' --mode reaction --threshold 0.6
chem-predict nitrosation 'CCNCC' --nitrite --ph 3.5
```

Los errores de entrada retornan código 2 con un mensaje breve. La expansión
truncada retorna JSON con `truncated=true` y advertencias; comprobarlas antes
de interpretar el resultado como una búsqueda completa.

## Evaluación de las fuentes solicitadas

| Fuente | Resultado de la revisión |
|---|---|
| [SynKit](https://github.com/TieuLongPhan/SynKit), [artículo](https://pubs.acs.org/jcisd8/article/65/24/13012/3688699/SynKit-A-Graph-Based-Python-Framework-for-Rule) | Ya existía adaptador ITS 1.6.x. Se agregó validación de mapeo completo y exportación de cambios compatible con JSON. |
| [ResolveMass](https://resolvemass.ca/process-parameters-nitrosamine-formation/) | Sus tablas no suministran protocolo/datos suficientes para parametrizar un modelo transferible. Se incorporaron cálculos explícitos con entradas y supuestos conocidos. |
| [Reiss et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11921035/), DOI 10.1021/acs.chemrestox.4c00435 | Describe modelos de susceptibilidad a nitrosación, distintos de CPCA. No se reproduce su GNN ni se inventan categorías a partir del screening actual. |
| [reaction-network](https://github.com/materialsproject/reaction-network) | Dirigido a síntesis inorgánica con entradas termodinámicas. Se implementó una red molecular bipartita propia sin incorporar esa dependencia ni atribuir energías a las rutas. |
| [SynAD](https://github.com/deepsynthesis/synad) | Su código contiene métodos de dominio con matrices de descriptores. El checkout revisado fija NumPy 1.22.4 y XGBoost 1.5.2; no se agrega al entorno >=3.11. La cobertura Morgan local es una alternativa básica, no una reproducción de SynAD. |
| [rxn_yields](https://github.com/rxn4chemistry/rxn_yields) | Se comprobó la estandarización por media/desvío en sus scripts y el contrato predict de rxnfp/SimpleTransformers. El adaptador exige los parámetros del checkpoint y mantiene la dependencia externa opcional. |

Revisión de APIs: 2026-09-30. Morgan se basa en la [API oficial de RDKit](https://www.rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html).
Los ejemplos de reglas y parámetros son demostraciones de funcionamiento;
la nueva infraestructura no incorpora un catálogo validado de degradación.
