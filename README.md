# Chem-Predict

Herramientas modulares en Python para explorar reacciones, degradación e impurezas a partir de estructuras y reglas explícitas. Cada producto candidato conserva la regla y los reactivos que lo originaron.

**Empezar con el [notebook completo de reacciones](notebooks/ChemPredict_reacciones_completo.ipynb)**: incluye estructuras, mezclas, cuatro escenarios, red multietapa, controles, gráficos y exportación. La [guía del notebook](docs/notebook.md) explica cómo ejecutarlo y adaptar los ejemplos.

## Capacidades actuales

| Necesidad | Funcionalidad | Alcance |
| --- | --- | --- |
| Normalizar estructuras | SMILES canónico con RDKit | No estandariza automáticamente sales ni estados de protonación |
| Enumerar productos | SMARTS y reglas JSON, uno o varios reactivos | Depende de las reglas suministradas |
| Explorar mezclas | Red multietapa y comparación de escenarios | Límites de profundidad, especies, pasos y combinaciones |
| Filtrar condiciones | Estrés, pH, temperatura, duración, oxígeno y luz | Ventanas de aplicabilidad de reglas |
| Revisar cobertura | Vecinos por similitud molecular o de reacción | Referencias y umbral definidos por el usuario |
| Explorar nitrosación | Precursores, productos directos, especiación y cota por nitrito | Cálculos separados de formación y evaluación de potencia |
| Calcular exposición térmica | Primer orden y Arrhenius por etapas | Requiere parámetros cinéticos aportados por el usuario |
| Estimar rendimientos | Adaptador de modelos externos | No incluye pesos entrenados ni un predictor por defecto |
| Inspeccionar reacciones | Dibujos, centro de reacción e ITS opcional de SynKit | El centro y el ITS necesitan mapeo atómico |
| Revisar alertas medchem | Port parcial RDKit de Lilly y adaptador oficial | El port no cubre todas las reglas originales |

Los resultados del motor son **estructuras candidatas**. La prioridad de una regla, la cantidad de rutas y la similitud no representan probabilidades, concentraciones ni rendimientos. Las reglas de demostración no constituyen una biblioteca experimental validada.

## Instalación

Python **3.11 o superior**. Ejecutar desde una terminal:

```bash
git clone https://github.com/juanjosecas/Chem-Predict.git
cd Chem-Predict
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[notebook]"
```

Para usar solamente la biblioteca y la CLI: `python -m pip install -e .`.
Para contribuir: `python -m pip install -e ".[dev,notebook]"`.
SynKit se instala por separado: `python -m pip install -e ".[synkit]"`.

## Primer cálculo

Este ejemplo usa las reglas demostrativas incluidas y proporciona el agua como correaccionante:

```bash
chem-predict mixture examples/rules_mixture_demo.json "CCOC(C)=O" "O" \
  --stress acid --ph 3 --temperature-c 25 --duration-h 24 --depth 2 \
  --output red_reacciones.json
```

Genera ácido acético y etanol como candidatos de hidrólisis. El JSON conserva reactivos, productos, reglas, condiciones, límites y advertencias. Revisar `truncated` y `warnings` antes de interpretar una red.

```python
from chem_predict import Conditions, DegradationEngine, RuleRegistry, StressType

engine = DegradationEngine(RuleRegistry.from_json("examples/rules_mixture_demo.json"))
network = engine.predict_mixture(
    {"acetato_de_etilo": "CCOC(C)=O", "agua": "O"},
    Conditions(stresses=frozenset({StressType.ACID}), ph=3,
               temperature_c=25, duration_h=24),
    max_depth=2,
)
for row in network.to_rows():
    print(row["rule_id"], row["reactants"], "->", row["products"])
```

El pH no asigna automáticamente la etiqueta `acid`. Oxígeno y luz filtran reglas; no añaden especies a la mezcla. Ver [condiciones y reglas](docs/rules-format.md).

## Notebook y ejemplos

```bash
python -m jupyterlab notebooks/ChemPredict_reacciones_completo.ipynb
```

Seleccionar el kernel del entorno instalado y ejecutar las celdas en orden. Mantener `INSTALL_DEPENDENCIES = False` para usar el checkout local. El notebook también permite instalar una versión fijada desde un entorno nuevo; ver [guía](docs/notebook.md).

El ejemplo breve `python examples/modular_workflows.py` funciona sin notebook, GPU ni modelos descargados.

## Documentación

| Guía | Contenido |
| --- | --- |
| [Índice](docs/README.md) | Rutas para empezar, ampliar reglas y contribuir |
| [Notebook completo](docs/notebook.md) | Ejecución, entradas editables, resultados y controles |
| [CLI](docs/cli.md) | Comandos, opciones y ejemplos reproducibles |
| [Flujos modulares](docs/modular-workflows.md) | Mezclas, escenarios, cobertura, cinética, nitrosación y rendimientos |
| [Formato de reglas](docs/rules-format.md) | Esquema JSON y semántica de condiciones |
| [Fuentes de reglas](docs/rule-sources.md) | Procedencia y curación |
| [Nitrosaminas](docs/nitrosamines.md) | Evaluación estructural y límites del cálculo |
| [Visualización](docs/visualization.md) | Moléculas, reacciones y centro de reacción |
| [Integraciones](docs/integrations.md) | Lilly y SynKit |
| [Arquitectura](docs/architecture.md) | Módulos y límites entre capas |
| [Dependencias](docs/dependency-policy.md) | Política de dependencias |

## Desarrollo y validación

```bash
python -m pip install -e ".[dev,notebook]"
pytest
python examples/modular_workflows.py
python scripts/validate_notebook.py
```

GitHub Actions ejecuta los tests en Python 3.11 y 3.13, comprueba la integración opcional de SynKit y ejecuta el notebook en un directorio temporal. Los tests verifican comportamiento del software; la validación química requiere referencias experimentales y reglas curadas para cada uso.

## Próximas ampliaciones

Ampliar bibliotecas curadas, validar contra datos experimentales, calibrar modelos de rendimiento y conectar nuevos motores mediante adaptadores. `properties/` y `purge/` son puntos de extensión; no ofrecen todavía modelos completos de propiedades o purga.

Licencia: [MIT](LICENSE).
