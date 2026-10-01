# Arquitectura

[Inicio](README.md) · [Flujos y ejemplos API](modular-workflows.md)

Chem-Predict mantiene el núcleo basado en RDKit y conecta herramientas mayores mediante adaptadores opcionales.

| Módulo | Responsabilidad |
| --- | --- |
| `core` | Condiciones, ventanas de aplicación, tipos de estrés y reglas |
| `chemistry` | Lectura/canonicalización, aplicación SMARTS y visualización |
| `rules` | Registro JSON y procedencia de reglas |
| `degradation` | Predicción individual, redes de mezcla, escenarios y cinética térmica |
| `applicability` | Similitud con un conjunto de referencia local |
| `nitrosamines` | Evaluación estructural, contexto, productos, especiación y cotas |
| `yields` | Interfaz para modelos externos de rendimiento |
| `medchem` | Port parcial de Lilly y backend oficial |
| `integrations` | Adaptador opcional SynKit para ITS y cambios de enlaces |
| `properties`, `purge` | Puntos de extensión pendientes |
| `cli` | Entrada de terminal y serialización de resultados |

## Flujo de cálculo

Las entradas se canonicalizan, el registro filtra reglas por condiciones y el motor aplica SMARTS con RDKit. `predict` considera una molécula y reglas de un reactivo; `predict_mixture` combina especies disponibles según los slots de cada regla y añade productos para las siguientes rondas.

La red registra entradas, profundidad mínima por especie, pasos con regla/reactivos/productos, condiciones, límites, advertencias y truncamiento. `to_dict()` exporta JSON, `to_rows()` ofrece filas de pasos y `to_graph_dict()` representa el grafo bipartito de especies y reacciones. Los slots de coproductos se conservan incluso si dos tienen el mismo SMILES.

No hay consumo de especies, concentraciones, generación automática de correaccionantes ni balance cinético de la mezcla. Los ciclos pueden conservarse como pasos sin provocar una enumeración infinita. La profundidad y los topes controlan la expansión combinatoria.

## Límites entre módulos

Las condiciones restringen reglas; el cálculo térmico es una utilidad independiente que exige parámetros cinéticos. La cobertura estructural tampoco modifica automáticamente los productos. Un modelo de rendimiento debe suministrarse de forma explícita: el adaptador no entrena ni descarga pesos por defecto.

SynKit se importa en la frontera de integración, requiere mapeo positivo y único en cada lado y no modifica el motor SMARTS. La visualización de centros también exige mapas; la red no los genera automáticamente.

Al agregar una función, conservar la procedencia de resultados, validar entradas en la frontera y mantener las dependencias grandes fuera del núcleo. Consultar la [política de dependencias](dependency-policy.md).
