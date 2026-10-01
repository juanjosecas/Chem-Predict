# Referencia de la CLI

[Inicio](README.md) · [Reglas](rules-format.md) · [Flujos Python](modular-workflows.md)

Tras instalar el paquete, usar `chem-predict --help` y `chem-predict COMANDO --help`. Ejecutar los ejemplos desde la raíz del repositorio y encerrar los SMILES/SMARTS entre comillas.

| Comando | Entrada | Salida |
| --- | --- | --- |
| `normalize` | Un SMILES | SMILES canónico |
| `apply` | SMARTS de reacción y reactivos separados | Productos unidos por puntos, una salida por línea |
| `predict` | Archivo de reglas y una molécula | TSV: regla, prioridad, productos |
| `mixture` | Archivo de reglas y especies disponibles | JSON de red multietapa |
| `domain` | Archivo de referencias y consulta | JSON de similitud y vecinos |
| `nitrosation` | SMILES y contexto opcional | JSON de evaluación y productos directos |

## Ejemplos mínimos

```bash
chem-predict normalize "OCC"
chem-predict apply "[C:1][Br:2].[O;H1;-1:3]>>[C:1][O;+0:3].[Br-:2]" "CCBr" "[OH-]"
chem-predict mixture examples/rules_mixture_demo.json "CCOC(C)=O" "O" \
  --stress acid --ph 3 --temperature-c 25 --duration-h 24 --depth 2
chem-predict nitrosation "CNC" --ph 3 --nitrite
```

`predict` procesa reglas de un solo reactivo; para las reglas multirreactivo del ejemplo usar `mixture`.

## Condiciones de `predict` y `mixture`

`--stress` puede repetirse. Las opciones numéricas son `--ph`, `--temperature-c` y `--duration-h`. Oxígeno y luz admiten `--oxygen` / `--no-oxygen` y `--light` / `--no-light`. Omitir una opción significa condición desconocida; una regla que la exige no se activa.

Las etiquetas de estrés se seleccionan explícitamente. `--ph 3` por sí solo no activa reglas restringidas a `acid`. Las condiciones no incorporan agua, oxidantes u otros reactivos al listado de especies.

## Límites de enumeración

| Opción de `mixture` | Valor por defecto |
| --- | ---: |
| `--depth` | 2 |
| `--max-species` | 500 |
| `--max-steps` | 1000 |
| `--max-combinations-per-rule` | 1000 |
| `--max-products-per-rule` | 100 |

`--output ruta.json` escribe el resultado en un archivo; el directorio padre debe existir. `--collect-errors` registra fallos de aplicación en `warnings` y continúa; sin esa opción los errores interrumpen el cálculo. Una red truncada puede terminar con código 0: revisar siempre `truncated` y `warnings` en el JSON.

En `apply`, `--max-products` vale 1000 por defecto; 0 elimina ese límite de RDKit. En `predict`, `--max-products-per-rule` vale 1000. Los límites de productos se aplican a las salidas de RDKit antes de la deduplicación.

## Cobertura estructural

Crear `referencias.smi` en UTF-8, con un SMILES por línea; se ignoran líneas vacías y comentarios que empiezan con `#`.

```bash
chem-predict domain referencias.smi "CCOC(C)=O" --threshold 0.5 --neighbors 3
```

Para referencias de reacción usar `--mode reaction`, con una reacción SMILES por línea y una consulta del mismo tipo. El umbral de similitud es configurable y no equivale a incertidumbre calibrada.

La CLI devuelve código 2 para errores de entrada habituales y muestra un mensaje en stderr.
