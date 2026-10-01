# Notebook completo: reacciones y mezclas

[Abrir el notebook](../notebooks/ChemPredict_reacciones_completo.ipynb) · [Instalación](../README.md#instalación) · [Documentación](README.md)

## Ejecutar en el checkout local

Desde la raíz del repositorio, con el entorno virtual activado:

```bash
python -m pip install -e ".[notebook]"
python -m jupyterlab notebooks/ChemPredict_reacciones_completo.ipynb
```

Seleccionar el kernel de ese entorno, mantener `INSTALL_DEPENDENCIES = False` y ejecutar **Run All**. El notebook incluye salidas de ejemplo para inspeccionarlo sin ejecutar. Una ejecución nueva reemplaza esas salidas.

En Colab u otro entorno nuevo, cargar el archivo `.ipynb`, cambiar `INSTALL_DEPENDENCIES = True` y ejecutar desde la primera celda. Esa opción instala el commit `3ca742e8a7cd68e94a23cfe45677cb045eee1d38`, fijado para reproducir el ejemplo. Si se quiere evaluar otro commit, editar `CHEM_PREDICT_REF`. Reiniciar el kernel si ya se habían importado otras versiones de las dependencias.

## Adaptar el análisis

| Entrada | Cómo usarla |
| --- | --- |
| `MIXTURE` | Diccionario nombre → SMILES de las especies disponibles, incluidos correaccionantes |
| `SCENARIOS` | Estrés, pH, temperatura y duración de cada escenario |
| `SCENARIO_TO_INSPECT` | Escenario usado en la inspección detallada |
| `RULES_JSON` | Ruta a reglas propias; `None` activa las cuatro reglas demostrativas |
| `NETWORK_LIMITS` | Topes de profundidad, especies, pasos, combinaciones y productos |
| `OUTPUT_DIR` | Directorio relativo para los archivos exportados |
| `MAX_STRUCTURES_TO_DRAW` | Número máximo de estructuras mostradas |

El motor utiliza estructuras explícitas: no convierte automáticamente nitrito en ácido nitroso ni ajusta protonación por pH. En el ejemplo se aporta `O=NO` como correaccionante de nitrosación.

Las mismas cuatro reglas están disponibles en [examples/rules_mixture_demo.json](../examples/rules_mixture_demo.json) para usarlas desde la CLI. Sus ventanas de condiciones sirven para ilustrar el filtrado; no son parámetros ajustados a datos experimentales.

## Recorrido y resultados

1. Revisar estructuras, fórmulas, masas y cargas de la mezcla.
2. Inspeccionar las reglas SMARTS y comprobar balance de átomos y carga.
3. Ejecutar un control de hidrólisis y visualizar una reacción con mapeo explícito.
4. Comparar cuatro escenarios y construir una red multietapa con coproductos.
5. Ejecutar controles negativos y revisar la similitud con un panel de referencia.
6. Explorar nitrosación, especiación, cota por nitrito y exposición térmica con parámetros ilustrativos.
7. Usar SynKit o un modelo de rendimiento si están disponibles y exportar resultados.

Con la configuración incluida, los escenarios producen estos conteos de referencia:

| Escenario | Pasos | Productos nuevos distintos |
| --- | ---: | ---: |
| Control | 0 | 0 |
| Ácido | 2 | 3 |
| Ácido, oxidante y nitrosante | 4 | 5 |
| Ácido caliente | 5 | 6 |

Estos conteos describen la enumeración del ejemplo, no tasas de reacción. La red conserva pasos y coproductos; no calcula consumo de reactivos ni balances de concentración. Revisar las advertencias y la señal de truncamiento.

## Archivos exportados

Se escriben tablas CSV, redes y reglas JSON, estructuras y gráficos SVG/PNG, y un manifiesto con versiones y configuración en `chem_predict_resultados/`, relativo al directorio de ejecución. Este directorio está ignorado por git. Cambiar `OUTPUT_DIR` para conservar ejecuciones separadas.

## Integraciones opcionales

SynKit se instala con `python -m pip install -e ".[synkit]"`; el notebook omite esa sección si falta la dependencia. El adaptador necesita reacciones completamente mapeadas. No asigna mapas automáticamente.

El ejemplo deja `yield_model = None`: para obtener rendimientos se debe aportar un modelo compatible y sus datos de calibración. No se generan rendimientos ficticios en su ausencia.

## Verificación y problemas frecuentes

```bash
python scripts/validate_notebook.py
```

El script valida el formato y ejecuta todas las celdas con un kernel nuevo en un directorio temporal. Requiere el extra `notebook` y no modifica el archivo versionado.

- **No se importa Chem-Predict:** seleccionar el kernel correcto e instalar con el mismo intérprete del kernel.
- **No aparecen productos:** revisar SMILES, formas de protonación, correaccionantes, etiquetas de estrés y ventanas de condiciones.
- **La red se trunca:** revisar `warnings` y ajustar los límites con cautela; las combinaciones crecen con el tamaño de la mezcla.
- **SMILES inválido:** corregir la estructura antes de continuar; no sustituirla silenciosamente por otra.
- **Una aserción falla al usar reglas propias:** los controles y conteos están ligados a la mezcla y reglas demostrativas; adaptar también esas comprobaciones al caso nuevo.
