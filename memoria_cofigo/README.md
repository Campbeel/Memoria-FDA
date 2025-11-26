# Memoria FDA + IBEX

Guia rapida para compilar y ejecutar las variantes de Feasible Diving adaptadas a IBEX. Incluye la descripcion de los ejecutables, parametros por defecto y scripts auxiliares.

## Prerrequisitos
- CMake >= 3.10 y un compilador C++17.
- IBEX compilado localmente. La ruta se configura en `CMakeLists.txt` via `IBEX_ROOT` (por defecto `/home/benjamin-mu-oz/ibex-lib`); ajusta esa variable a tu instalacion.
- WRAPPERS de intervalos ubicados en `../interval_lib_wrapper` (relativos a `memoria_cofigo/`). Si estan en otro sitio, corrige las rutas en `CMakeLists.txt`.
- Coin-OR CLP y CoinUtils (headers en `/usr/include/coin`, libs en `/usr/local/lib`). Vienen al compilar IBEX con LP.
- FILIB si quieres evaluacion con intervalos (ruta usada en `CMakeLists.txt`: `/home/benjamin-mu-oz/filib-src`).
- Por defecto se compila con AddressSanitizer (`-fsanitize=address`) para cazar corrupciones; quitale esas flags en `CMakeLists.txt` si buscas maximo rendimiento.

## Compilacion
Desde la raiz del repo (`Memoria-FDA`):

```bash
cmake -S memoria_cofigo -B memoria_cofigo/build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build memoria_cofigo/build -j
```

El directorio `memoria_cofigo/build` contendra los ejecutables `fda_bridge_ibex`, `fda_experimentos` y `mini_ibexsolve`.

## Ejecutables principales
### fda_bridge_ibex
Puente que toma un problema IBEX (`.bch`) y ejecuta variantes de Feasible Diving; imprime un CSV por stdout.

Uso:
```bash
./fda_bridge_ibex problema.bch variante [num_runs]
```
Variantes disponibles:
- `base`: FD clasico (heuristica sobre valor medio de la funcion).
- `base_bb`: FD embebido en un B&B superficial.
- `depth_k` / `depth_k_bb`: temperatura jerarquica T*k con limite de profundidad.
- `depth_k_rand` / `depth_k_rand_bb`: igual que depth_k pero con factor aleatorio en T.
- `vol_k` / `vol_k_bb`: corte por volumen relativo `tau_V = eps_V * exp(-beta*depth)`.
- `vol_k_rand` / `vol_k_rand_bb`: variante de volumen con ruido en la temperatura y en `tau_V`.

Parametros internos (por defecto en el codigo):
- `eps_box=1e-9`, `max_iters=100000`.
- `T0=10.0`, `k=0.5`, `d_max=100`.
- `eps_V=1e-4`, `beta=0.1`.
- En las variantes `_bb`: `max_bb_nodes` entre 1e6 y 5e6 segun el caso, `max_iters_fdive=100000`.
- Si la caja inicial es no acotada, las variantes `depth_*` hacen fallback a `base`.

Salida CSV (una fila por run):
```
run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal
```
Ejemplo desde `memoria_cofigo/build`:
```bash
./fda_bridge_ibex ../casos/medium/alkyl.bch depth_k 10 > resultados/alkyl_depth_k.csv
```

### fda_experimentos
Motor sintetico (sin IBEX) que replica la logica de las variantes para experimentos rapidos.

Uso:
```bash
./fda_experimentos <variante> <reps> <max_nodes> <T0> <k> [tauV]
```
`variante` en {`base`, `depth_k`, `depth_k_rand`, `vol_k`, `vol_k_rand`}. Genera CSV por stdout con las mismas columnas que `fda_bridge_ibex`. `tauV` es el umbral de volumen relativo para las variantes de volumen (p.ej. 0.10).

### mini_ibexsolve
Wrapper minimo sobre `ibexsolve`/`DefaultOptimizer` para validar un `.bch`.

Uso:
```bash
./mini_ibexsolve problema.bch [--trace]
```
Intenta optimizar si el problema tiene objetivo; si falla, cae a `DefaultSolver`. Muestra estado, celdas y primera solucion encontrada.

## Casos de prueba incluidos
- `casos/medium/`: coleccion de problemas IBEX de dificultad media (`alkyl`, `bearing`, etc).
- `casos/hard/`: instancias mas exigentes (`chem-1`, `ex6_1_3`, ...).
Ejecuta los binarios desde `memoria_cofigo/build` o ajusta las rutas relativas a `../casos/...`.

## Scripts utiles (en memoria_cofigo/build)
- `run_all_fda.sh [case ...]`: corre todas las variantes sobre un subconjunto de casos y fusiona el CSV en `resultados_ibex/<caso>.csv`. Usa rutas `../casos/<case>.bch`.
- `run_ibex_parallel.sh`: lanza en paralelo (hasta `MAX_JOBS`) las combinaciones definidas en el script, guarda CSVs en `resultados/`.
- `run_memoria_tests.sh`: version parametrizable con lista fija de problemas/variantes, corre en paralelo hasta `MAX_PARALLEL`.
- `run_all_ibexsolve.sh [case ...]`: ejecuta `ibexsolve` sobre los casos seleccionados y guarda la salida completa en `memoria_cofigo/build/ibexsolve_<case>.txt`.

## Flujo recomendado
1) Ajusta rutas en `CMakeLists.txt` para IBEX/filib/CLP si difieren de tu entorno.  
2) Compila con CMake (seccion anterior).  
3) Desde `memoria_cofigo/build`, ejecuta `./fda_bridge_ibex ../casos/<ruta>.bch <variante> [runs]` para obtener CSVs.  
4) Para corridas masivas, usa `run_all_fda.sh` o `run_ibex_parallel.sh`.  
5) Si necesitas validar un `.bch` antes de lanzarlo, prueba `./mini_ibexsolve archivo.bch`.
