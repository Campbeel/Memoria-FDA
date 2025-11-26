# Memoria FDA + IBEX

Guía rápida para compilar y ejecutar las variantes de Feasible Diving (FDA) integradas con IBEX.

## Prerrequisitos
- IBEX ya compilado en `~/ibex-lib` (este repo).
- CMake ≥ 3.10
- Compilador C++17
- Dependencias de IBEX (gaol/filib, Clp/CoinUtils instalados con IBEX).

## Compilación
Desde `memoria_cofigo`:
```bash
cmake -S . -B build
cmake --build build -j4
```
Esto genera los binarios:
- `fda_run_base`
- `fda_run_depth_k`
- `fda_run_depth_k_rand`
- `fda_run_vol_k`
- `fda_run_vol_k_rand`
- `fda_bridge_ibex` (ejecutor multiparámetro)

## Ejecución de variantes individuales
Formato: `./<binario> ruta/al/problema.bch [num_runs]`

Ejemplos:
```bash
cd memoria_cofigo/build

# Base
./fda_run_base ../casos/medium/alkyl.bch 3

# Profundidad + temperatura determinista
./fda_run_depth_k ../casos/medium/alkyl.bch 3

# Profundidad + temperatura con ruido
./fda_run_depth_k_rand ../casos/medium/alkyl.bch 3

# Volumen + temperatura determinista
./fda_run_vol_k ../casos/medium/alkyl.bch 3

# Volumen + temperatura con ruido
./fda_run_vol_k_rand ../casos/medium/alkyl.bch 3
```

La salida es CSV: `run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal`.

## Ejecución con bridge (todas las variantes en un binario)
`fda_bridge_ibex` permite elegir variante por nombre:
```bash
./fda_bridge_ibex ../casos/medium/alkyl.bch depth_k 5
```
Variantes disponibles: `base`, `base_bb`, `depth_k`, `depth_k_bb`, `depth_k_rand`, `depth_k_rand_bb`, `vol_k`, `vol_k_bb`, `vol_k_rand`, `vol_k_rand_bb`.

## Parámetros internos (hardcodeados)
- `eps_box = 1e-9`
- `max_iters = 100000`
- Profundidad/temperatura: `T0=10`, `k=0.5`, `d_max=100`
- Volumen: `eps_V=1e-4`, `beta=0.1`, `depth_vol>=5` antes de cortar por volumen

Ajusta en los `.cpp` si necesitas otros valores por defecto.
