# Memoria FDA + IBEX

Guía rápida para compilar y ejecutar Ibex con los modos FD (`--fd-mode`) integrados en `ibex_opt_full`.

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
- `ibex_opt_full` (original de Ibex con soporte `--fd-mode`)
- `ibex_opt_base` (optimizador base sin LP/XTaylor; útil si `ibexopt` se cae con CLP)
- `ibex_menu` (menú interactivo que llama a `ibex_opt_full` con el modo elegido)

### Rutas estándar (ajusta en tu entorno)
Por defecto el `CMakeLists.txt` asume:
- `IBEX_ROOT=~/ibex-lib`
- `IBEX_BUILD_DIR=~/ibex-lib/build-soplex`
- `SOPLEX_LIB=/usr/lib/x86_64-linux-gnu/libsoplex.a`

En otro host, sobreescribe al configurar:
```bash
cmake -S . -B build \
  -DIBEX_ROOT=/ruta/a/ibex \
  -DIBEX_BUILD_DIR=/ruta/a/ibex/build-soplex \
  -DSOPLEX_LIB=/ruta/a/libsoplex.a
```

### Cómo obtener SoPlex (recomendado) dentro de tu árbol de Ibex
1) Descarga SoPlex (ZIB Academic License) desde https://soplex.zib.de/.
2) Compílalo y obtén `libsoplex.a`. Ejemplo rápido:
```bash
tar xf soplex-*.tar.gz
cd soplex-*
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
cp build/libsoplex.a /ruta/a/ibex/build-soplex/lib/   # o donde prefieras
```
3) Reconfigura Ibex con SoPlex:
```bash
cd /ruta/a/ibex
cmake -S . -B build-soplex -DLP_LIB=soplex -DUSE_SOPLEX=ON
cmake --build build-soplex -j4
```
4) Reconfigura `memoria_cofigo` apuntando a esas rutas (ver comando arriba).

Si no quieres/puedes usar SoPlex, puedes compilar Ibex con CLP y pasar `-DSOPLEX_LIB=` para omitirlo, pero perderás la mitigación del bug de CLP.

## Ejecución directa (ibex_opt_full)
Formato general: `./ibex_opt_full problema.bch [opciones]`

Modos FD disponibles:
- Base (sin flag): `./ibex_opt_full ../casos/medium/ex6_2_14.bch`
- Profundidad determinista: `./ibex_opt_full ../casos/medium/ex6_2_14.bch --fd-mode=depth_k`
- Profundidad aleatoria: `./ibex_opt_full ../casos/medium/ex6_2_14.bch --fd-mode=depth_k_rand`
- Volumen determinista: `./ibex_opt_full ../casos/medium/ex6_2_14.bch --fd-mode=vol_k`
- Volumen aleatoria: `./ibex_opt_full ../casos/medium/ex6_2_14.bch --fd-mode=vol_k_rand`

La salida es el reporte estándar de IbexOpt.

### Optimización base (sin CLP)
`ibex_opt_base` evita el loup XTaylor/CLP y suele ser más robusto:
```bash
./ibex_opt_base ../casos/medium/alkyl.bch [rel_eps_f] [abs_eps_f]
```
Los eps son opcionales (por defecto usa los de Ibex). Para problemas grandes puede tardar más que `ibexopt`, pero no se cae por el bug de CLP.

## Menú interactivo
`ibex_menu` permite elegir variante y problema (.bch) sin recordar las flags:
```bash
./ibex_menu
```
Las opciones del menú disparan `ibex_opt_full` con el `--fd-mode` correspondiente.
