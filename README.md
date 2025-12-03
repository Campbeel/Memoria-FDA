# Memoria FDA – Ejecutables y datos listos para mover

Carpeta autocontenida con los ejecutables variantes (`base`, `depth_k`, `vol_k`, `*_rand`), los buffers térmicos y todos los CSV de resultados agregados.

## Requisitos
- CMake ≥ 3.10
- Un build de Ibex (ruta en `CMakeLists.txt` por defecto `IBEX_ROOT="/home/benjamin-mu-oz/ibex-lib"`), con SoPlex disponible.
- Compilador C++17

> Nota: las clases `ibex_Cell.{h,cpp}` con temperatura están copiadas en `memoria_cofigo/ibex_cell/`. Si usas otra instalación de Ibex, sobrescribe su `src/cell/` con esos dos archivos o usa directamente este repo como `IBEX_ROOT`.

## Compilar
Dentro de `memoria_cofigo/`:
```bash
mkdir -p build
cd build
cmake ..
make -j
```
Los binarios quedan en `memoria_cofigo/build/`:
- `ibex_opt_full` (base)
- `ibex_opt_full_depth_k`, `ibex_opt_full_depth_k_rand`
- `ibex_opt_full_vol_k`, `ibex_opt_full_vol_k_rand`

## Ejecutar
```bash
./ibex_opt_full_vol_k path/al/archivo.bch --random-seed 123
./ibex_opt_full_depth_k path/al/archivo.bch
```
Parámetros útiles:
- `--quiet` para suprimir salidas.
- `--output resultado.cov` para guardar el COV.
- `--random-seed N` para reproducibilidad.

## Datos de resultados
- Individuales: `memoria_cofigo/build/results/results_*.csv`
- Agregados: `memoria_cofigo/build/results/super/results_*_all.csv`

## Transportar la carpeta
Desde la raíz del repo:
```bash
tar -czf memoria_cofigo_ready.tgz memoria_cofigo
```
Si necesitas llevar también las clases de temperatura de Ibex:
```bash
tar -czf ibex_temp_support.tgz memoria_cofigo/ibex_cell/ibex_Cell.h memoria_cofigo/ibex_cell/ibex_Cell.cpp
```
