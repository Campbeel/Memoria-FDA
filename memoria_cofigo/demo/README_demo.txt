Demo fuera del árbol principal
==============================

Esta carpeta permite compilar los ejecutables de prueba sin tocar el CMake del
proyecto grande. Usa las mismas fuentes y las mismas rutas de IBEX.

Cómo compilar:
1) Desde `memoria_cofigo/demo`:
   cmake -S . -B build
   cmake --build build --target mini_ibexsolve fda_bridge_ibex_demo -j4

   Si IBEX no está en /home/benjamin-mu-oz/ibex-lib, pasa IBEX_ROOT:
   cmake -S . -B build -DIBEX_ROOT=/ruta/a/ibex

2) Ejecutables resultantes en `demo/build/`:
   - mini_ibexsolve: solver mínimo estilo ibexsolve.
   - fda_bridge_ibex_demo: las variantes FDA (base, depth, vol, etc.).

Ejemplos:
   ./build/mini_ibexsolve ../casos/hard/chem-1.bch
   ./build/fda_bridge_ibex_demo ../casos/medium/alkyl.bch base_bb 3

Las rutas a .bch son relativas a `memoria_cofigo/demo` (usa ../casos/...).
