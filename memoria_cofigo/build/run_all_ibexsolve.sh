#!/usr/bin/env bash
set -e

# Ejecuta ibexsolve sobre los casos de memoria_cofigo/casos
# y guarda toda la salida (incluyendo warnings) en build.
# Uso:
#   ./run_all_ibexsolve.sh                  # corre todos los casos
#   ./run_all_ibexsolve.sh ackley_5 alkyl   # corre solo esos

cd "$(dirname "$0")/../.."  # raíz del repo (ibex-lib)

OUT_DIR="memoria_cofigo/build"

run_case() {
  local case="$1"   # nombre sin extensión, ej: ackley_5
  local bch="memoria_cofigo/casos/${case}.bch"
  local out="${OUT_DIR}/ibexsolve_${case}.txt"

  if [[ ! -f "$bch" ]]; then
    echo ">> [WARN] No existe $bch, se omite."
    return
  fi

  echo "=== ibexsolve: ${case} ==="
  # script captura todo lo que sale por pantalla, incluso si hay segfault al final
  script -c "ibexsolve ${bch}" "${out}"
}

if [[ "$#" -eq 0 ]]; then
  # casos por defecto (mismos que en run_all_fda.sh)
  run_case "ackley_5"
  run_case "alkyl"
  run_case "rastrigin"
  run_case "griewank"
  run_case "schaffer"
else
  for case in "$@"; do
    run_case "$case"
  done
fi
