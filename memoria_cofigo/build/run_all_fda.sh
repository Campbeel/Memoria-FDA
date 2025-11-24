#!/usr/bin/env bash
set -e

# Ejecuta fda_bridge_ibex sobre los casos de memoria_cofigo/casos
# Uso:
#   ./run_all_fda.sh                # corre todos los casos
#   ./run_all_fda.sh ackley_5 alkyl # corre solo esos

cd "$(dirname "$0")"  # memoria_cofigo/build

OUT_DIR="resultados_ibex"
mkdir -p "$OUT_DIR"

run_case() {
  local case="$1"   # nombre con ruta relativa dentro de casos, ej: medium/alkyl
  local bch="../casos/${case}.bch"
  local slug="${case//\//_}"
  local merged="${OUT_DIR}/${slug}.csv"

  if [[ ! -f "$bch" ]]; then
    echo ">> [WARN] No existe $bch, se omite."
    return
  fi

  echo "=== FDA: ${case} ==="
  # generamos un solo CSV por caso, con todas las variantes apiladas
  # (mantenemos "variant" en la segunda columna para distinguirlas)
  rm -f "$merged"

  if [[ "$case" == "bearing" ]]; then
    # bearing: depth_k y depth_k_rand son inestables (InvalidIntervalVectorOp),
    # así que solo comparamos base + variantes de volumen.
    ./fda_bridge_ibex "$bch" base   10 > "$merged"
    ./fda_bridge_ibex "$bch" vol_k        10 | tail -n +2 >> "$merged"
    ./fda_bridge_ibex "$bch" vol_k_rand   10 | tail -n +2 >> "$merged"
  else
    # casos normales: todas las variantes
    ./fda_bridge_ibex "$bch" base         10 > "$merged"
    ./fda_bridge_ibex "$bch" depth_k      10 | tail -n +2 >> "$merged"
    ./fda_bridge_ibex "$bch" depth_k_rand 10 | tail -n +2 >> "$merged"
    ./fda_bridge_ibex "$bch" vol_k        10 | tail -n +2 >> "$merged"
    ./fda_bridge_ibex "$bch" vol_k_rand   10 | tail -n +2 >> "$merged"
  fi
}

if [[ "$#" -eq 0 ]]; then
  # casos por defecto (solo medium/hard: los easy dan problemas en IBEX)
  run_case "medium/alkyl"
  run_case "medium/bearing"
  run_case "medium/ex6_2_14"
  run_case "hard/chem-1"
  run_case "hard/ex6_1_3"
else
  for case in "$@"; do
    run_case "$case"
  done
fi
