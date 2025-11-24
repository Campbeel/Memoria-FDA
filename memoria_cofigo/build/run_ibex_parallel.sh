#!/usr/bin/env bash
set -e

# Problemas (ruta relativa dentro de casos/, sin .bch)
PROBLEMS=(medium/alkyl medium/bearing medium/ex6_2_14 hard/chem-1 hard/ex6_1_3)

# Variantes de Feasible Diving
VARIANTS=(base depth_k depth_k_rand vol_k vol_k_rand)

RUNS=10        # número de repeticiones por (problema, variante)
MAX_JOBS=5     # cuántos procesos simultáneos quieres

mkdir -p resultados

for p in "${PROBLEMS[@]}"; do
  for v in "${VARIANTS[@]}"; do
    safe_name="${p//\//_}"
    out="resultados/${safe_name}_${v}.csv"
    echo "Lanzando problema=${p}, variante=${v} -> ${out}"

    ./fda_bridge_ibex "../casos/${p}.bch" "${v}" "${RUNS}" > "${out}" &

    # limitar número de jobs en paralelo
    while (( $(jobs -r | wc -l) >= MAX_JOBS )); do
      sleep 0.5
    done
  done
done

wait
echo ">>> Todas las corridas han finalizado. Resultados en 'resultados/'."
