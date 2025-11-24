#!/usr/bin/env bash
set -euo pipefail

# Directorio donde estamos
BUILD_DIR="$(pwd)"
CASES_DIR="${BUILD_DIR}/../casos"
OUT_DIR="${BUILD_DIR}/resultados"

mkdir -p "${OUT_DIR}"

# Problemas y variantes
PROBLEMS=(ackley_5 alkyl bearing chem haverly)
VARIANTS=(base depth_k depth_k_rand vol_k vol_k_rand)

NUM_RUNS=10         # 10 iteraciones por problema/variante
MAX_PARALLEL=10     # máximo 10 procesos en paralelo

job_count=0

for prob in "${PROBLEMS[@]}"; do
  for var in "${VARIANTS[@]}"; do
    in_bch="${CASES_DIR}/${prob}.bch"
    out_csv="${OUT_DIR}/${prob}_${var}.csv"

    echo ">>> Lanzando ${prob} - ${var} (${NUM_RUNS} runs) -> ${out_csv}"

    # Ejecutar en background
    ./fda_bridge_ibex "${in_bch}" "${var}" "${NUM_RUNS}" > "${out_csv}" &

    ((job_count++))

    # Cada vez que llegamos a MAX_PARALLEL, esperamos a que terminen
    if (( job_count % MAX_PARALLEL == 0 )); then
      echo ">>> Alcanzado límite de ${MAX_PARALLEL} jobs. Esperando..."
      wait
      echo ">>> Continuando con el siguiente bloque de trabajos."
    fi
  done
done

# Esperar a que terminen los últimos jobs
wait
echo ">>> TODAS las pruebas han terminado. Resultados en: ${OUT_DIR}"
