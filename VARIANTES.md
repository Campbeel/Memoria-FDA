# Variantes FDA + IBEX (explicación corta)

Todas las variantes comparten utilidades en `fda_bridge_common.{h,cpp}`:
- `OptContext`: sistema IBEX original + contratista `CtcHC4`.
- `contract_with_goal`: contrae la caja actual y evalúa el objetivo (guarda intervalo `goal_bounds`).
- `bisect_box`: parte la caja por la variable más ancha.
- `heur_diving_score`: heurística = `-lb(goal)` si hay objetivo.
- `adaptive_k`: ajusta k si mejora la cota superior global.

La estructura del *dive* es similar: contracción + elección de hijo por heurística + actualización de temperatura/volumen/profundidad. El campo `optimal` se marca si hay una UB finita al final del dive/BB.

## Base (`FD_base.cpp`)
- Dive simple: contracción HC4, actualiza UB con el objetivo, bisecta y elige el hijo con mejor score.
- Parada: caja < `eps_box` (solo si ya bajaste un nivel) o `max_iters`.
- `run_fda_base_bb`: B&B que encola cajas, llama al dive en cada nodo y actualiza UB global.

## Depth_k (`FDA_depth_k.cpp`)
- Igual que base, pero con temperatura determinista: `T_child = k*T/2`.
- Límite de profundidad `d_max`; si se alcanza, extiende el presupuesto para seguir explorando más hondo.
- Variante `*_bb`: B&B que propaga T en los nodos.

## Depth_k_rand (`FDA_depth_k_rand.cpp`)
- Como `depth_k`, pero añade ruido a `T_child` y perturba `d_max` por ejecución.
- Variante `*_bb`: mismo esquema B&B con T heredada.

## Vol_k (`FDA_vol_k.cpp`)
- Usa volumen relativo: `Vrel = volume(box)/V0` (V0 se actualiza al nodo elegido).
- Cortes por volumen solo después de `depth_vol>=5`; umbral `tauV = eps_V * exp(-beta*depth)`.
- Temperatura determinista `T_child = k*T/2`.
- Variante `*_bb`: B&B con T heredada.

## Vol_k_rand (`FDA_vol_k_rand.cpp`)
- Igual que `vol_k`, pero `tauV` se perturba con un factor aleatorio y `T_child` incluye ruido.
- Variante `*_bb`: B&B con T heredada.

## Convenciones de salida
CSV: `run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal`
- `best_value`: mejor UB hallada.
- `optimal`: 1 si `best_value` es finita al terminar la corrida.
