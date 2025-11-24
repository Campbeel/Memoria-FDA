cmake_minimum_required(VERSION 3.10)
project(memoria_fda CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Habilita AddressSanitizer para cazar corrupción de memoria en fda_bridge_ibex.
add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
link_libraries(-fsanitize=address)

# === Ruta base de IBEX (ajusta si la tienes en otra carpeta) ===
set(IBEX_ROOT "/home/benjamin-mu-oz/ibex-lib")

# === Recolectar automáticamente TODOS los subdirectorios de IBEX ===
file(GLOB IBEX_ALL_DIRS
    "${IBEX_ROOT}/src"
    "${IBEX_ROOT}/src/*"
    "${IBEX_ROOT}/build/src"
    "${IBEX_ROOT}/build/src/*"
    "${IBEX_ROOT}/3rd"
    "${IBEX_ROOT}/3rd/*"
    "${IBEX_ROOT}/plugins"
    "${IBEX_ROOT}/plugins/*"

    # 💡 Añadir la carpeta del wrapper LP de IBEX
    "${IBEX_ROOT}/lp_lib_wrapper"
    "${IBEX_ROOT}/lp_lib_wrapper/*"
)

set(IBEX_INCLUDE_DIRS "")
foreach(d ${IBEX_ALL_DIRS})
    if (IS_DIRECTORY "${d}")
        list(APPEND IBEX_INCLUDE_DIRS "${d}")
    endif()
endforeach()

# === Includes: IBEX + tus wrappers + CLP + FILIB ===
include_directories(
    ${IBEX_INCLUDE_DIRS}

    # Wrappers FDA (relativos a memoria_cofigo/)
    ../interval_lib_wrapper
    ../interval_lib_wrapper/direct
    ../interval_lib_wrapper/bias
    ../interval_lib_wrapper/gaol
    ../interval_lib_wrapper/filib

    # CLP (ClpSimplex.hpp)
    /usr/include/coin

    # FILIB (interval/interval.hpp) si lo tienes aquí
    /home/benjamin-mu-oz/filib-src
)

# Librerías compiladas de IBEX (+ gaol / Mathlib instaladas por IBEX)
link_directories(
    ${IBEX_ROOT}/build/src
    /usr/local/lib
    /usr/local/lib/ibex/3rd
)

# === Ejecutable de experimentos sintéticos ===
add_executable(fda_experimentos
    fda_experimentos.cpp
)
target_link_libraries(fda_experimentos ibex)

# === Ejecutable puente IBEX + FDA ===
add_executable(fda_bridge_ibex
    fda_bridge_common.cpp
    FD_base.cpp
    FDA_depth_k.cpp
    FDA_depth_k_rand.cpp
    FDA_vol_k.cpp
    FDA_vol_k_rand.cpp
    fda_bridge_ibex.cpp
)
# enlazar también contra CLP y CoinUtils (solver LP que usa IBEX)
target_link_libraries(fda_bridge_ibex ibex gaol ultim Clp CoinUtils)
