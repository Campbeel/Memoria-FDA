cmake_minimum_required(VERSION 3.10)
project(memoria_demo CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
file(GLOB IBEX_SRC_SUBDIRS "${IBEX_ROOT}/src/*" "${IBEX_ROOT}/build/src/*")
include_directories(
    ${IBEX_ROOT}/src
    ${IBEX_ROOT}/build/src
    ${IBEX_SRC_SUBDIRS}
    ${IBEX_ROOT}/3rd
    ${IBEX_ROOT}/plugins
    ${IBEX_ROOT}/lp_lib_wrapper
    ${IBEX_ROOT}/lp_lib_wrapper/include
    ${IBEX_ROOT}/lp_lib_wrapper/clp
    ${IBEX_ROOT}/lp_lib_wrapper/soplex
    ${IBEX_ROOT}/lp_lib_wrapper/none
    ${IBEX_ROOT}/lp_lib_wrapper/cplex
    ${IBEX_ROOT}/interval_lib_wrapper
    ${IBEX_ROOT}/interval_lib_wrapper/direct
    ${IBEX_ROOT}/interval_lib_wrapper/bias
    ${IBEX_ROOT}/interval_lib_wrapper/gaol
    ${IBEX_ROOT}/interval_lib_wrapper/filib

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

# Ejecutable ibexsolve mínimo
add_executable(mini_ibexsolve
    ../mini_ibexsolve.cpp
)
target_link_libraries(mini_ibexsolve ibex gaol ultim Clp CoinUtils)

# Ejecutable FDA (reutiliza las fuentes del proyecto principal)
add_executable(fda_bridge_ibex_demo
    ../fda_bridge_common.cpp
    ../FD_base.cpp
    ../FDA_depth_k.cpp
    ../FDA_depth_k_rand.cpp
    ../FDA_vol_k.cpp
    ../FDA_vol_k_rand.cpp
    ../fda_bridge_ibex.cpp
)
target_link_libraries(fda_bridge_ibex_demo ibex gaol ultim Clp CoinUtils)
