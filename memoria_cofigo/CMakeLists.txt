cmake_minimum_required(VERSION 3.10)
project(memoria_fda CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# === Ruta base de IBEX (ajusta si la tienes en otra carpeta) ===
set(IBEX_ROOT "/home/benjamin-mu-oz/ibex-lib")
# Build de Ibex con SoPlex (para evitar crash de CLP)
set(IBEX_BUILD_DIR "${IBEX_ROOT}/build-soplex")
# Usa la lib SoPlex instalada en el sistema (evita depender del build interno).
set(SOPLEX_LIB "/usr/lib/x86_64-linux-gnu/libsoplex.a")

# === Recolectar automáticamente TODOS los subdirectorios de IBEX ===
file(GLOB IBEX_ALL_DIRS
    "${IBEX_ROOT}/src"
    "${IBEX_ROOT}/src/*"
    "${IBEX_ROOT}/build/src"
    "${IBEX_ROOT}/build/src/*"
    "${IBEX_BUILD_DIR}/src"
    "${IBEX_BUILD_DIR}/src/*"
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
    ${IBEX_BUILD_DIR}/src
    /usr/local/lib
    /usr/local/lib/ibex/3rd
)

set(MEMORIA_LIBS ibex gaol ultim ${SOPLEX_LIB} soplex gmp gmpxx z)

add_executable(ibex_opt_base ibex_opt_base.cpp)
target_link_libraries(ibex_opt_base ${MEMORIA_LIBS} Clp CoinUtils)

# Ejecutable "ibexopt" original (clonado de src/bin/ibexopt.cpp)
add_executable(ibex_opt_full ibex_opt_full.cpp)
target_link_libraries(ibex_opt_full ${MEMORIA_LIBS} Clp CoinUtils)

add_executable(ibex_menu run_menu.cpp)
target_link_libraries(ibex_menu pthread)
