#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <fstream>
#include <string>
#include <cmath>

// ==================== Estructuras simples de ejemplo ====================

// Nodo genérico de buceo (aquí solo guardamos profundidad, volumen estimado y score)
struct Node {
    int depth;
    double volume_rel;   // volumen relativo respecto al inicial (V/V0)
    double heur_score;   // valor heurístico base (Feasible Diving clásico)
};

// Resultado de una ejecución
struct RunResult {
    bool reached_optimum;
    double best_value;
    int nodes_visited;
    double elapsed_seconds;
    int max_depth;
    double avg_depth;
};


// Generador global para rand(0,1)
std::mt19937_64 rng(123456);
std::uniform_real_distribution<double> unif01(0.0, 1.0);

// ==================== Simulación mínima de "evaluar nodo" ====================
// IMPORTANTE: Aquí es donde luego reemplazas por tus llamadas reales a IBEX

// Simula que el valor de la función objetivo mejora a medida que profundizamos
double evaluate_node_value(const Node& n) {
    // Ejemplo: valor empeora con profundidad (solo para tener "algo" que medir)
    return 100.0 / (1.0 + n.depth) + 0.1 * unif01(rng);
}

// Simula la generación de subnodos (ramificación)
// En tu código real, aquí iría: branch + contraction (IBEX)
std::vector<Node> generate_children(const Node& parent) {
    std::vector<Node> children;
    int num_children = 2; // por simplicidad, binario

    for (int j = 0; j < num_children; ++j) {
        Node c;
        c.depth = parent.depth + 1;
        // ejemplo de volumen que se reduce (no exactamente 1/2 para no ser tan determinista)
        double factor = 0.4 + 0.3 * unif01(rng); // entre 0.4 y 0.7
        c.volume_rel = parent.volume_rel * factor;
        c.heur_score = evaluate_node_value(c); // heurística base
        children.push_back(c);
    }
    return children;
}

// ==================== Implementación 0: Base sin temperatura ni cortes ====================

RunResult run_base(int max_nodes) {
    RunResult res{};
    res.reached_optimum = false;
    res.best_value = 1e100;
    res.nodes_visited = 0;
    res.max_depth = 0;
    res.avg_depth = 0.0;
    long long depth_accum = 0;   

    auto t0 = std::chrono::high_resolution_clock::now();

    Node current;
    current.depth = 0;
    current.volume_rel = 1.0;
    current.heur_score = evaluate_node_value(current);

    while (res.nodes_visited < max_nodes) {
        ++res.nodes_visited;

        depth_accum += current.depth;
        if (current.depth > res.max_depth) {
            res.max_depth = current.depth;
        }

        double val = evaluate_node_value(current);
        if (val < res.best_value) {
            res.best_value = val;
        }

        // Criterio de poda simple de ejemplo
        if (current.depth > 50) break;

        // Ramificación y selección por heurística clásica
        auto children = generate_children(current);
        if (children.empty()) break;

        int best_j = 0;
        double best_score = children[0].heur_score;
        for (int j = 1; j < (int)children.size(); ++j) {
            if (children[j].heur_score > best_score) {
                best_score = children[j].heur_score;
                best_j = j;
            }
        }
        current = children[best_j];
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
    // En un problema real, aquí podrías comprobar si best_value está cerca del óptimo conocido
    res.reached_optimum = true;
    
    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    } else {
        res.avg_depth = 0.0;
    }

    return res;
}

// ==================== Variante 1: profundidad + T * k ====================

RunResult run_depth_Tk(int max_nodes, int d_max, double T0, double k) {
    RunResult res{};
    res.reached_optimum = false;
    res.best_value = 1e100;
    res.nodes_visited = 0;
    res.max_depth = 0;
    res.avg_depth = 0.0;
    long long depth_accum = 0;   

    auto t0 = std::chrono::high_resolution_clock::now();

    Node current{0, 1.0, 0.0};
    current.heur_score = evaluate_node_value(current);
    double T = T0;

    while (res.nodes_visited < max_nodes) {
        ++res.nodes_visited;

        depth_accum += current.depth;
        if (current.depth > res.max_depth) {
            res.max_depth = current.depth;
        }

        double val = evaluate_node_value(current);
        if (val < res.best_value) {
            res.best_value = val;
        }

        // PODA simple
        if (current.depth >= d_max) break;

        auto children = generate_children(current);
        if (children.empty()) break;

        // Temperatura jerárquica
        for (auto& c : children) {
            double T_child = k * T / 2.0;
            c.heur_score = c.heur_score + T_child;
        }

        int best_j = 0;
        double best_score = children[0].heur_score;
        for (int j = 1; j < (int)children.size(); ++j) {
            if (children[j].heur_score > best_score) {
                best_score = children[j].heur_score;
                best_j = j;
            }
        }

        current = children[best_j];
        T = k * T / 2.0;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
    res.reached_optimum = true; // en tu código real, compara con óptimo conocido
    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    } else {
        res.avg_depth = 0.0;
    }
    return res;
}

// ==================== Variante 2: profundidad + T * k * rand(0,1) ====================

RunResult run_depth_Tk_rand(int max_nodes, int d_max, double T0, double k) {
    RunResult res{};
    res.reached_optimum = false;
    res.best_value = 1e100;
    res.nodes_visited = 0;
    res.max_depth = 0;
    res.avg_depth = 0.0;
    long long depth_accum = 0;

    auto t0 = std::chrono::high_resolution_clock::now();

    Node current{0, 1.0, 0.0};
    current.heur_score = evaluate_node_value(current);
    double T = T0;

    while (res.nodes_visited < max_nodes) {
        ++res.nodes_visited;

        depth_accum += current.depth;
        if (current.depth > res.max_depth) {
            res.max_depth = current.depth;
        }

        double val = evaluate_node_value(current);
        if (val < res.best_value) {
            res.best_value = val;
        }

        if (current.depth >= d_max) break;

        auto children = generate_children(current);
        if (children.empty()) break;

        for (auto& c : children) {
            double r = unif01(rng);
            double T_child = k * T / 2.0 * r;
            c.heur_score = c.heur_score + T_child;
        }

        int best_j = 0;
        double best_score = children[0].heur_score;
        for (int j = 1; j < (int)children.size(); ++j) {
            if (children[j].heur_score > best_score) {
                best_score = children[j].heur_score;
                best_j = j;
            }
        }

        // Temperatura para el siguiente nivel (usamos el hijo elegido)
        double r_star = unif01(rng);
        T = k * T / 2.0 * r_star;
        current = children[best_j];
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
    res.reached_optimum = true;
    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    } else {
        res.avg_depth = 0.0;
    }
    return res;
}

// ==================== Variante 3: volumen + T * k ====================

RunResult run_vol_Tk(int max_nodes, double T0, double k, double tauV) {
    RunResult res{};
    res.reached_optimum = false;
    res.best_value = 1e100;
    res.nodes_visited = 0;
    res.max_depth = 0;
    res.avg_depth = 0.0;
    long long depth_accum = 0;

    auto t0 = std::chrono::high_resolution_clock::now();

    Node current{0, 1.0, 0.0};
    current.heur_score = evaluate_node_value(current);
    double T = T0;

    while (res.nodes_visited < max_nodes) {
        ++res.nodes_visited;

        depth_accum += current.depth;
        if (current.depth > res.max_depth) {
            res.max_depth = current.depth;
        }

        double val = evaluate_node_value(current);
        if (val < res.best_value) {
            res.best_value = val;
        }

        // Corte por volumen relativo
        if (current.volume_rel <= tauV) break;

        auto children = generate_children(current);
        if (children.empty()) break;

        for (auto& c : children) {
            double T_child = k * T / 2.0;
            c.heur_score = c.heur_score + T_child;
        }

        int best_j = 0;
        double best_score = children[0].heur_score;
        for (int j = 1; j < (int)children.size(); ++j) {
            if (children[j].heur_score > best_score) {
                best_score = children[j].heur_score;
                best_j = j;
            }
        }

        current = children[best_j];
        T = k * T / 2.0;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
    res.reached_optimum = true;
    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    } else {
        res.avg_depth = 0.0;
    }
    return res;
}

// ==================== Variante 4: volumen + T * k * rand(0,1) ====================

RunResult run_vol_Tk_rand(int max_nodes, double T0, double k, double tauV) {
    RunResult res{};
    res.reached_optimum = false;
    res.best_value = 1e100;
    res.nodes_visited = 0;
    res.max_depth = 0;
    res.avg_depth = 0.0;
    long long depth_accum = 0;

    auto t0 = std::chrono::high_resolution_clock::now();

    Node current{0, 1.0, 0.0};
    current.heur_score = evaluate_node_value(current);
    double T = T0;

    while (res.nodes_visited < max_nodes) {
        ++res.nodes_visited;

        depth_accum += current.depth;
        if (current.depth > res.max_depth) {
            res.max_depth = current.depth;
        }

        double val = evaluate_node_value(current);
        if (val < res.best_value) {
            res.best_value = val;
        }

        if (current.volume_rel <= tauV) break;

        auto children = generate_children(current);
        if (children.empty()) break;

        for (auto& c : children) {
            double r = unif01(rng);
            double T_child = k * T / 2.0 * r;
            c.heur_score = c.heur_score + T_child;
        }

        int best_j = 0;
        double best_score = children[0].heur_score;
        for (int j = 1; j < (int)children.size(); ++j) {
            if (children[j].heur_score > best_score) {
                best_score = children[j].heur_score;
                best_j = j;
            }
        }

        double r_star = unif01(rng);
        T = k * T / 2.0 * r_star;
        current = children[best_j];
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    res.elapsed_seconds = std::chrono::duration<double>(t1 - t0).count();
    res.reached_optimum = true;
    if (res.nodes_visited > 0) {
        res.avg_depth = static_cast<double>(depth_accum) / res.nodes_visited;
    } else {
        res.avg_depth = 0.0;
    }
    return res;
}

// ============================================================================
// MAIN COMPLETO PARA EJECUCIÓN DE EXPERIMENTOS
// ============================================================================

int main(int argc, char** argv) {

    // ------------------------------------------------------------------------
    // 1. Validación de argumentos
    // ------------------------------------------------------------------------
    if (argc < 6) {
        std::cerr << "Uso: ./fda_experimentos <variante> <reps> <max_nodes> <T0> <k> [tauV]\n";
        std::cerr << "Ejemplos:\n";
        std::cerr << "  ./fda_experimentos base       10 50000 0 0\n";
        std::cerr << "  ./fda_experimentos depth_k    10 50000 1.0 0.7\n";
        std::cerr << "  ./fda_experimentos depth_k_rand 10 50000 1.0 0.7\n";
        std::cerr << "  ./fda_experimentos vol_k      10 50000 1.0 0.7 0.10\n";
        std::cerr << "  ./fda_experimentos vol_k_rand 10 50000 1.0 0.7 0.10\n";
        return 1;
    }

    std::string variant = argv[1];
    int reps            = std::stoi(argv[2]);
    int max_nodes       = std::stoi(argv[3]);
    double T0           = std::stod(argv[4]);
    double k_param      = std::stod(argv[5]);
    double tauV         = (argc >= 7 ? std::stod(argv[6]) : 0.10);

    int d_max = 30;  // Puedes ajustarlo más adelante según tus pruebas

    // ------------------------------------------------------------------------
    // 2. Encabezado para CSV o salida estructurada
    // ------------------------------------------------------------------------
    std::cout << "run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal\n";

    // ------------------------------------------------------------------------
    // 3. Bucle de repeticiones
    // ------------------------------------------------------------------------
    for (int r = 0; r < reps; r++) {

        RunResult res;

        if (variant == "base") {
            res = run_base(max_nodes);

        } else if (variant == "depth_k") {
            res = run_depth_Tk(max_nodes, d_max, T0, k_param);

        } else if (variant == "depth_k_rand") {
            res = run_depth_Tk_rand(max_nodes, d_max, T0, k_param);

        } else if (variant == "vol_k") {
            res = run_vol_Tk(max_nodes, T0, k_param, tauV);

        } else if (variant == "vol_k_rand") {
            res = run_vol_Tk_rand(max_nodes, T0, k_param, tauV);

        } else {
            std::cerr << "ERROR: Variante no reconocida: " << variant << "\n";
            return 1;
        }

        // --------------------------------------------------------------------
        // 4. Imprimir resultado en formato CSV
        // --------------------------------------------------------------------
        std::cout   << r << ","
                    << variant << ","
                    << res.best_value << ","
                    << res.nodes_visited << ","
                    << res.elapsed_seconds << ","
                    << res.max_depth << ","
                    << res.avg_depth << ","
                    << (res.reached_optimum ? 1 : 0)
                    << "\n";
    }

    return 0;
}
