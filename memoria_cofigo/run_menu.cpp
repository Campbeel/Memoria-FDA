// run_menu.cpp
// Pequeño menú para lanzar variantes FDA o ibex_opt_full en paralelo.

#include <atomic>
#include <algorithm>
#include <chrono>
#include <limits>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <mutex>
#include <string>
#include <sstream>
#include <thread>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

struct VariantInfo {
    string name;
    string binary;
    string fd_mode; // opcional: modo FD para ibex_opt_full
};

static const VariantInfo kVariants[] = {
    {"base",        "./ibex_opt_full",      ""            },
    {"vol_k",       "./ibex_opt_full",      "vol_k"       },
    {"vol_k_rand",  "./ibex_opt_full",      "vol_k_rand"  },
    {"depth_k",     "./ibex_opt_full",      "depth_k"     },
    {"depth_k_rand","./ibex_opt_full",      "depth_k_rand"}
};

// Ejecuta ibex_opt_full y extrae nodos, tiempo, loup; formato CSV.
string run_ibex_base(const string& cmd, const string& variant, int run_id) {
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "";
    char buf[4096];
    string line;
    double best = 0.0;
    long nodes = -1;
    double elapsed = 0.0;
    bool optimal = false;
    while (fgets(buf, sizeof(buf), pipe)) {
        line = buf;
        if (line.find("nodes (cells):") != string::npos) {
            sscanf(line.c_str(), "nodes (cells):%ld", &nodes);
        } else if (line.find("cpu time used:") != string::npos) {
            sscanf(line.c_str(), " cpu time used:%lf", &elapsed);
        } else if (line.find("f* in") != string::npos) {
            double lb, ub;
            if (sscanf(line.c_str(), " f* in\t[%lf,%lf]", &lb, &ub) == 2) {
                best = ub;
            }
        } else if (line.find("optimization successful") != string::npos) {
            optimal = true;
        }
    }
    pclose(pipe);
    // CSV: run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal
    stringstream ss;
    ss << run_id << "," << variant << "," << best << "," << nodes << "," << elapsed
       << ",NA,NA," << (optimal ? 1 : 0);
    return ss.str();
}

int main() {
    cout << "Selecciona variante:\n";
    for (size_t i = 0; i < sizeof(kVariants)/sizeof(kVariants[0]); ++i) {
        cout << " " << (i+1) << ". " << kVariants[i].name << "\n";
    }
    cout << " 6. Correr todo (5 variantes x 20 problemas x 10 repeticiones)\n";

    int choice = 0;
    cout << "Opcion: ";
    cin >> choice;
    if (choice < 1 || choice > 6) {
        cerr << "Opcion invalida\n";
        return 1;
    }

    auto list_bchs = [](const string& folder) {
        vector<fs::path> files;
        fs::path base = fs::path("..") / "casos" / folder;
        if (fs::exists(base) && fs::is_directory(base)) {
            for (auto& p : fs::directory_iterator(base)) {
                if (p.path().extension() == ".bch") files.push_back(p.path());
            }
        }
        sort(files.begin(), files.end());
        vector<string> out;
        out.reserve(files.size());
        for (auto& f : files) out.push_back(f.string());
        return out;
    };

    if (choice == 6) {
        cout << "Elegir variantes a recorrer (ej: 1,3,5). Enter = todas:\n";
        for (size_t i = 0; i < sizeof(kVariants)/sizeof(kVariants[0]); ++i) {
            cout << " " << (i+1) << ") " << kVariants[i].name << "\n";
        }
        cout << "Seleccion: ";
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        string sel_line;
        getline(cin, sel_line);

        vector<const VariantInfo*> selected;
        auto add_variant = [&](size_t idx) {
            if (idx < sizeof(kVariants)/sizeof(kVariants[0])) {
                selected.push_back(&kVariants[idx]);
            }
        };
        if (sel_line.empty()) {
            for (size_t i = 0; i < sizeof(kVariants)/sizeof(kVariants[0]); ++i) add_variant(i);
        } else {
            string token;
            stringstream ss(sel_line);
            while (getline(ss, token, ',')) {
                if (token.empty()) continue;
                size_t pos = token.find_first_not_of(" \t");
                if (pos != string::npos) token = token.substr(pos);
                pos = token.find_last_not_of(" \t");
                if (pos != string::npos) token = token.substr(0, pos+1);
                int idx = atoi(token.c_str());
                if (idx >= 1 && idx <= (int)(sizeof(kVariants)/sizeof(kVariants[0]))) {
                    add_variant((size_t)(idx-1));
                }
            }
            if (selected.empty()) {
                cerr << "Sin selección válida; usando todas las variantes.\n";
                for (size_t i = 0; i < sizeof(kVariants)/sizeof(kVariants[0]); ++i) add_variant(i);
            }
        }

        vector<string> problems = list_bchs("medium");
        auto hard = list_bchs("hard");
        problems.insert(problems.end(), hard.begin(), hard.end());
        if (problems.empty()) {
            cerr << "No se encontraron problemas en casos/medium ni casos/hard\n";
            return 1;
        }

        fs::create_directories("results");
        struct Task {
            const VariantInfo* vinfo;
            string problem;
            int run_id;
        };
        vector<Task> tasks;
        tasks.reserve(problems.size() * selected.size() * 10);
        for (const VariantInfo* vp : selected) {
            for (const string& prob : problems) {
                for (int r = 0; r < 10; ++r) {
                    tasks.push_back({vp, prob, r});
                }
            }
        }

        const unsigned max_parallel = 10;
        atomic<size_t> next{0};
        mutex m;

        auto worker_all = [&]() {
            while (true) {
                size_t idx = next.fetch_add(1);
                if (idx >= tasks.size()) break;
                const Task& t = tasks[idx];
                string cmd = t.vinfo->binary + " " + t.problem;
                if (!t.vinfo->fd_mode.empty()) {
                    cmd += " --fd-mode=" + t.vinfo->fd_mode;
                }
                string line = run_ibex_base(cmd, t.vinfo->name, t.run_id);

                string prob_name = fs::path(t.problem).filename().string();
                string csv = "results/results_" + t.vinfo->name + "_all.csv";
                lock_guard<mutex> lk(m);
                bool exists = fs::exists(csv);
                ofstream o(csv, ios::app);
                if (!exists) {
                    o << "run,variant,problem,best_value,nodes,elapsed,max_depth,avg_depth,optimal\n";
                }
                if (!line.empty()) {
                    size_t first = line.find(',');
                    size_t second = line.find(',', first+1);
                    if (first!=string::npos && second!=string::npos) {
                        line.insert(second, "," + prob_name);
                    }
                    o << line;
                    if (line.back() != '\n') o << "\n";
                } else {
                    o << t.run_id << "," << t.vinfo->name << "," << prob_name
                      << ",nan,-1,-1,NA,NA,0\n";
                }
            }
        };

        vector<thread> pool;
        pool.reserve(max_parallel);
        auto start_all = chrono::steady_clock::now();
        for (unsigned i = 0; i < max_parallel; ++i) {
            pool.emplace_back(worker_all);
        }
        for (thread& th : pool) th.join();
        auto end_all = chrono::steady_clock::now();
        double elapsed_all = chrono::duration<double>(end_all - start_all).count();
        cout << "Ejecuciones completas (" << tasks.size() << " corridas) en "
             << elapsed_all << " s. Revisa results/results_<variant>_all.csv\n";
        return 0;
    }

    const VariantInfo& vinfo = kVariants[choice-1];
    string problem;

    cout << "Selecciona origen del .bch:\n";
    cout << " 1) casos/medium\n";
    cout << " 2) casos/hard\n";
    cout << " 3) ruta manual\n";
    int src_choice = 0;
    cout << "Opcion: ";
    cin >> src_choice;

    auto pick_file_from = [&](const string& folder) -> string {
        vector<string> files = list_bchs(folder);
        if (files.empty()) {
            cerr << "No se encontraron .bch en " << (fs::path("..") / "casos" / folder) << "\n";
            return "";
        }
        cout << "Archivos en " << (fs::path("..") / "casos" / folder) << ":\n";
        for (size_t i=0;i<files.size();++i) {
            cout << " " << (i+1) << ") " << fs::path(files[i]).filename().string() << "\n";
        }
        int opt=0;
        cout << "Elige: ";
        cin >> opt;
        if (opt<1 || (size_t)opt>files.size()) {
            cerr << "Opcion invalida\n";
            return "";
        }
        return files[opt-1];
    };

    if (src_choice==1) {
        problem = pick_file_from("medium");
    } else if (src_choice==2) {
        problem = pick_file_from("hard");
    } else {
        cout << "Ruta al .bch: ";
        cin >> problem;
    }
    if (problem.empty()) return 1;

    unsigned hw = thread::hardware_concurrency();
    unsigned max_parallel = hw > 2 ? hw - 2 : 1;
    cout << "Max paralelo sugerido: " << max_parallel << "\n";
    int runs = 1;
    cout << "Cuantas corridas? (max " << max_parallel << "): ";
    cin >> runs;
    if (runs < 1) runs = 1;
    if ((unsigned) runs > max_parallel) runs = max_parallel;
    unsigned par = runs;

    fs::create_directories("results");
    string base_name = fs::path(problem).stem().string();
    string csv = "results/results_" + vinfo.name + "_" + base_name + ".csv";
    ofstream ofs(csv, ios::out);
    ofs << "run,variant,problem,best_value,nodes,elapsed,max_depth,avg_depth,optimal\n";
    ofs.close();

    mutex m;
    vector<thread> pool;
    int launched = 0;

    auto worker = [&](int idx) {
        string cmd = vinfo.binary + " " + problem;
        if (!vinfo.fd_mode.empty()) {
            cmd += " --fd-mode=" + vinfo.fd_mode;
        }
        string line = run_ibex_base(cmd, vinfo.name, idx);
        lock_guard<mutex> lk(m);
        ofstream o(csv, ios::app);
        if (!line.empty()) {
            string prob_name = fs::path(problem).filename().string();
            size_t first = line.find(',');
            size_t second = line.find(',', first+1);
            if (first!=string::npos && second!=string::npos) {
                line.insert(second, "," + prob_name);
            }
            o << line;
            if (line.back() != '\n') o << "\n";
        } else {
            o << idx << "," << vinfo.name << "," << fs::path(problem).filename().string()
              << ",NA,NA,NA,NA,NA,0\n";
        }
    };

    auto start = chrono::steady_clock::now();
    while (launched < runs) {
        while (pool.size() < par && launched < runs) {
            int idx = launched;
            pool.emplace_back(worker, idx);
            launched++;
        }
        for (auto it = pool.begin(); it != pool.end();) {
            it->join();
            it = pool.erase(it);
        }
    }
    auto end = chrono::steady_clock::now();
    double elapsed = chrono::duration<double>(end - start).count();
    cout << "Hecho. CSV en " << csv << ". Tiempo total: " << elapsed << " s\n";
    return 0;
}
