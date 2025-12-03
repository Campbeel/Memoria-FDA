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
#include <random>

using namespace std;
namespace fs = std::filesystem;

struct VariantInfo {
    string name;
    string binary;
    string fd_mode; // opcional: modo FD para ibex_opt_full
    double timeout = 0.0; // timeout opcional en segundos
};

static const VariantInfo kVariants[] = {
// Timeouts en segundos; 120s = 2 min.
{"base",        "./ibex_opt_full",      "",           120.0 },
{"vol_k",       "./ibex_opt_full",      "vol_k",      120.0 },
{"vol_k_rand",  "./ibex_opt_full",      "vol_k_rand", 120.0 },
{"depth_k",     "./ibex_opt_full",      "depth_k",    120.0 },
{"depth_k_rand","./ibex_opt_full",      "depth_k_rand",120.0}
};

struct RunResult {
    string csv_line;
    long triggers_total = -1;
    long triggers_depth = -1;
    long triggers_vol = -1;
    long vol_eval = -1;
    long vol_nonfinite = -1;
};

// Ejecuta ibex_opt_full y extrae nodos, tiempo, loup; formato CSV.
RunResult run_ibex_base(const string& cmd, const string& variant, int run_id) {
    FILE* pipe = popen(cmd.c_str(), "r");
    RunResult res;
    if (!pipe) return res;
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
        } else if (line.find("fd_triggers_total:") != string::npos) {
            try { res.triggers_total = stol(line.substr(line.find(':')+1)); } catch(...) {}
        } else if (line.find("fd_triggers_depth:") != string::npos) {
            try { res.triggers_depth = stol(line.substr(line.find(':')+1)); } catch(...) {}
        } else if (line.find("fd_triggers_vol:") != string::npos) {
            try { res.triggers_vol = stol(line.substr(line.find(':')+1)); } catch(...) {}
        } else if (line.find("fd_vol_eval:") != string::npos) {
            try { res.vol_eval = stol(line.substr(line.find(':')+1)); } catch(...) {}
        } else if (line.find("fd_vol_nonfinite:") != string::npos) {
            try { res.vol_nonfinite = stol(line.substr(line.find(':')+1)); } catch(...) {}
        } else if (line.find("fd_triggers:") != string::npos || line.find("triggers_csv:") != string::npos) {
            // Compatibilidad con versiones anteriores: tratamos como total.
            size_t pos = line.find("fd_triggers:");
            int offset = 12;
            if (pos == string::npos) {
                pos = line.find("triggers_csv:");
                offset = 13;
            }
            if (pos != string::npos) {
                try {
                    long t = stol(line.substr(pos + offset));
                    res.triggers_total = t;
                } catch(...) {}
            }
        }
    }
    pclose(pipe);
        // CSV: run,variant,best_value,nodes,elapsed,max_depth,avg_depth,optimal,triggers
    stringstream ss;
    ss << run_id << "," << variant << "," << best << "," << nodes << "," << elapsed
       << ",NA,NA," << (optimal ? 1 : 0) << ",";
    if (res.triggers_total>=0) ss << res.triggers_total; else ss << "0";
    ss << ",";
    if (res.triggers_depth>=0) ss << res.triggers_depth; else ss << "0";
    ss << ",";
    if (res.triggers_vol>=0) ss << res.triggers_vol; else ss << "0";
    ss << ",";
    if (res.vol_eval>=0) ss << res.vol_eval; else ss << "0";
    ss << ",";
    if (res.vol_nonfinite>=0) ss << res.vol_nonfinite; else ss << "0";
    res.csv_line = ss.str();
    return res;
}

int main() {
    auto now_str = []() {
        auto now = chrono::system_clock::now();
        time_t tt = chrono::system_clock::to_time_t(now);
        char buf[64];
        strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", localtime(&tt));
        return string(buf);
    };

    cout << "Selecciona variante:\n";
    for (size_t i = 0; i < sizeof(kVariants)/sizeof(kVariants[0]); ++i) {
        cout << " " << (i+1) << ". " << kVariants[i].name << "\n";
    }
    cout << " 6. Correr todo (5 variantes x 30 problemas x 10 repeticiones)\n";

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

    auto prompt_seed_base = []() -> long long {
        cout << "Semilla base (Enter o 0 = aleatoria distinta por corrida): ";
        string line;
        getline(cin >> ws, line);
        if (line.empty()) return -1;
        long long s = atoll(line.c_str());
        return s > 0 ? s : -1;
    };

    auto make_seed_for_run = [&](long long seed_base) {
        return [seed_base](int run_id) -> long long {
            if (seed_base > 0) return seed_base + run_id;
            static std::random_device rd;
            static std::mt19937_64 gen(rd());
            uint64_t salt = static_cast<uint64_t>(run_id) * 0x9e3779b97f4a7c15ULL;
            return static_cast<long long>(gen() ^ salt);
        };
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

        vector<string> medium = list_bchs("medium");
        vector<string> hard = list_bchs("hard");
        // Si existe ex17_2_7.bch en cualquier carpeta, agrégalo.
        for (const auto& p : {fs::path("..")/"casos"/"medium"/"ex17_2_7.bch",
                              fs::path("..")/"casos"/"hard"/"ex17_2_7.bch"}) {
            if (fs::exists(p)) medium.push_back(p.string());
        }
        // Excluir casos muy lentos si molestan en el barrido.
        const vector<string> exclude = {"ex14_2_7.bch"};
        auto apply_exclude = [&](vector<string>& v) {
            v.erase(remove_if(v.begin(), v.end(), [&](const string& p) {
                string fname = fs::path(p).filename().string();
                return find(exclude.begin(), exclude.end(), fname) != exclude.end();
            }), v.end());
        };
        apply_exclude(medium);
        apply_exclude(hard);
        sort(medium.begin(), medium.end());
        medium.erase(unique(medium.begin(), medium.end()), medium.end());
        sort(hard.begin(), hard.end());
        hard.erase(unique(hard.begin(), hard.end()), hard.end());
        if (medium.empty() && hard.empty()) {
            cerr << "No se encontraron problemas en casos/medium ni casos/hard\n";
            return 1;
        }

        cout << "¿Qué subconjunto quieres correr?\n";
        cout << " 1) Solo medium\n";
        cout << " 2) Solo hard\n";
        cout << " 3) Medium y luego hard\n";
        int subset = 3;
        cout << "Opcion: ";
        cin >> subset;
        bool run_medium = (subset == 1 || subset == 3);
        bool run_hard   = (subset == 2 || subset == 3);
        if (run_medium && medium.empty()) run_medium = false;
        if (run_hard && hard.empty())     run_hard = false;
        if (!run_medium && !run_hard) {
            cerr << "Nada que correr (verifica la selección).\n";
            return 1;
        }

        long long seed_base = prompt_seed_base();
        auto seed_for_run = make_seed_for_run(seed_base);

        fs::create_directories("results");
    auto run_problem = [&](const string& prob_path) {
        string prob_name = fs::path(prob_path).filename().string();
        // Excluir problemas lentos conocidos
        static const std::vector<std::string> kSkip = {"schwefel5.bch", "schwefel5-abs.bch", "ex8_5_2_1.bch", "ex7_3_4.bch"};
        if (std::find(kSkip.begin(), kSkip.end(), prob_name) != kSkip.end()) {
            cout << "  [skip] " << prob_name << "\n";
            return;
        }
            for (const VariantInfo* vp : selected) {
                // CSV por problema y variante (se reinicia en cada problema).
                string csv = "results/results_" + vp->name + "_" + fs::path(prob_name).stem().string() + ".csv";
                string trig_csv = "results/triggers_" + vp->name + "_" + fs::path(prob_name).stem().string() + ".csv";
                if (fs::exists(csv)) fs::remove(csv);
                if (fs::exists(trig_csv)) fs::remove(trig_csv);
                ofstream o(csv, ios::out);
                o << "run,variant,problem,best_value,nodes,elapsed,max_depth,avg_depth,optimal,triggers,triggers_depth,triggers_vol,vol_eval,vol_nonfinite\n";
                o.close();
                ofstream ot(trig_csv, ios::out);
                ot << "run,variant,problem,triggers_total,triggers_depth,triggers_vol,vol_eval,vol_nonfinite\n";
                ot.close();

                cout << " Ejecutando " << vp->name << " sobre " << prob_name << "...\n";
                const int runs = 10;
                unsigned hw = thread::hardware_concurrency();
                unsigned max_parallel = hw ? std::max(1u, hw/2) : 4u;
                // respetar un paralelo máximo de cores físicos si se detecta SMT (aprox hw/2)
                if (hw > 0) {
                    unsigned phys_approx = std::max(1u, hw/2);
                    if (max_parallel > phys_approx) max_parallel = phys_approx;
                }
                if (max_parallel > 10) max_parallel = 10;
                // Si es un problema hard, limitar aún más el paralelismo para evitar contención.
                bool is_hard = (prob_path.find("/hard/") != string::npos);
                if (is_hard && max_parallel > 2) max_parallel = 2;
                struct Task { int r; };
                vector<Task> tasks;
                tasks.reserve(runs);
                for (int r=0; r<runs; ++r) tasks.push_back({r});
                atomic<size_t> next{0};
                mutex m;
                auto worker = [&]() {
                    while (true) {
                        size_t idx = next.fetch_add(1);
                        if (idx >= tasks.size()) break;
                        int r = tasks[idx].r;
                        string cmd = vp->binary + " " + prob_path;
                        long long seed = seed_for_run(r); // semilla controlada por usuario
                        cmd += " --random-seed=" + std::to_string(seed);
                        if (!vp->fd_mode.empty()) cmd += " --fd-mode=" + vp->fd_mode;
                        if (vp->timeout > 0.0)   cmd += " --timeout=" + std::to_string(vp->timeout);
                        RunResult res = run_ibex_base(cmd, vp->name, r);
                        lock_guard<mutex> lk(m);
                        ofstream out(csv, ios::app);
                        ofstream outt(trig_csv, ios::app);
                        if (!res.csv_line.empty()) {
                            string line = res.csv_line;
                            size_t first = line.find(',');
                            size_t second = line.find(',', first+1);
                            if (first!=string::npos && second!=string::npos) {
                                line.insert(second, "," + prob_name);
                            }
                            out << line;
                            if (line.back() != '\n') out << "\n";
                            outt << r << "," << vp->name << "," << prob_name << ","
                                 << (res.triggers_total>=0?res.triggers_total:0) << ","
                                 << (res.triggers_depth>=0?res.triggers_depth:0) << ","
                                 << (res.triggers_vol>=0?res.triggers_vol:0) << ","
                                 << (res.vol_eval>=0?res.vol_eval:0) << ","
                                 << (res.vol_nonfinite>=0?res.vol_nonfinite:0) << "\n";
                        } else {
                            out << r << "," << vp->name << "," << prob_name
                                << ",nan,-1,-1,NA,NA,0,0,0,0,0,0\n";
                            outt << r << "," << vp->name << "," << prob_name << ",0,0,0,0,0\n";
                        }
                    }
                };
                vector<thread> pool;
                for (unsigned i=0;i<max_parallel;++i) pool.emplace_back(worker);
                for (auto& th: pool) th.join();
            }
        };

        cout << "Inicio batch: " << now_str() << "\n";
        // Primero medium (si aplica), luego hard, secuencialmente por problema.
        if (run_medium) {
            for (const string& p : medium) run_problem(p);
        }
        if (run_hard) {
            for (const string& p : hard)   run_problem(p);
        }
        cout << "Fin batch: " << now_str() << "\n";
        cout << "Listo. Revisa results/results_<variante>_<problema>.csv\n";
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

    long long seed_base = prompt_seed_base();
    auto seed_for_run = make_seed_for_run(seed_base);

    fs::create_directories("results");
    string base_name = fs::path(problem).stem().string();
    string csv = "results/results_" + vinfo.name + "_" + base_name + ".csv";
    string trig_csv = "results/triggers_" + vinfo.name + "_" + base_name + ".csv";
    if (fs::exists(csv)) fs::remove(csv);
    if (fs::exists(trig_csv)) fs::remove(trig_csv);
    ofstream ofs(csv, ios::out);
    ofs << "run,variant,problem,best_value,nodes,elapsed,max_depth,avg_depth,optimal,triggers,triggers_depth,triggers_vol,vol_eval,vol_nonfinite\n";
    ofs.close();
    ofstream ot(trig_csv, ios::out);
    ot << "run,variant,problem,triggers_total,triggers_depth,triggers_vol,vol_eval,vol_nonfinite\n";
    ot.close();

    mutex m;
    vector<thread> pool;
    int launched = 0;

    auto worker = [&](int idx) {
        string cmd = vinfo.binary + " " + problem;
        long long seed = seed_for_run(idx);
        cmd += " --random-seed=" + std::to_string(seed);
        if (!vinfo.fd_mode.empty()) {
            cmd += " --fd-mode=" + vinfo.fd_mode;
        }
        if (vinfo.timeout > 0.0) {
            cmd += " --timeout=" + std::to_string(vinfo.timeout);
        }
        RunResult res = run_ibex_base(cmd, vinfo.name, idx);
        lock_guard<mutex> lk(m);
        ofstream o(csv, ios::app);
        ofstream ot(trig_csv, ios::app);
        if (!res.csv_line.empty()) {
            string line = res.csv_line;
            string prob_name = fs::path(problem).filename().string();
            size_t first = line.find(',');
            size_t second = line.find(',', first+1);
            if (first!=string::npos && second!=string::npos) {
                line.insert(second, "," + prob_name);
            }
            o << line;
            if (line.back() != '\n') o << "\n";
            ot << idx << "," << vinfo.name << "," << prob_name << ","
               << (res.triggers_total>=0?res.triggers_total:0) << ","
               << (res.triggers_depth>=0?res.triggers_depth:0) << ","
               << (res.triggers_vol>=0?res.triggers_vol:0) << ","
               << (res.vol_eval>=0?res.vol_eval:0) << ","
               << (res.vol_nonfinite>=0?res.vol_nonfinite:0) << "\n";
        } else {
            o << idx << "," << vinfo.name << "," << fs::path(problem).filename().string()
              << ",NA,NA,NA,NA,NA,0,0,0,0,0,0\n";
            ot << idx << "," << vinfo.name << "," << fs::path(problem).filename().string() << ",0,0,0,0,0\n";
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
