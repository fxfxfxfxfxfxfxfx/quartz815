#include "quartz/context/context.h"
#include "quartz/gate/gate_utils.h"
#include "quartz/pybind/pybind.h"
#include "quartz/simulator/schedule.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace quartz;

// Function to parse the result of
// compute_qubit_layout_with_hyper_stage_heuristic and count the number of hyper
// stages by detecting changes in the last qubit id.
int count_hyper_stages_from_layout(
    const std::vector<std::vector<int>> &layout) {
  if (layout.empty())
    return 0;
  int count = 1;
  int last_qubit_id = layout[0].back();
  for (size_t i = 1; i < layout.size(); ++i) {
    int current_qubit_id = layout[i].back();
    if (current_qubit_id != last_qubit_id) {
      ++count;
      last_qubit_id = current_qubit_id;
    }
  }
  return count;
}

// Extracted helper function to compute gate statistics and first unexecuted
// gate
static void compute_gate_statistics(
    CircuitSeq *seq, const std::unordered_map<CircuitGate *, bool> &executed,
    const std::vector<bool> &local_qubit, std::vector<int> &local_gates,
    std::vector<int> &global_gates, std::vector<bool> &first_unexecuted_gate) {
  int num_qubits = seq->get_num_qubits();
  local_gates.assign(num_qubits, 0);
  global_gates.assign(num_qubits, 0);
  first_unexecuted_gate.assign(num_qubits, false);

  bool first = true;
  for (auto &gate : seq->gates) {
    if (gate->gate->is_quantum_gate() && !executed.at(gate.get())) {
      bool local = true;
      if (!gate->gate->is_diagonal()) {
        int num_remaining_control_qubits = gate->gate->get_num_control_qubits();
        for (auto &output : gate->output_wires) {
          if (output->is_qubit()) {
            num_remaining_control_qubits--;
            if (num_remaining_control_qubits < 0 &&
                !local_qubit[output->index]) {
              local = false;
            }
          }
        }
      }
      int num_remaining_control_qubits = gate->gate->get_num_control_qubits();
      for (auto &output : gate->output_wires) {
        if (output->is_qubit()) {
          num_remaining_control_qubits--;
          if (local) {
            local_gates[output->index]++;
          } else {
            global_gates[output->index]++;
          }
          if (first && num_remaining_control_qubits < 0) {
            first_unexecuted_gate[output->index] = true;
          }
        }
      }
      first = false;
    }
  }
}

// Extracted helper function to update executed map and executable vector
static bool update_executed_and_executable(
    CircuitSeq *seq,
    const std::unordered_set<CircuitGate *>
        *gates_in_hyperstage,  // can be nullptr
    std::unordered_map<CircuitGate *, bool> &executed,
    std::vector<bool> &executable, const std::vector<bool> &local_qubit) {
  bool all_done = true;
  for (auto &gate : seq->gates) {
    if (gates_in_hyperstage && gates_in_hyperstage->count(gate.get()) == 0) {
      continue;  // skip gates not in the hyperstage
    }
    if (gate->gate->is_quantum_gate() && !executed[gate.get()]) {
      bool ok = true;  // indicates if the gate can be executed
      for (auto &output : gate->output_wires) {
        if (!executable[output->index]) {
          ok = false;
        }
      }
      if (!gate->gate->is_diagonal()) {
        int num_remaining_control_qubits = gate->gate->get_num_control_qubits();
        for (auto &output : gate->output_wires) {
          if (output->is_qubit()) {
            num_remaining_control_qubits--;
            if (num_remaining_control_qubits < 0 &&
                !local_qubit[output->index]) {
              ok = false;
            }
          }
        }
      }
      if (ok) {
        // execute
        executed[gate.get()] = true;
      } else {
        // not executable, block the qubits
        all_done = false;
        for (auto &output : gate->output_wires) {
          executable[output->index] = false;
        }
      }
    }
  }
  return all_done;
}

void get_stages_by_heuristics(
    CircuitSeq *seq, int num_local_qubits,
    std::vector<std::vector<bool>> &local_qubits, int &num_swaps,
    std::unordered_set<CircuitGate *> &gates_in_hyperstage,
    std::vector<std::vector<int>> &res, std::unordered_set<int> &frozen_qubits,
    std::unordered_set<int> &prev_frozen_qubits) {
  res.clear();
  int num_qubits = seq->get_num_qubits();
  std::unordered_map<CircuitGate *, bool> executed;
  for (auto &gate : seq->gates) {
    executed[gate.get()] = false;
  }
  std::vector<bool> local_qubit(num_qubits, false);
  int num_stages = 0;
  int iter = 0;
  while (true) {
    std::vector<bool> executable(num_qubits, true);
    bool all_done = update_executed_and_executable(
        seq, &gates_in_hyperstage, executed, executable, local_qubit);
    if (all_done) {
      break;
    }
    num_stages++;
    std::vector<bool> first_unexecuted_gate;
    std::vector<int> local_gates, global_gates;
    compute_gate_statistics(seq, executed, local_qubit, local_gates,
                            global_gates, first_unexecuted_gate);

    auto cmp = [&](int a, int b) {
      if (iter == 0 && prev_frozen_qubits.count(a) &&
          prev_frozen_qubits.count(b)) {
        return a < b;  // both are frozen, use index as tiebreaker
      }
      if (iter == 0 && prev_frozen_qubits.count(a))
        return true;  // a is frozen, b is not, a comes first
      if (iter == 0 && prev_frozen_qubits.count(b))
        return false;  // b is frozen, a is not, b comes first

      if (frozen_qubits.count(a) && frozen_qubits.count(b))
        return a < b;  // both are frozen, use index as tiebreaker
      if (frozen_qubits.count(a))
        return false;  // a is frozen, b is not, b comes first
      if (frozen_qubits.count(b))
        return true;  // b is frozen, a is not, a comes first
      if (first_unexecuted_gate[b])
        return false;
      if (first_unexecuted_gate[a])
        return true;
      if (global_gates[a] != global_gates[b]) {
        return global_gates[a] > global_gates[b];
      }
      if (local_gates[a] != local_gates[b]) {
        return local_gates[a] > local_gates[b];
      }
      return a < b;
    };
    std::vector<int> candidate_indices(num_qubits, 0);
    for (int i = 0; i < num_qubits; i++) {
      candidate_indices[i] = i;
      local_qubit[i] = false;
    }
    std::sort(candidate_indices.begin(), candidate_indices.end(), cmp);
    res.push_back(candidate_indices);
    std::cout << "Stage " << num_stages << ": {";
    for (int i = 0; i < num_qubits; i++) {
      std::cout << candidate_indices[i];
      if (i < num_qubits - 1) {
        std::cout << ", ";
      }
      if (i < num_local_qubits) {
        local_qubit[candidate_indices[i]] = true;
      }
    }
    std::cout << "}" << std::endl;
    local_qubits.push_back(local_qubit);
    if (num_stages != 1) {
      for (int i = 0; i < num_local_qubits; i++) {
        if (!local_qubits[num_stages - 2][candidate_indices[i]]) {
          num_swaps++;
        }
      }
    }
    iter++;
  }
  std::cout << num_stages << " stages." << std::endl;
}

void get_hyper_stages(
    CircuitSeq *seq, int num_frozen_qubits,
    std::vector<std::vector<bool>> &local_qubits, int &num_swaps,
    std::vector<std::vector<int>> &result,
    std::vector<std::unordered_set<CircuitGate *>> &executed_gates_per_stage) {
  result.clear();
  executed_gates_per_stage.clear();
  int num_qubits = seq->get_num_qubits();
  int num_local_qubits = num_qubits - num_frozen_qubits;
  std::unordered_map<CircuitGate *, bool> executed;
  for (auto &gate : seq->gates) {
    executed[gate.get()] = false;
  }
  std::vector<bool> local_qubit(num_qubits, false);
  int num_stages = 0;
  while (true) {
    std::vector<bool> executable(num_qubits, true);
    bool all_done = update_executed_and_executable(seq, nullptr, executed,
                                                   executable, local_qubit);
    if (all_done) {
      break;
    }
    num_stages++;
    std::vector<bool> first_unexecuted_gate;
    std::vector<int> local_gates, global_gates;
    compute_gate_statistics(seq, executed, local_qubit, local_gates,
                            global_gates, first_unexecuted_gate);

    auto cmp = [&](int a, int b) {
      if (first_unexecuted_gate[b])
        return false;
      if (first_unexecuted_gate[a])
        return true;
      if (global_gates[a] != global_gates[b]) {
        return global_gates[a] > global_gates[b];
      }
      if (local_gates[a] != local_gates[b]) {
        return local_gates[a] > local_gates[b];
      }
      return a < b;
    };
    std::vector<int> candidate_indices(num_qubits, 0);
    for (int i = 0; i < num_qubits; i++) {
      candidate_indices[i] = i;
      local_qubit[i] = false;
    }
    std::sort(candidate_indices.begin(), candidate_indices.end(), cmp);
    result.push_back(candidate_indices);
    std::cout << "Stage " << num_stages << ": {";
    for (int i = 0; i < num_qubits; i++) {
      std::cout << candidate_indices[i];
      if (i < num_qubits - 1) {
        std::cout << ", ";
      }
      if (i < num_local_qubits) {
        local_qubit[candidate_indices[i]] = true;
      }
    }
    std::cout << "}" << std::endl;

    // Collect executable gates in this hyper stage
    std::vector<bool> executable_for_collect_gates(num_qubits, true);
    std::unordered_set<CircuitGate *> executed_this_stage;
    for (auto &gate : seq->gates) {
      if (gate->gate->is_quantum_gate() && !executed[gate.get()]) {
        bool ok = true;
        for (auto &output : gate->output_wires) {
          if (!executable_for_collect_gates[output->index]) {
            ok = false;
          }
        }
        if (!gate->gate->is_diagonal()) {
          int num_remaining_control_qubits =
              gate->gate->get_num_control_qubits();
          for (auto &output : gate->output_wires) {
            if (output->is_qubit()) {
              num_remaining_control_qubits--;
              if (num_remaining_control_qubits < 0 &&
                  !local_qubit[output->index]) {
                ok = false;
              }
            }
          }
        }
        if (ok) {
          executed_this_stage.insert(gate.get());
        } else {
          for (auto &output : gate->output_wires) {
            executable_for_collect_gates[output->index] = false;
          }
        }
      }
    }
    if (executed_this_stage.empty()) {
      std::cout << "executed_this_stage.size(): "
                << (int)executed_this_stage.size() << std::endl;
    }
    executed_gates_per_stage.push_back(std::move(executed_this_stage));

    local_qubits.push_back(local_qubit);
    if (num_stages != 1) {
      for (int i = 0; i < num_local_qubits; i++) {
        if (!local_qubits[num_stages - 2][candidate_indices[i]]) {
          num_swaps++;
        }
      }
    }
  }
  std::cout << num_stages << " stages." << std::endl;
}

// Extracted function as requested
// NOTE: This function is deprecated and have some bugs
// should use compute_qubit_layout_with_hyper_stage_heuristic in schedule.cpp
// instead
void get_hyper_q_stages(CircuitSeq *seq, int num_frozen_qubits, int num_q,
                        int NUM_LOCAL_QUBITS, FILE *fout) {
  std::vector<std::vector<bool>> local_qubits_by_heuristics;
  int num_swaps = 0;
  std::vector<std::vector<int>> hyper_stages;
  std::vector<std::unordered_set<CircuitGate *>> executed_gates_per_stage;

  // first get hyper stages
  get_hyper_stages(seq, num_frozen_qubits, local_qubits_by_heuristics,
                   num_swaps, hyper_stages, executed_gates_per_stage);
  std::cout << "heuristics_result.size(): " << hyper_stages.size() << std::endl;
  std::cout << "executed_gates_per_stage.size(): "
            << executed_gates_per_stage.size() << std::endl;
  std::vector<int> n_swaps;
  n_swaps.push_back(num_swaps);

  // then get stages in each hyper stage
  std::unordered_set<int> prev_frozen_qubits;
  int num_stages = 0;
  for (int i = 0; i < hyper_stages.size(); i++) {
    std::vector<std::vector<int>> stages;
    std::vector<std::vector<bool>> local_qubits;
    int num_swaps = 0;
    // Get the last `frozen_qubits` qubits from heuristics_result[i]
    std::unordered_set<int> frozen_qubits;
    for (int j = 0; j < num_frozen_qubits; j++) {
      frozen_qubits.insert(hyper_stages[i][num_q - 1 - j]);
    }
    get_stages_by_heuristics(seq, NUM_LOCAL_QUBITS, local_qubits, num_swaps,
                             executed_gates_per_stage[i], stages, frozen_qubits,
                             prev_frozen_qubits);
    prev_frozen_qubits = frozen_qubits;
    num_stages += stages.size();
  }

  // Write to file and stdout
  fprintf(fout, "%d, %d", (int)hyper_stages.size(), num_stages);
  fflush(fout);
  std::cout << (int)hyper_stages.size() << ", " << num_stages;
  std::cout << std::endl;
  fprintf(fout, "\n");
  fflush(fout);
}

void run_hyper_stage_partition(CircuitSeq *seq, int num_frozen_qubits,
                               int num_q, int NUM_LOCAL_QUBITS, FILE *fout,
                               Context *ctx) {
  std::vector<std::vector<bool>> local_qubits_by_heuristics;
  std::vector<std::vector<int>> hyper_stages;
  std::vector<std::unordered_set<CircuitGate *>> executed_gates_per_stage;

  // first get hyper stages
  auto res = compute_qubit_layout_with_hyper_stage_heuristic(
      *seq, NUM_LOCAL_QUBITS, 1, ctx);

  // Use the new function to count hyper stages
  int hyper_stage_count = count_hyper_stages_from_layout(res);
  fprintf(fout, "%d, %d", hyper_stage_count, (int)res.size());
  fflush(fout);
  std::cout << hyper_stage_count << ", " << (int)res.size();
  std::cout << std::endl;
  fprintf(fout, "\n");
  fflush(fout);
}

void print_usage(const char *prog_name) {
  std::cout << "Usage: " << prog_name
            << " -q <qasm_file> -l <num_local_qubits> -f <num_frozen_qubits>\n";
  std::cout << "  -q <qasm_file>           Path to QASM file\n";
  std::cout << "  -l <num_local_qubits>    Number of local qubits\n";
  std::cout << "  -f <num_frozen_qubits>   Number of frozen qubits\n";
  std::cout << "  -h                       Show this help message\n";
}

int main(int argc, char *argv[]) {
  auto start = std::chrono::steady_clock::now();
  init_python_interpreter();
  PythonInterpreter interpreter;
  Context ctx({GateType::input_qubit, GateType::input_param, GateType::h,
               GateType::x, GateType::ry, GateType::u2, GateType::u3,
               GateType::cx, GateType::cz, GateType::cp, GateType::swap,
               GateType::rz, GateType::p, GateType::ccx, GateType::rx});

  std::string qasm_file = "";
  int num_local_qubits = -1;
  int num_frozen_qubits = -1;

  // Parse command line arguments
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-q") == 0 && i + 1 < argc) {
      qasm_file = argv[++i];
    } else if (strcmp(argv[i], "-l") == 0 && i + 1 < argc) {
      num_local_qubits = std::stoi(argv[++i]);
    } else if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
      num_frozen_qubits = std::stoi(argv[++i]);
    } else if (strcmp(argv[i], "-h") == 0) {
      print_usage(argv[0]);
      return 0;
    } else {
      std::cerr << "Unknown or incomplete argument: " << argv[i] << std::endl;
      print_usage(argv[0]);
      return 1;
    }
  }

  if (qasm_file.empty() || num_local_qubits < 0 || num_frozen_qubits < 0) {
    std::cerr << "Missing required arguments.\n";
    print_usage(argv[0]);
    return 1;
  }

  FILE *fout = fopen("hyper_stage_result.csv", "w");
  if (!fout) {
    std::cerr << "Failed to open output file.\n";
    return 1;
  }

  // Load circuit from QASM file
  auto seq = CircuitSeq::from_qasm_file(&ctx, qasm_file);
  if (!seq) {
    std::cerr << "Failed to load QASM file: " << qasm_file << std::endl;
    fclose(fout);
    return 1;
  }

  int num_q = seq->get_num_qubits();

  fprintf(fout, "%d, ", num_q);
  std::cout << num_q << ", ";

  // Call the extracted function
  run_hyper_stage_partition(seq.get(), num_frozen_qubits, num_q,
                            num_local_qubits, fout, &ctx);

  auto res = compute_qubit_layout_with_snuqs_heuristic(
      *seq, num_local_qubits + 1, num_frozen_qubits, &ctx);
  fprintf(fout, "%d, ", num_q);
  fprintf(fout, "%d, ", (int)res.size());
  fflush(fout);
  std::cout << num_q << ", ";
  std::cout << (int)res.size() << ", ";
  std::cout << std::endl;
  fprintf(fout, "\n");
  fflush(fout);

  fclose(fout);
  auto end = std::chrono::steady_clock::now();
  std::cout
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " seconds." << std::endl;
  return 0;
}
