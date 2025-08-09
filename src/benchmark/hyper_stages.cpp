#include "quartz/context/context.h"
#include "quartz/gate/gate_utils.h"
#include "quartz/pybind/pybind.h"
#include "quartz/simulator/schedule.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace quartz;

void get_stages_by_heuristics(
    CircuitSeq *seq, int num_local_qubits,
    std::vector<std::vector<bool>> &local_qubits, int &num_swaps,
    std::unordered_set<CircuitGate *> &gates_in_hyperstage) {
  int num_qubits = seq->get_num_qubits();
  std::unordered_map<CircuitGate *, bool> executed;
  // No initial configuration -- all qubits are global.
  std::vector<bool> local_qubit(num_qubits, false);
  int num_stages = 0;
  while (true) {
    bool all_done = true;
    std::vector<bool> executable(num_qubits, true);
    for (auto &gate : seq->gates) {
      if (gates_in_hyperstage.count(gate.get()) == 0) {
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
    if (all_done) {
      break;
    }
    num_stages++;
    // count global and local gates
    std::vector<bool> first_unexecuted_gate(
        num_qubits, false);  // for tiebreaker, means that the first unexecuted
                             // gate on this qubit
    // how many local and global gates on each qubit
    std::vector<int> local_gates(num_qubits, 0);
    std::vector<int> global_gates(num_qubits, 0);
    // how many local and global gates on each qubit

    bool first = true;
    for (auto &gate : seq->gates) {
      if (gate->gate->is_quantum_gate() && !executed[gate.get()]) {
        bool local = true;
        if (!gate->gate->is_diagonal()) {
          int num_remaining_control_qubits =
              gate->gate->get_num_control_qubits();
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
      // Use the qubit index as a final tiebreaker.
      return a < b;
    };
    std::vector<int> candidate_indices(num_qubits, 0);
    for (int i = 0; i < num_qubits; i++) {
      candidate_indices[i] = i;
      local_qubit[i] = false;
    }
    std::sort(candidate_indices.begin(), candidate_indices.end(), cmp);
    std::cout << "Stage " << num_stages << ": {";
    for (int i = 0; i < num_local_qubits; i++) {
      local_qubit[candidate_indices[i]] = true;
      std::cout << candidate_indices[i];
      if (i < num_local_qubits - 1) {
        std::cout << ", ";
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
  // No initial configuration -- all qubits are global.
  std::vector<bool> local_qubit(num_qubits, false);
  int num_stages = 0;
  int iter = 0;
  std::unordered_set<CircuitGate *> newly_executed_gates;
  while (true) {
    bool all_done = true;
    std::vector<bool> executable(num_qubits, true);
    for (auto &gate : seq->gates) {
      if (gate->gate->is_quantum_gate() && !executed[gate.get()]) {
        bool ok = true;  // indicates if the gate can be executed
        for (auto &output : gate->output_wires) {
          if (!executable[output->index]) {
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
          // execute
          newly_executed_gates.insert(gate.get());
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
    if (iter > 0) {  // the first and second iteration generate the final list
                     // of gates that are executed by in the first stage
      std::unordered_set<CircuitGate *> newly_executed_gates_copy =
          newly_executed_gates;
      executed_gates_per_stage.push_back(std::move(newly_executed_gates_copy));
      newly_executed_gates.clear();
    }
    if (all_done) {
      break;
    }
    num_stages++;
    // count global and local gates
    std::vector<bool> first_unexecuted_gate(
        num_qubits, false);  // for tiebreaker, means that the first unexecuted
                             // gate on this qubit
    // how many local and global gates on each qubit
    std::vector<int> local_gates(num_qubits, 0);
    std::vector<int> global_gates(num_qubits, 0);
    // how many local and global gates on each qubit

    bool first = true;
    for (auto &gate : seq->gates) {
      if (gate->gate->is_quantum_gate() && !executed[gate.get()]) {
        bool local = true;
        if (!gate->gate->is_diagonal()) {
          int num_remaining_control_qubits =
              gate->gate->get_num_control_qubits();
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
      // Use the qubit index as a final tiebreaker.
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
    for (int i = 0; i < num_local_qubits; i++) {
      local_qubit[candidate_indices[i]] = true;
      std::cout << candidate_indices[i];
      if (i < num_local_qubits - 1) {
        std::cout << ", ";
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

int main() {
  auto start = std::chrono::steady_clock::now();
  init_python_interpreter();
  PythonInterpreter interpreter;
  Context ctx({GateType::input_qubit, GateType::input_param, GateType::h,
               GateType::x, GateType::ry, GateType::u2, GateType::u3,
               GateType::cx, GateType::cz, GateType::cp, GateType::swap,
               GateType::rz, GateType::p, GateType::ccx, GateType::rx});
  FILE *fout = fopen("heuristic_result.csv", "w");
  // 31 or 42 total qubits, 0-23 global qubits
  // std::vector<int> num_qubits = {28, 29, 28, 29, 31, 32, 33,
  std::vector<int> num_qubits = {30, 32, 34, 36, 38, 40, 42, 44, 46, 48};
  int num_frozen_qubits = 1;
  constexpr int kMaxGlobalQubitsFor31 = 16;
  std::vector<int> num_global_qubits;
  for (int i = 0; i <= 24; i++) {
    num_global_qubits.push_back(i);
  }
  for (int num_q : num_qubits) {
    // requires running test_remove_swap first
    auto seq = CircuitSeq::from_qasm_file(
        &ctx, (std::string("./circuit/qiskit-random/rqc") + "_" +
               std::to_string(num_q) + ".qasm"));

    fprintf(fout, "%d, ", num_q);
    std::vector<int> n_swaps;
    std::vector<std::vector<bool>> local_qubits_by_heuristics;
    int num_swaps = 0;
    std::vector<std::vector<int>> heuristics_result;
    std::vector<std::unordered_set<CircuitGate *>> executed_gates_per_stage;
    get_hyper_stages(seq.get(), num_frozen_qubits, local_qubits_by_heuristics,
                     num_swaps, heuristics_result, executed_gates_per_stage);
    std::cout << "heuristics_result.size(): " << heuristics_result.size()
              << std::endl;
    std::cout << "executed_gates_per_stage.size(): "
              << executed_gates_per_stage.size() << std::endl;
    n_swaps.push_back(num_swaps);
    fprintf(fout, "%d, ", (int)heuristics_result.size());
    fflush(fout);
    fprintf(fout, "\n");
    fflush(fout);
  }

  fclose(fout);
  auto end = std::chrono::steady_clock::now();
  std::cout
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " seconds." << std::endl;
  return 0;
}
