#include "quartz/context/context.h"
#include "quartz/gate/gate_utils.h"
#include "quartz/parser/qasm_parser.h"
#include "quartz/pybind/pybind.h"
#include "quartz/simulator/schedule.h"
#include "quartz/tasograph/tasograph.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <future>
#include <iostream>
#include <signal.h>
#include <stdexcept>
#include <sys/types.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <vector>

using namespace quartz;

// Timeout for ILP in seconds (change this value for testing)
constexpr int ILP_TIMEOUT_SECONDS = 10;  // e.g., 10 seconds for testing

// Temporary file for process-based result passing
const char *ILP_RESULT_FILE = "ilp_result.tmp";

// const std::string QASM_FILE_PREFIX = "./circuit/qiskit-random/rqc_";
const std::string QASM_FILE_PREFIX = "./circuit/qiskit-random/rqc_depth_1";

int num_stages_by_heuristics(CircuitSeq *seq, int num_local_qubits,
                             std::vector<std::vector<bool>> &local_qubits,
                             int &num_swaps) {
  int num_qubits = seq->get_num_qubits();
  std::unordered_map<CircuitGate *, bool> executed;
  // No initial configuration -- all qubits are global.
  std::vector<bool> local_qubit(num_qubits, false);
  int num_stages = 0;
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
  return num_stages;
}

// Helper function to serialize result to file
void write_ilp_result(const std::vector<std::vector<int>> &result) {
  std::ofstream ofs(ILP_RESULT_FILE, std::ios::binary | std::ios::trunc);
  int n = result.size();
  ofs.write(reinterpret_cast<const char *>(&n), sizeof(int));
  for (const auto &vec : result) {
    int m = vec.size();
    ofs.write(reinterpret_cast<const char *>(&m), sizeof(int));
    ofs.write(reinterpret_cast<const char *>(vec.data()), m * sizeof(int));
  }
  ofs.close();
}

// Helper function to deserialize result from file
std::vector<std::vector<int>> read_ilp_result() {
  std::vector<std::vector<int>> result;
  std::ifstream ifs(ILP_RESULT_FILE, std::ios::binary);
  if (!ifs)
    return result;
  int n = 0;
  ifs.read(reinterpret_cast<char *>(&n), sizeof(int));
  for (int i = 0; i < n; ++i) {
    int m = 0;
    ifs.read(reinterpret_cast<char *>(&m), sizeof(int));
    std::vector<int> vec(m);
    ifs.read(reinterpret_cast<char *>(vec.data()), m * sizeof(int));
    result.push_back(std::move(vec));
  }
  ifs.close();
  return result;
}

// Helper function to run compute_qubit_layout_with_ilp in a separate process
template <typename... Args>
std::optional<std::vector<std::vector<int>>>
run_with_timeout_process(std::chrono::seconds timeout, Args &&...args) {
  // Remove old result file if exists
  std::remove(ILP_RESULT_FILE);

  pid_t pid = fork();
  if (pid < 0) {
    std::cerr << "fork() failed!" << std::endl;
    return std::nullopt;
  }
  if (pid == 0) {
    // Child process
    auto result = compute_qubit_layout_with_ilp(std::forward<Args>(args)...);
    write_ilp_result(result);
    _exit(0);
  } else {
    // Parent process
    int status = 0;
    int waited = 0;
    while (waited < timeout.count()) {
      pid_t result = waitpid(pid, &status, WNOHANG);
      if (result == 0) {
        sleep(1);
        ++waited;
      } else {
        break;
      }
    }
    if (waited >= timeout.count()) {
      kill(pid, SIGKILL);
      waitpid(pid, &status, 0);
      return std::nullopt;
    } else {
      // Read result from file
      auto result = read_ilp_result();
      return result;
    }
  }
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
  // std::vector<int> num_qubits = {30, 32};
  int num_local_qubits = 28;

  // Test staging by heuristics
  for (int num_q : num_qubits) {
    // requires running test_remove_swap first
    auto seq = CircuitSeq::from_qasm_file(
        &ctx, (QASM_FILE_PREFIX + "_" + std::to_string(num_q) + ".qasm"));

    fprintf(fout, "%d, ", num_q);
    std::vector<int> n_swaps;
    std::vector<std::vector<bool>> local_qubits_by_heuristics;
    int num_swaps = 0;
    int heuristics_result = num_stages_by_heuristics(
        seq.get(), num_local_qubits, local_qubits_by_heuristics, num_swaps);
    n_swaps.push_back(num_swaps);
    fprintf(fout, "%d, ", heuristics_result);
    fflush(fout);
    fprintf(fout, "\n");
    fflush(fout);
  }

  // Test hyper-staging by heuristics
  for (int num_q : num_qubits) {
    // requires running test_remove_swap first
    auto seq = CircuitSeq::from_qasm_file(
        &ctx, (QASM_FILE_PREFIX + "_" + std::to_string(num_q) + ".qasm"));

    fprintf(fout, "%d, ", num_q);
    std::vector<int> n_swaps;
    std::vector<std::vector<bool>> local_qubits_by_heuristics;
    int num_swaps = 0;
    int heuristics_result = num_stages_by_heuristics(
        seq.get(), num_q - 1, local_qubits_by_heuristics, num_swaps);
    n_swaps.push_back(num_swaps);
    fprintf(fout, "%d, ", heuristics_result);
    fflush(fout);
    fprintf(fout, "\n");
    fflush(fout);
  }

  // Test staging by ILP
  for (int num_q : num_qubits) {
    auto seq = CircuitSeq::from_qasm_file(
        &ctx, (QASM_FILE_PREFIX + "_" + std::to_string(num_q) + ".qasm"));

    fprintf(fout, "%d, ", num_q);
    int answer_start_with = 1;
    std::vector<int> n_swaps;
    std::vector<std::vector<int>> local_qubits;
    int num_global_q = num_q - num_local_qubits;

    // Run compute_qubit_layout_with_ilp in a separate process with a timeout
    std::optional<std::vector<std::vector<int>>> ilp_result_opt =
        run_with_timeout_process(std::chrono::seconds(ILP_TIMEOUT_SECONDS),
                                 *seq, num_local_qubits,
                                 std::min(2, num_q - num_local_qubits), &ctx,
                                 &interpreter, answer_start_with);

    int ilp_result = -1;
    int num_swaps = 0;
    if (ilp_result_opt.has_value() && !ilp_result_opt.value().empty()) {
      local_qubits = ilp_result_opt.value();
      ilp_result = (int)local_qubits.size();
      std::vector<bool> prev_local(num_q, false);
      for (int j = 0; j < ilp_result; j++) {
        std::cout << "Stage " << j << ": ";
        local_qubits[j].resize(num_q - num_global_q);
        for (int k : local_qubits[j]) {
          std::cout << k << " ";
          if (j > 0) {
            if (!prev_local[k]) {
              num_swaps++;
            }
          }
        }
        prev_local.assign(num_q, false);
        for (int k : local_qubits[j]) {
          prev_local[k] = true;
        }
        std::cout << std::endl;
      }
    } else {
      std::cerr << "Error: compute_qubit_layout_with_ilp exceeded timeout ("
                << ILP_TIMEOUT_SECONDS << " seconds)." << std::endl;
      ilp_result = -1;
    }
    n_swaps.push_back(num_swaps);
    fprintf(fout, "%d, ", ilp_result);
    fflush(fout);
    answer_start_with = ilp_result;
    fprintf(fout, "\n");
    fflush(fout);
  }

  // Test hyper-staging by ILP
  for (int num_q : num_qubits) {
    auto seq = CircuitSeq::from_qasm_file(
        &ctx, (QASM_FILE_PREFIX + "_" + std::to_string(num_q) + ".qasm"));
    fprintf(fout, "%d, ", num_q);
    int answer_start_with = 1;
    std::vector<int> n_swaps;
    std::vector<std::vector<int>> local_qubits;
    int num_global_q = num_q - num_local_qubits;

    // Run compute_qubit_layout_with_ilp in a separate process with a timeout
    std::optional<std::vector<std::vector<int>>> ilp_result_opt =
        run_with_timeout_process(std::chrono::seconds(ILP_TIMEOUT_SECONDS),
                                 *seq, num_q - 1, 0, &ctx, &interpreter,
                                 answer_start_with);

    int ilp_result = -1;
    int num_swaps = 0;
    if (ilp_result_opt.has_value() && !ilp_result_opt.value().empty()) {
      local_qubits = ilp_result_opt.value();
      ilp_result = (int)local_qubits.size();
      std::vector<bool> prev_local(num_q, false);
      for (int j = 0; j < ilp_result; j++) {
        std::cout << "Stage " << j << ": ";
        local_qubits[j].resize(num_q - num_global_q);
        for (int k : local_qubits[j]) {  // print the local qubits
          std::cout << k << " ";
          if (j > 0) {
            if (!prev_local[k]) {
              num_swaps++;
            }
          }
        }
        prev_local.assign(num_q, false);
        for (int k : local_qubits[j]) {
          prev_local[k] = true;
        }
        std::cout << std::endl;
      }
    } else {
      std::cerr << "Error: compute_qubit_layout_with_ilp exceeded timeout ("
                << ILP_TIMEOUT_SECONDS << " seconds)." << std::endl;
      ilp_result = -1;
    }
    n_swaps.push_back(num_swaps);
    fprintf(fout, "%d, ", ilp_result);
    fflush(fout);
    answer_start_with = ilp_result;
    fprintf(fout, "\n");
    fflush(fout);
  }

  fclose(fout);
  // Remove temporary result file if exists
  std::remove(ILP_RESULT_FILE);

  auto end = std::chrono::steady_clock::now();
  std::cout
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " seconds." << std::endl;
  return 0;
}
