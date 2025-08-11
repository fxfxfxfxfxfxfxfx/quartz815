#include <algorithm>
#include <climits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Extract global qubits from a stage
std::vector<int> get_global_qubits(const std::vector<int>& stage,
                                   int num_local_qubits) {
  return std::vector<int>(stage.begin() + num_local_qubits, stage.end());
}

// Extract local qubits from a stage
std::vector<int> get_local_qubits(const std::vector<int>& stage,
                                  int num_local_qubits) {
  return std::vector<int>(stage.begin(), stage.begin() + num_local_qubits);
}

// Find common elements between two vectors
std::vector<int> find_common_elements(const std::vector<int>& vec1,
                                      const std::vector<int>& vec2) {
  std::unordered_set<int> set1(vec1.begin(), vec1.end());
  std::vector<int> common;

  for (int elem : vec2) {
    if (set1.count(elem)) {
      common.push_back(elem);
    }
  }
  return common;
}

// Find all possible frozen qubits for a range of stages
std::vector<int>
find_common_global_qubits(const std::vector<std::vector<int>>& stages,
                          int start_idx, int end_idx, int num_local_qubits) {
  if (start_idx > end_idx)
    return {};

  // Find qubits that appear in ALL stages' global qubits in this range
  std::unordered_map<int, int> frequency;
  int range_size = end_idx - start_idx + 1;

  for (int i = start_idx; i <= end_idx; i++) {
    auto global_qubits = get_global_qubits(stages[i], num_local_qubits);
    std::unordered_set<int> stage_qubits(global_qubits.begin(),
                                         global_qubits.end());
    for (int qubit : stage_qubits) {
      frequency[qubit]++;
    }
  }

  std::vector<int> common_qubits;
  for (const auto& pair : frequency) {
    if (pair.second == range_size) {  // Appears in ALL stages
      common_qubits.push_back(pair.first);
    }
  }

  return common_qubits;
}

// Calculate transition score: how many frozen qubits can become local qubits in
// next hyper-stage
int calculate_transition_score(const std::vector<int>& frozen_qubits,
                               const std::vector<std::vector<int>>& stages,
                               int next_hyper_start, int num_local_qubits) {
  if (next_hyper_start >= stages.size())
    return 0;

  auto next_local_qubits =
      get_local_qubits(stages[next_hyper_start], num_local_qubits);
  std::unordered_set<int> next_local_set(next_local_qubits.begin(),
                                         next_local_qubits.end());

  int score = 0;
  for (int qubit : frozen_qubits) {
    if (next_local_set.count(qubit)) {
      score++;
    }
  }
  return score;
}

// Find optimal frozen qubits considering both objectives
// Returns empty vector if not enough common qubits available
std::vector<int> find_optimal_frozen_qubits_with_transition(
    const std::vector<std::vector<int>>& stages, int start_idx, int end_idx,
    int next_hyper_start, int num_local_qubits, int num_frozen_qubits) {
  // Get all possible frozen qubits (qubits that appear in all stages of this
  // range)
  std::vector<int> candidates =
      find_common_global_qubits(stages, start_idx, end_idx, num_local_qubits);

  // Check if we have enough candidates to form a valid hyper-stage
  if (candidates.size() < num_frozen_qubits) {
    return {};  // Invalid hyper-stage - not enough common qubits
  }

  // If we have exactly the right number, use all candidates
  if (candidates.size() == num_frozen_qubits) {
    return candidates;
  }

  // If we have more candidates than needed, choose the best combination
  // considering transition to next hyper-stage
  std::vector<int> best_frozen;

  if (next_hyper_start < stages.size()) {
    auto next_local_qubits =
        get_local_qubits(stages[next_hyper_start], num_local_qubits);
    std::unordered_set<int> next_local_set(next_local_qubits.begin(),
                                           next_local_qubits.end());

    // Separate candidates into two groups: helpful for transition and others
    std::vector<int> transition_helpful;
    std::vector<int> others;

    for (int qubit : candidates) {
      if (next_local_set.count(qubit)) {
        transition_helpful.push_back(qubit);
      } else {
        others.push_back(qubit);
      }
    }

    // Greedily select: first prioritize transition-helpful qubits
    int helpful_to_use =
        std::min((int)transition_helpful.size(), num_frozen_qubits);
    for (int i = 0; i < helpful_to_use; i++) {
      best_frozen.push_back(transition_helpful[i]);
    }

    // Fill remaining slots with other candidates
    int remaining = num_frozen_qubits - helpful_to_use;
    for (int i = 0; i < remaining && i < others.size(); i++) {
      best_frozen.push_back(others[i]);
    }
  } else {
    // No next hyper-stage, just take first num_frozen_qubits
    best_frozen = std::vector<int>(candidates.begin(),
                                   candidates.begin() + num_frozen_qubits);
  }

  return best_frozen;
}

// Reorder stage to put frozen qubits at the end of global qubits
std::vector<int>
reorder_stage_with_frozen_at_end(const std::vector<int>& stage,
                                 const std::vector<int>& frozen_qubits,
                                 int num_local_qubits) {
  auto local_qubits = get_local_qubits(stage, num_local_qubits);
  auto global_qubits = get_global_qubits(stage, num_local_qubits);

  std::vector<int> reordered_stage = local_qubits;
  std::unordered_set<int> frozen_set(frozen_qubits.begin(),
                                     frozen_qubits.end());

  // First, add non-frozen global qubits
  for (int qubit : global_qubits) {
    if (frozen_set.find(qubit) == frozen_set.end()) {
      reordered_stage.push_back(qubit);
    }
  }

  // Then, add frozen qubits at the end
  for (int qubit : frozen_qubits) {
    // Make sure this frozen qubit actually exists in this stage
    if (std::find(global_qubits.begin(), global_qubits.end(), qubit) !=
        global_qubits.end()) {
      reordered_stage.push_back(qubit);
    }
  }

  return reordered_stage;
}

// Main optimization function
std::vector<std::vector<int>>
optimize_stages(const std::vector<std::vector<int>>& stages,
                int num_local_qubits, int num_frozen_qubits) {
  if (stages.empty())
    return {};

  int n = stages.size();
  std::vector<std::vector<int>> result = stages;  // Start with original stages

  // DP state: dp_cost[i] = minimum hyper-stages needed from stage i to end
  // dp_frozen[i] = optimal frozen qubits for hyper-stage starting at i
  // dp_end[i] = end position of hyper-stage starting at i
  std::vector<int> dp_cost(n + 1, 0);  // dp_cost[n] = 0 (base case)
  std::vector<std::vector<int>> dp_frozen(n);
  std::vector<int> dp_end(n);

  // Fill DP table backwards
  for (int i = n - 1; i >= 0; i--) {
    int best_cost = INT_MAX;
    std::vector<int> best_frozen;
    int best_end = i;

    // Try all possible hyper-stage endings starting from position i
    for (int j = i; j < n; j++) {
      // Check if we can form a valid hyper-stage from i to j
      std::vector<int> candidate_frozen =
          find_optimal_frozen_qubits_with_transition(
              stages, i, j, j + 1, num_local_qubits, num_frozen_qubits);

      // If candidate_frozen is empty, this range is invalid
      if (candidate_frozen.empty() ||
          candidate_frozen.size() != num_frozen_qubits) {
        continue;  // Skip this invalid range
      }

      // Double-check validity: verify that all stages from i to j can provide
      // these frozen qubits
      bool valid = true;
      for (int k = i; k <= j && valid; k++) {
        auto stage_global = get_global_qubits(stages[k], num_local_qubits);
        auto common = find_common_elements(candidate_frozen, stage_global);
        if (common.size() < candidate_frozen.size()) {
          valid = false;
        }
      }

      if (valid) {
        int cost = 1 + dp_cost[j + 1];

        // Add transition bonus (prefer solutions with better transitions)
        int transition_score = calculate_transition_score(
            candidate_frozen, stages, j + 1, num_local_qubits);
        // Use transition score as tie-breaker (smaller cost is better, larger
        // transition score is better)

        if (cost < best_cost ||
            (cost == best_cost &&
             transition_score > calculate_transition_score(best_frozen, stages,
                                                           best_end + 1,
                                                           num_local_qubits))) {
          best_cost = cost;
          best_frozen = candidate_frozen;
          best_end = j;
        }
      }
    }

    dp_cost[i] = best_cost;
    dp_frozen[i] = best_frozen;
    dp_end[i] = best_end;
  }

  // Reconstruct the solution
  int current_pos = 0;
  while (current_pos < n) {
    std::vector<int> current_frozen = dp_frozen[current_pos];
    int hyper_stage_end = dp_end[current_pos];

    // Reorder all stages in this hyper-stage
    for (int k = current_pos; k <= hyper_stage_end; k++) {
      result[k] = reorder_stage_with_frozen_at_end(stages[k], current_frozen,
                                                   num_local_qubits);
    }

    current_pos = hyper_stage_end + 1;
  }

  return result;
}
