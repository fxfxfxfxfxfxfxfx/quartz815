#!/usr/bin/env python
# coding: utf-8

import argparse
from qiskit_aer import AerSimulator
from qiskit.circuit.library import EfficientSU2
import numpy as np
from scipy.optimize import minimize
from qiskit import transpile, qasm3
from qiskit.quantum_info import Statevector, state_fidelity
from qiskit.circuit import QuantumCircuit, ParameterVector
import matplotlib.pyplot as plt
import sys


def conv_circuit(params):
    target = QuantumCircuit(2)
    target.rz(-np.pi / 2, 1)
    target.cx(1, 0)
    target.rz(params[0], 0)
    target.ry(params[1], 1)
    target.cx(0, 1)
    target.ry(params[2], 1)
    target.cx(1, 0)
    target.rz(np.pi / 2, 0)
    return target


def conv_layer(num_qubits, param_prefix):
    qc = QuantumCircuit(num_qubits, name="Convolutional Layer")
    qubits = list(range(num_qubits))
    param_index = 0
    params = ParameterVector(param_prefix, length=num_qubits * 3)
    for q1, q2 in zip(qubits[0::2], qubits[1::2]):
        qc = qc.compose(conv_circuit(params[param_index : (param_index + 3)]), [q1, q2])
        qc.barrier()
        param_index += 3
    for q1, q2 in zip(qubits[1::2], qubits[2::2] + [0]):
        qc = qc.compose(conv_circuit(params[param_index : (param_index + 3)]), [q1, q2])
        qc.barrier()
        param_index += 3

    qc_inst = qc.to_instruction()

    qc = QuantumCircuit(num_qubits)
    qc.append(qc_inst, qubits)
    return qc


def pool_circuit(params):
    target = QuantumCircuit(2)
    target.rz(-np.pi / 2, 1)
    target.cx(1, 0)
    target.rz(params[0], 0)
    target.ry(params[1], 1)
    target.cx(0, 1)
    target.ry(params[2], 1)

    return target


def pool_layer(sources, sinks, param_prefix):
    num_qubits = max(max(sources, default=0), max(sinks, default=0)) + 1
    qc = QuantumCircuit(num_qubits, name="Pooling Layer")
    param_index = 0
    params = ParameterVector(param_prefix, length=len(sources) * 3)
    for source, sink in zip(sources, sinks):
        qc = qc.compose(
            pool_circuit(params[param_index : (param_index + 3)]), [source, sink]
        )
        qc.barrier()
        param_index += 3

    qc_inst = qc.to_instruction()

    qc = QuantumCircuit(num_qubits)
    qc.append(qc_inst, range(num_qubits))
    return qc


def build_ansatz(num_qubits):
    ansatz = QuantumCircuit(num_qubits, name="Ansatz")
    # First Convolutional Layer
    ansatz.compose(conv_layer(num_qubits, "c1"), list(range(num_qubits)), inplace=True)
    # First Pooling Layer
    half = num_qubits // 2
    sources1 = list(range(half))
    sinks1 = list(range(half, num_qubits))
    if len(sources1) > 0 and len(sinks1) > 0:
        ansatz.compose(
            pool_layer(sources1, sinks1, "p1"), list(range(num_qubits)), inplace=True
        )
    # Second Convolutional Layer
    if half > 0:
        ansatz.compose(
            conv_layer(half, "c2"), list(range(half, num_qubits)), inplace=True
        )
        # Second Pooling Layer
        if half // 2 > 0:
            sources2 = list(range(half // 2))
            sinks2 = list(range(half // 2, half))
            if len(sources2) > 0 and len(sinks2) > 0:
                ansatz.compose(
                    pool_layer(sources2, sinks2, "p2"),
                    list(range(half, num_qubits)),
                    inplace=True,
                )
    return ansatz


def main():
    parser = argparse.ArgumentParser(
        description="Quantum Convolutional Neural Network (QCNN)"
    )
    parser.add_argument(
        "--num_qubits",
        type=int,
        default=4,
        help="Number of qubits for the QCNN circuit",
    )
    args = parser.parse_args()

    num_qubits = args.num_qubits

    simulator_aer = AerSimulator()

    ansatz = build_ansatz(num_qubits)

    def cost_function(params):
        # Create a new circuit with the ansatz and assign the parameters
        bound_circuit = ansatz.assign_parameters(params)
        bound_circuit.save_statevector()
        # Get the statevector from the circuit
        qc_aer = transpile(bound_circuit, backend=simulator_aer)
        result = simulator_aer.run(qc_aer).result()
        output_state = result.get_statevector()

        # Calculate the probability of measuring 1 at the last qubit
        # For generality, sum over all basis states where the last qubit is 1
        n = num_qubits
        prob = -sum(np.abs(output_state[i]) ** 2 for i in range(2**n) if (i & 1) == 1)
        return prob

    # Number of parameters in the ansatz
    num_params = ansatz.num_parameters

    # Initialize parameters randomly
    initial_params = np.random.uniform(0, 2 * np.pi, num_params)

    # Perform the optimization
    # result = minimize(
    #     cost_function,
    #     initial_params,
    #     method="COBYLA",
    #     options={"maxiter": 1000, "disp": True},
    # )
    #
    # # Extract optimized parameters
    # optimized_params = result.x
    #
    # print("Optimization Success:", result.success)
    # print("Final Infidelity:", result.fun)

    # Bind the optimized parameters
    final_circuit = ansatz.assign_parameters(initial_params)
    # final_circuit = ansatz
    # Save the qasm file of the circuit
    # qasm_str = qasm3.dumps(final_circuit.decompose())
    qasm_str = final_circuit.decompose().qasm(filename=f"qcnn_{num_qubits}.qasm")
    # with open(f"qcnn_{num_qubits}.qasm", "w") as f:
    #     f.write(qasm_str)
    # final_circuit.save_statevector()
    # # Execute the circuit
    # qc_aer = transpile(final_circuit, backend=simulator_aer)
    # result = simulator_aer.run(qc_aer).result()
    # final_state = result.get_statevector()
    #
    # # Calculate fidelity
    # n = num_qubits
    # prob = sum(np.abs(final_state[i]) ** 2 for i in range(2**n) if (i & 1) == 1)
    # print("Final Fidelity:", prob)
    #
    # # Optional: Visualize the final state
    # state = Statevector(final_state)
    # # state.draw('bloch')
    #
    # print("Final State:", np.round(final_state, 2))
    # print(qasm_str)


if __name__ == "__main__":
    main()
