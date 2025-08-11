"""
Generate Hardware-Efficient Ansatz quantum circuit

http://www.nature.com/articles/nature23879

The circuit structure is repeated blocks of
1. Single qubit rotations
2. Entangling gates (CNOTs)

"""

from qiskit import QuantumCircuit, transpile
import random


def generate_hea(num_qubits, num_blocks):
    """
    Generate a Hardware-Efficient Ansatz circuit.

    Args:
        num_qubits (int): Number of qubits in the circuit.
        num_blocks (int): Number of blocks to repeat.

    Returns:
        QuantumCircuit: The generated quantum circuit.
    """
    qc = QuantumCircuit(num_qubits)

    for _ in range(num_blocks):
        # Single qubit rotations
        for qubit in range(num_qubits):
            angle_rx = random.uniform(
                0, 2 * 3.141592653589793
            )  # Random angle between 0 and 2pi
            angle_ry = random.uniform(0, 2 * 3.141592653589793)
            qc.rx(angle_rx, qubit)
            qc.ry(angle_ry, qubit)

        # Entangling gates (CNOTs)
        # # randomized connections
        # qubits = list(range(num_qubits))
        # random.shuffle(qubits)
        for i in range(0, num_qubits - 1, 2):
            # qc.cx(qubits[i], qubits[i + 1])
            qc.cx(i, i + 1)
            # qc.h(qubits[i])  # Apply X gate to the control qubit

    # return transpile(qc, basis_gates=["u", "cx", "rz"])
    return qc


def save_qc_to_qasm(qc, filename):
    """
    Save the quantum circuit to a QASM file.

    Args:
        qc (QuantumCircuit): The quantum circuit to save.
        filename (str): The name of the file to save the circuit to.
    """
    qc.qasm(filename=filename)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate and save a Hardware-Efficient Ansatz quantum circuit."
    )
    parser.add_argument(
        "--num-qubits", type=int, default=5, help="Number of qubits in the circuit"
    )
    parser.add_argument(
        "--num-blocks", type=int, default=3, help="Number of blocks to repeat"
    )
    args = parser.parse_args()

    num_qubits = args.num_qubits
    num_blocks = args.num_blocks
    qc = generate_hea(num_qubits, num_blocks)
    save_qc_to_qasm(
        qc, f"../../circuit/qiskit-hea/hea_blocks_{args.num_blocks}_{num_qubits}.qasm"
    )
    print(
        f"Hardware-Efficient Ansatz circuit with {num_qubits} qubits and {num_blocks} blocks saved to QASM file."
    )
