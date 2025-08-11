from qiskit import transpile
import argparse
from qiskit.circuit.library import QuantumVolume


def generate_complex_qc(num_qubits: int, depth: int) -> QuantumVolume:
    """
    Generate a complex quantum circuit using the QuantumVolume class from Qiskit.

    Args:
        num_qubits (int): The number of qubits in the quantum circuit.
        depth (int): The depth of the quantum circuit.

    Returns:
        QuantumVolume: A quantum circuit with the specified number of qubits and depth.
    """
    return QuantumVolume(num_qubits=num_qubits, depth=depth, seed=42)


def save_qc_to_qasm(qc: QuantumVolume, filename: str) -> None:
    """
    Save the quantum circuit to a QASM file.

    Args:
        qc (QuantumVolume): The quantum circuit to save.
        filename (str): The name of the file to save the circuit to.
    """
    qc.qasm(filename=filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate and save a quantum circuit to QASM."
    )
    parser.add_argument(
        "--num-qubits", type=int, default=5, help="Number of qubits in the circuit"
    )
    parser.add_argument(
        "--depth-factor", type=int, default=1, help="Depth of the circuit"
    )
    args = parser.parse_args()

    num_qubits = args.num_qubits
    depth = args.depth_factor * num_qubits
    qc = generate_complex_qc(num_qubits, depth)
    qc = transpile(qc, basis_gates=["u", "cz", "x", "rz"])
    save_qc_to_qasm(
        qc,
        f"../../circuit/qiskit-qv/qv_depth_{args.depth_factor}_{num_qubits}.qasm",
    )
    print(
        f"Quantum circuit with {num_qubits} qubits and depth {depth} saved to QASM file."
    )
