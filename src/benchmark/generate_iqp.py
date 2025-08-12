"""
Generate Instantaneous quantum polynomial circuit using qiskit
"""

import numpy as np
from qiskit import transpile
import argparse
from qiskit.circuit.library import IQP


def generate_iqp(n: int):
    """
    Generate an IQP circuit using the IQP class from Qiskit.

    Args:
        n (int): Dimension of the interactions matrix

    Returns:
        IQP: An IQP circuit with the specified number of qubits and depth.
    """
    mat = np.random.rand(n, n)
    sym_mat = (mat + mat.T) / 2  # Make it symmetric
    circ = IQP(sym_mat)
    return transpile(circ, basis_gates=["u", "cz", "x", "rz"])


def save_qc_to_qasm(qc: IQP, filename: str) -> None:
    """
    Save the quantum circuit to a QASM file.

    Args:
        qc (IQP): The quantum circuit to save.
        filename (str): The name of the file to save the circuit to.
    """
    qc.qasm(filename=filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate and save an IQP circuit to QASM."
    )
    parser.add_argument(
        "num_qubits", type=int, default=5, help="Number of qubits in the circuit"
    )
    args = parser.parse_args()

    num_qubits = args.num_qubits
    qc = generate_iqp(num_qubits)
    save_qc_to_qasm(
        qc,
        f"../../circuit/qiskit-iqp/iqp_{num_qubits}.qasm",
    )
    print(f"IQP circuit with {num_qubits} qubits saved to QASM file.")
