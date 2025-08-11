from qiskit.circuit.random.utils import random_circuit
from qiskit import transpile
import argparse


def main(num_qubit: int, depth_factor: int):
    circ = random_circuit(
        num_qubit, depth_factor * num_qubit, max_operands=2, measure=False, reset=False
    )
    circ = transpile(
        circ,
        basis_gates=["u", "cz", "x", "rz"],
    )
    circ.qasm(
        filename=f"../../circuit/qiskit-random/rqc_depth_{depth_factor}_{num_qubit}.qasm"
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Generate random quantum circuits.")
    parser.add_argument("num_qubit", type=int, help="Number of qubits")
    parser.add_argument(
        "--depth-factor",
        type=int,
        default=5,
        help="Depth factor for the circuit (default: 5)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(args.num_qubit, args.depth_factor)
