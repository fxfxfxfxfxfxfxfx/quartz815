import sys
from qiskit.circuit.random.utils import random_circuit
from qiskit import transpile


def main(num_qubit: int):
    circ = random_circuit(
        num_qubit, 5 * num_qubit, max_operands=2, measure=False, reset=False
    )
    circ = transpile(
        circ,
        basis_gates=["u", "cz", "x", "rz"],
    )
    circ.qasm(filename=f"../../circuit/qiskit-random/rqc_{num_qubit}.qasm")


if __name__ == "__main__":
    if len(sys.argv) == 2:
        num_qubit = int(sys.argv[1])
        main(num_qubit)
    else:
        for n in range(30, 50, 2):
            main(n)
