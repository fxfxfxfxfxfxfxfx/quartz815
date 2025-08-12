import sys
from qiskit import QuantumCircuit, transpile

QASM_PREFIX = "../../circuit/veriq-bench/combinational/grover/grover_"
QASM_PREFIX_OUT = "../../circuit/misc/grover_"


if __name__ == "__main__":
    num_qubits = sys.argv[1]
    qasm_file = QASM_PREFIX + num_qubits + ".qasm"

    qc = QuantumCircuit.from_qasm_file(qasm_file)
    qc = transpile(qc, basis_gates=["cx", "u3"])

    output_file = QASM_PREFIX_OUT + num_qubits + "_transpiled.qasm"
    qc.qasm(filename=output_file)
