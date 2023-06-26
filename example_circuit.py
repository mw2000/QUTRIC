import stim

def generate_2x2_toric_code():
    circuit = stim.Circuit()

    # Qubits are laid out in a 2D grid like this:
    # 0---1
    # |   |
    # 2---3

    # Vertex (star) stabilizers
    circuit.append_operation("MPP", [
        stim.target_x(0),
        stim.target_combiner(),
        stim.target_x(1),
        stim.target_combiner(),
        stim.target_x(2),
    ])
    circuit.append_operation("MPP", [
        stim.target_x(1),
        stim.target_combiner(),
        stim.target_x(0),
        stim.target_combiner(),
        stim.target_x(3),
    ])
    circuit.append_operation("MPP", [
        stim.target_x(2),
        stim.target_combiner(),
        stim.target_x(0),
        stim.target_combiner(),
        stim.target_x(3),
    ])
    circuit.append_operation("MPP", [
        stim.target_x(3),
        stim.target_combiner(),
        stim.target_x(1),
        stim.target_combiner(),
        stim.target_x(2),
    ])

    # Plaquette (face) stabilizers
    circuit.append_operation("MPP", [
        stim.target_z(0),
        stim.target_combiner(),
        stim.target_z(1),
        stim.target_combiner(),
        stim.target_z(2),
    ])
    circuit.append_operation("MPP", [
        stim.target_z(1),
        stim.target_combiner(),
        stim.target_z(0),
        stim.target_combiner(),
        stim.target_z(3),
    ])
    circuit.append_operation("MPP", [
        stim.target_z(2),
        stim.target_combiner(),
        stim.target_z(0),
        stim.target_combiner(),
        stim.target_z(3),
    ])
    circuit.append_operation("MPP", [
        stim.target_z(3),
        stim.target_combiner(),
        stim.target_z(1),
        stim.target_combiner(),
        stim.target_z(2),
    ])

    return circuit

circuit = generate_2x2_toric_code()

print("--- Circuit ---")
print(circuit)

# Simulate the circuit
sampler = circuit.compile_sampler()
result = sampler.sample(shots=10)

# Print the result
print("--- Results ---")
print(result)