import stim

# Define a circuit
circuit = stim.Circuit()


circuit = stim.Circuit("""
    H 0
    CX 0 1
    X_ERROR(0.2) 0 1
    M 0 1
    DETECTOR rec[-1] rec[-2]
""")

print("--- Circuit ---")
print(circuit)

# Simulate the circuit
sampler = circuit.compile_sampler()
result = sampler.sample(shots=10)

# Print the result
print("--- Results ---")
print(result)