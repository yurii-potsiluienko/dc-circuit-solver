#include "circuit_solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


int main() {
    // Read the circuit from user input
    circuit c = input_circuit();
    // Alternatively, read the circuit from a file
    // circuit c = read_circuit_from_file("test_circuit.txt");
    print_circuit(c);
    // Solve the circuit
    double * voltages = malloc((c.num_nodes - 1) * sizeof(double));
    assert(voltages != NULL);
    double * currents = malloc((c.num_resistors + c.num_voltage_sources) * sizeof(double));
    assert(currents != NULL);
    
    solve_circuit(c, voltages, currents);
    //
    print_solution(c, voltages, currents);
    // Free the memory
    free(voltages);
    free(currents);
    return 0;
}