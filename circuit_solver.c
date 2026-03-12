#include "circuit_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// helper function to solve the linear system Ax = b using Gaussian elimination
void solve(int n, double *A, double *b, double *x);

// Function to read a circuit from user input
circuit input_circuit() {
    printf(" =================================================\n");
    printf("                   CIRCUIT INPUT                  \n");
    printf(" =================================================\n");
    circuit c;
    printf(" Enter the number of nodes: ");
    scanf("%d", &c.num_nodes);
    assert(c.num_nodes > 0 && c.num_nodes <= MAX_NODES);
    printf("\n");
    printf(" Enter the number of resistors: ");
    scanf("%d", &c.num_resistors);
    assert(c.num_resistors >= 0 && c.num_resistors <= MAX_RES);
    printf("\n");
    printf(" Enter the number of voltage sources: ");
    scanf("%d", &c.num_voltage_sources);
    assert(c.num_voltage_sources >= 0 && c.num_voltage_sources <= MAX_VSRC);
    printf("\n");
    printf(" Enter the number of current sources: ");
    scanf("%d", &c.num_current_sources);
    assert(c.num_current_sources >= 0 && c.num_current_sources <= MAX_ISRC);
    printf("\n");
    printf(" Enter the resistors:\n");
    printf("\n");
    for (int i = 0; i < c.num_resistors; i++) {
        printf(" Resistor %d:\n", i + 1);
        printf("  Positive node: ");
        scanf("%d", &c.resistors[i].node_pos);
        printf("  Negative node: ");
        scanf("%d", &c.resistors[i].node_neg);
        printf("  Resistance: ");
        scanf("%lf", &c.resistors[i].resistance);
        printf("\n");
    }
    printf(" Enter the voltage sources:\n");
    printf("\n");
    for (int i = 0; i < c.num_voltage_sources; i++) {
        printf(" Voltage source %d:\n", i + 1);
        printf("  Positive node: ");
        scanf("%d", &c.voltage_sources[i].node_pos);
        printf("  Negative node: ");
        scanf("%d", &c.voltage_sources[i].node_neg);
        printf("  Voltage: ");
        scanf("%lf", &c.voltage_sources[i].voltage);
        printf("\n");
    }
    printf(" Enter the current sources:\n");
    printf("\n");
    for (int i = 0; i < c.num_current_sources; i++) {
        printf(" Current source %d:\n", i + 1);
        printf("  Positive node: ");
        scanf("%d", &c.current_sources[i].node_pos);
        printf("  Negative node: ");
        scanf("%d", &c.current_sources[i].node_neg);
        printf("  Current: ");
        scanf("%lf", &c.current_sources[i].current);
        printf("\n");
    }
    printf("\n");
    return c;
}

// Function to read a circuit from a file
circuit read_circuit_from_file(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        exit(1);
    }
    circuit c;
    fscanf(file, "%d", &c.num_nodes);
    fscanf(file, "%d", &c.num_resistors);
    fscanf(file, "%d", &c.num_voltage_sources);
    fscanf(file, "%d", &c.num_current_sources);
    for (int i = 0; i < c.num_resistors; i++) {
        fscanf(file, "%d %d %lf", &c.resistors[i].node_pos, &c.resistors[i].node_neg, &c.resistors[i].resistance);
    }
    for (int i = 0; i < c.num_voltage_sources; i++) {
        fscanf(file, "%d %d %lf", &c.voltage_sources[i].node_pos, &c.voltage_sources[i].node_neg, &c.voltage_sources[i].voltage);
    }
    for (int i = 0; i < c.num_current_sources; i++) {
        fscanf(file, "%d %d %lf", &c.current_sources[i].node_pos, &c.current_sources[i].node_neg, &c.current_sources[i].current);
    }
    fclose(file);
    return c;
}

// Function to solve the circuit using Gaussian elimination
void solve_circuit(circuit c, double * voltages, double *currents) {
    int num_voltage_nodes = c.num_nodes - 1;
    int n = num_voltage_nodes + c.num_voltage_sources;
    double * A = malloc(n * n * sizeof(double));
    assert(A != NULL);
    double * b = malloc(n * sizeof(double));
    assert(b != NULL);
    // Initialize A and b to zero
    for (int i = 0; i < n * n; i++) {
        A[i] = 0.0;
    }
    for (int i = 0; i < n; i++) {
        b[i] = 0.0;
    }
    // Fill the conductance submatrix
    for (int i = 0; i < c.num_resistors; i++) {
        int n1 = c.resistors[i].node_pos;
        int n2 = c.resistors[i].node_neg;
        double g = 1.0 / c.resistors[i].resistance;

        int i1 = n1 - 1; // index of node n1 in the voltage nodes (excluding ground)
        int i2 = n2 - 1; // index of node n2 in the voltage nodes (excluding ground)

        if (n1 != 0) {
            A[i1 * n + i1] += g; // A[i1][i1] += g
        }
        if (n2 != 0) {
            A[i2 * n + i2] += g; // A[i2][i2] += g
        }
        if (n1 != 0 && n2 != 0) {
            A[i1 * n + i2] -= g; // A[i1][i2] -= g
            A[i2 * n + i1] -= g; // A[i2][i1] -= g
        }
    }
    // Fill the voltage constraints
    for (int k = 0; k < c.num_voltage_sources; k++) {
        int n1 = c.voltage_sources[k].node_pos;
        int n2 = c.voltage_sources[k].node_neg;
        
        int row = num_voltage_nodes + k; // index of the voltage source constraint row
        if (n1 != 0) {
            int i1 = n1 - 1; // index of node n1 in the voltage nodes (excluding ground)
            A[i1 * n + row] += 1.0; // A[i1][row] += 1
            A[row * n + i1] += 1.0; // A[row][i1] += 1
        }
        if (n2 != 0) {
            int i2 = n2 - 1; // index of node n2 in the voltage nodes (excluding ground)
            A[i2 * n + row] -= 1.0; // A[i2][row] -= 1
            A[row * n + i2] -= 1.0; // A[row][i2] -= 1
        }
    }
    // Fill the voltage source constraints
    for (int k = 0; k < c.num_voltage_sources; k++) {
        b[num_voltage_nodes + k] = c.voltage_sources[k].voltage;
    }
    // Solve the system Ax = b
    double *x = malloc(n * sizeof(double));
    assert(x != NULL);
    solve(n, A, b, x);

    // Copy the node voltages to the output array
    // Ignore the voltage of the ground node (node 0) which is always 0
    for (int i = 0; i < num_voltage_nodes; i++) {
        voltages[i] = x[i];
    }
    // Calculate the currents through the resistors using Ohm's law and copy to the output array
    for (int i = 0; i < c.num_resistors; i++) {
        int n1 = c.resistors[i].node_pos;
        int n2 = c.resistors[i].node_neg;
        double v1 = (n1 == 0) ? 0.0 : voltages[n1 - 1];
        double v2 = (n2 == 0) ? 0.0 : voltages[n2 - 1];
        currents[i] = (v1 - v2) / c.resistors[i].resistance;
    }
    // Copy the currents through the voltage sources to the output array
    for (int k = 0; k < c.num_voltage_sources; k++) {
        currents[c.num_resistors + k] = x[num_voltage_nodes + k];
    }
    free(A);
    free(b);
    free(x);
}

// Function to print the circuit
void print_circuit(circuit c) {
    printf(" =================================================\n");
    printf("                  CIRCUIT SUMMARY                 \n");
    printf(" =================================================\n");
    printf(" Number of nodes          : %d\n", c.num_nodes);
    printf(" Number of resistors      : %d\n", c.num_resistors);
    printf(" Number of voltage sources: %d\n", c.num_voltage_sources);
    printf(" Number of current sources: %d\n", c.num_current_sources);
    printf("\n");

    if (c.num_resistors > 0) {
        printf(" Resistors:\n");
        printf(" +------------+----------+----------+------------+\n");
        printf(" | Resistor   | Node+    | Node-    | Resistance |\n");
        printf(" +------------+----------+----------+------------+\n");
        for (int i = 0; i < c.num_resistors; i++) {
            printf(" | %-10d | %-8d | %-8d | %-10.2f |\n",
                i + 1,
                c.resistors[i].node_pos,
                c.resistors[i].node_neg,
                c.resistors[i].resistance
            );
        }
        printf(" +------------+----------+----------+------------+\n\n");
    }

    if (c.num_voltage_sources > 0) {
        printf("Voltage Sources:\n");
        printf(" +------------+----------+----------+------------+\n");
        printf(" | Source     | Node+    | Node-    | Voltage(V) |\n");
        printf(" +------------+----------+----------+------------+\n");
        for (int i = 0; i < c.num_voltage_sources; i++) {
            printf(" | %-10d | %-8d | %-8d | %-10.2f |\n",
                i + 1,
                c.voltage_sources[i].node_pos,
                c.voltage_sources[i].node_neg,
                c.voltage_sources[i].voltage
            );
        }
        printf(" +------------+----------+----------+------------+\n\n");
    }

    if (c.num_current_sources > 0) {
        printf("Current Sources:\n");
        printf(" +------------+----------+----------+------------+\n");
        printf(" | Source     | Node+    | Node-    | Current(A) |\n");
        printf(" +------------+----------+----------+------------+\n");
        for (int i = 0; i < c.num_current_sources; i++) {
            printf(" | %-10d | %-8d | %-8d | %-10.2f |\n",
                i + 1,
                c.current_sources[i].node_pos,
                c.current_sources[i].node_neg,
                c.current_sources[i].current
            );
        }
        printf( "+------------+----------+----------+------------+\n\n");
    }

    printf(" =================================================\n");
}

// Function to print the solution
void print_solution(circuit c, double *voltages, double *currents) {
    printf(" =================================================\n");
    printf("                 CIRCUIT SOLUTION                 \n");
    printf(" =================================================\n");

    // Node voltages
    if (c.num_nodes > 1) {
        printf(" Node Voltages (relative to ground):\n");
        printf(" +-----------------------+-----------------------+\n");
        printf(" | Node                  | Voltage(V)            |\n");
        printf(" +-----------------------+-----------------------+\n");
        for (int i = 1; i < c.num_nodes; i++) { // node 0 is ground
            printf(" | %-21d | %-21.2f |\n", i, voltages[i - 1]);
        }
        printf(" +-----------------------+-----------------------+\n\n");
    }

    // Currents through resistors
    if (c.num_resistors > 0) {
        printf("Currents through Resistors:\n");
        printf(" +-----------------------+-----------------------+\n");
        printf(" | Resistor              | Current(A)            |\n");
        printf(" +-----------------------+-----------------------+\n");
        for (int i = 0; i < c.num_resistors; i++) {
            printf(" | %-21d | %-21.2f |\n", i + 1, currents[i]);
        } 
        printf(" +-----------------------+-----------------------+\n\n");
    }

    // Currents through voltage sources
    if (c.num_voltage_sources > 0) {
        printf(" Currents through Voltage Sources:\n");
        printf(" +-----------------------+-----------------------+\n");
        printf(" | Source                | Current(A)            |\n");
        printf(" +-----------------------+-----------------------+\n");
        for (int k = 0; k < c.num_voltage_sources; k++) {
            printf(" | %-21d | %-21.2f |\n", k + 1, currents[c.num_resistors + k]);
        }
         printf(" +-----------------------+-----------------------+\n\n");
    }

    printf(" =================================================\n");
}


// helper function to solve the linear system Ax = b using Gaussian elimination
void solve(int n, double *A, double *b, double *x)
{
    for (int k = 0; k < n; k++) {

        // Pivot
        double pivot = A[k*n + k];

        for (int j = k; j < n; j++)
            A[k*n + j] /= pivot;

        b[k] /= pivot;

        // Eliminate rows below
        for (int i = k + 1; i < n; i++) {

            double factor = A[i*n + k];

            for (int j = k; j < n; j++)
                A[i*n + j] -= factor * A[k*n + j];

            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; i--) {

        x[i] = b[i];

        for (int j = i + 1; j < n; j++)
            x[i] -= A[i*n + j] * x[j];
    }
}