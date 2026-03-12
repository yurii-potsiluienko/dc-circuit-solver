#ifndef CIRCUIT_SOLVER
#define CIRCUIT_SOLVER

#define MAX_NODES 50
#define MAX_RES 100
#define MAX_VSRC 20
#define MAX_ISRC 20


typedef struct {
    int node_pos;
    int node_neg;
    double resistance;
} resistor;

typedef struct {
    int node_pos;
    int node_neg;
    double voltage;
} voltage_source;

typedef struct {
    int node_pos;
    int node_neg;
    double current;
} current_source;

typedef struct {
    int num_nodes;
    int num_resistors;
    int num_voltage_sources;
    int num_current_sources;
    resistor resistors[MAX_RES];
    voltage_source voltage_sources[MAX_VSRC];
    current_source current_sources[MAX_ISRC];
} circuit;

// Function to read a circuit from user input
circuit input_circuit();

// Function to read a circuit from a file
circuit read_circuit_from_file(const char *filename);

// Function to solve the circuit using Gaussian elimination
void solve_circuit(circuit c, double * voltages, double *currents);

// Function to print the circuit
void print_circuit(circuit c);

// Function to print the solution
void print_solution(circuit c, double *voltages, double *currents);

#endif