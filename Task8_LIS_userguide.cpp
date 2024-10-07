// Task 8: Export the Eigen matrix A2 and vector w in the .mtx format. 
// Using a suitable iterative solver and preconditioner technique available in the LIS library
// compute the approximate solution to the linear system A2 * x = w prescribing a tolerance of 1e−9.

// Esportazione della matrice A2 e del vettore w nel formato .mtx 
Eigen::saveMarket(A2, "A2.mtx");
std::cout << "Matrix A2 saved to A2.mtx" << std::endl;
Eigen::saveMarketVector(w, "w.mtx");
std::cout << "Vector w saved to w.mtx" << std::endl;

// Inizializza LIS
lis_initialize(&argc, &argv);

// Dichiarazioni LIS per matrice, vettori e solver
LIS_MATRIX A;
LIS_VECTOR b, x;
LIS_SOLVER solver;
int iter_count;
double residual;

// Carica la matrice A2 da "A2.mtx" e il vettore w da "w.mtx"
lis_matrix_create(0, &A);
lis_input_matrix(A, "A2.mtx");

lis_vector_create(0, &b);
lis_input_vector(b, "w.mtx");

// Crea il vettore soluzione x (inizialmente vuoto)
lis_vector_duplicate(b, &x);

// Crea il solver e imposta il metodo iterativo e il precondizionatore
lis_solver_create(&solver);
lis_solver_set_option("-i cg -p jacobi", solver);  // Metodo CG con precondizionatore Jacobi

// Imposta la tolleranza a 1e-9
lis_solver_set_option("-tol 1.0e-9", solver);

// Risolvi il sistema lineare A2 * x = w
lis_solve(A, b, x, solver);

// Ottieni il numero di iterazioni e il residuo finale
lis_solver_get_iter(solver, &iter_count);
lis_solver_get_residualnorm(solver, &residual);

// Stampa i risultati
std::cout << "Numero di iterazioni: " << iter_count << std::endl;
std::cout << "Residuo finale: " << residual << std::endl;

// Libera la memoria
lis_solver_destroy(solver);
lis_matrix_destroy(A);
lis_vector_destroy(b);
lis_vector_destroy(x);

// Finalizza LIS
lis_finalize();
