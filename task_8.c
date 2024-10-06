#include <stdio.h>
#include <string.h>

#include <lis_config.h>
#include <lis.h>

LIS_INT main(int argc, char *argv[]) {
	LIS_MATRIX A0, A;
	LIS_VECTOR x, b, u;
	LIS_SOLVER solver;
	int nprocs;
	int my_rank;
	LIS_INT nsol;
	LIS_INT rhs;
	LIS_INT err, iter, iter_double, iter_quad;
	double time, itime, ptime, p_c_time, p_i_time;
	LIS_REAL resid;

	LIS_DEBUG_FUNC_IN;

	lis_initialize(&argc, &argv);

	const LIS_Comm comm = LIS_COMM_WORLD;

#ifdef USE_MPI
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &my_rank);
#else
	nprocs = 1;
	my_rank = 0;
#endif

	if (argc < 5) {
		lis_printf(comm, "Usage: %s matrix_filename rhs_setting solution_filename rhistory_filename [options]\n",
				   argv[0]);
		CHKERR(1);
	}

	const LIS_INT len = (LIS_INT)strlen(argv[2]);
	if (len == 1) {
		if (argv[2][0] == '0' || argv[2][0] == '1' || argv[2][0] == '2') {
			rhs = atoi(argv[2]);
		} else {
			rhs = -1;
		}
	} else {
		rhs = -1;
	}

	char *const matrix_filename = argv[1];
	char *const solution_filename = argv[3];
	char *const history_filename = argv[4];

	lis_printf(comm, "\n");
	lis_printf(comm, "number of processes = %d\n", nprocs);

#ifdef _OPENMP
	lis_printf(comm, "max number of threads = %d\n", omp_get_num_procs());
	lis_printf(comm, "number of threads = %d\n", omp_get_max_threads());
#endif

	/* read matrix and vectors from file */
	err = lis_matrix_create(comm, &A);
	CHKERR(err);
	err = lis_vector_create(comm, &b);
	CHKERR(err);
	err = lis_vector_create(comm, &x);
	CHKERR(err);
	err = lis_input(A, b, x, matrix_filename);
	/*
	lis_output_matrix(A,LIS_FMT_MM,"matrix.out");
	lis_output_vector(b,LIS_FMT_MM,"vector.out");
	*/
	CHKERR(err);

	err = lis_matrix_duplicate(A, &A0);
	CHKERR(err);
	lis_matrix_set_type(A0, LIS_MATRIX_CSR);
	err = lis_matrix_convert(A, A0);
	CHKERR(err);
	lis_matrix_destroy(A);
	A = A0;

	err = lis_vector_duplicate(A, &u);
	CHKERR(err);
	if (lis_vector_is_null(b)) {
		lis_vector_destroy(b);
		lis_vector_duplicate(A, &b);
		CHKERR(err);
		if (rhs == 0) {
			CHKERR(1);
		} else if (rhs == 1) {
			err = lis_vector_set_all(1.0, b);
		} else {
			err = lis_vector_set_all(1.0, u);
			lis_matvec(A, u, b);
		}
	}
	if (rhs == -1) {
		lis_vector_destroy(b);
		err = lis_vector_create(comm, &b);
		CHKERR(err);
		lis_input_vector(b, argv[2]);
	}

	if (lis_vector_is_null(x)) {
		lis_vector_destroy(x);
		err = lis_vector_duplicate(A, &x);
		CHKERR(err);
	}

	err = lis_solver_create(&solver);
	CHKERR(err);
	lis_solver_set_option("-i cg -p jacobi -tol 1e-9 -print mem", solver);
	err = lis_solver_set_optionC(solver);
	CHKERR(err);

	err = lis_solve(A, b, x, solver);

	CHKERR(err);
	lis_solver_get_iterex(solver, &iter, &iter_double, &iter_quad);
	lis_solver_get_timeex(solver, &time, &itime, &ptime, &p_c_time, &p_i_time);
	lis_solver_get_residualnorm(solver, &resid);
	lis_solver_get_solver(solver, &nsol);

	char solvername[128];
	lis_solver_get_solvername(nsol, solvername);

	/* write results */
	lis_printf(comm, "%s: number of iterations = %D\n", solvername, iter);
#ifndef _LONG__DOUBLE
	lis_printf(comm, "%s:   double             = %D\n", solvername, iter_double);
	lis_printf(comm, "%s:   quad               = %D\n", solvername, iter_quad);
#endif
	lis_printf(comm, "%s: elapsed time         = %e sec.\n", solvername, time);
	lis_printf(comm, "%s:   preconditioner     = %e sec.\n", solvername, ptime);
	lis_printf(comm, "%s:     matrix creation  = %e sec.\n", solvername, p_c_time);
	lis_printf(comm, "%s:   linear solver      = %e sec.\n", solvername, itime);
	lis_printf(comm, "%s: relative residual    = %e\n\n", solvername, resid);

	/* write solution */
	lis_output_vector(x, LIS_FMT_MM, solution_filename);

	/* write residual history */
	lis_solver_output_rhistory(solver, history_filename);

	lis_solver_destroy(solver);
	lis_vector_destroy(x);
	lis_vector_destroy(u);
	lis_vector_destroy(b);
	lis_matrix_destroy(A);

	lis_finalize();

	LIS_DEBUG_FUNC_OUT;

	return 0;
}