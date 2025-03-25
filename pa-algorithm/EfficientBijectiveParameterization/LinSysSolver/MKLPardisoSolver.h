#pragma once
#include "LinSysSolver.hpp"
#include "mkl_pardiso.h"
#include "mkl_types.h"
using namespace std;

class MKLPardisoSolver:public  LinSysSolver
{
public:
	MKLPardisoSolver();
	~MKLPardisoSolver();

	void pardiso_init();
	bool factorize();
	void pardiso_solver();
	void free_numerical_factorization_memory();

};