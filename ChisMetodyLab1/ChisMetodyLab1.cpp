#include <iostream>
#include <fstream>
#include <iomanip>
#include "functions.h"

int main()
{
	if (CreateMatrixFiles() == -1)
		return 0;
	precision** al = NULL;
	precision* b = NULL;
	precision* x = NULL;
	int n = 0;
	int p = 0;
	ManageInput(al, b, n, p);
	PrintVariables(al, b, n);
	CalculateDecomposition(al, b, n, p);
	PrintVariables(al, b, n);
	//DecompositionCheck(al, n);
	SolveForX(al, b, x, n);
	for (int i = 0; i < n; i++)
		delete[] al[i];
	delete[] al;
	delete[] b;
	delete[] x;
}
