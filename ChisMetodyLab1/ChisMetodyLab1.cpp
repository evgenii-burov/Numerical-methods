﻿#include <iostream>
#include <fstream>
#include <iomanip>
#include "functions.h"
#include <random>

int main()
{
	srand(100);
	if (CreateMatrixFiles() == -1)
		return 0;
	precision** al = NULL;
	precision* b = NULL;
	int n = 0;
	int p = 0;
	ManageInput(al, b, n, p);
	PrintVariables(al, b, n, p);
	PrintFullMatrix(al, n, p);
	CalculateDecomposition(al, n, p);
	PrintVariables(al, b, n, p);
	PrintFullMatrix(al, n, p);
	//DecompositionCheck(al, n);
	SolveForX(al, b, n, p);
	//CalculateError(al, b, n, p);
	for (int i = 0; i < n; i++)
	{
		delete[] al[i];
	}
	delete[] al;
	delete[] b;
}
