#pragma once

using precision = double;

int CreateMatrixFiles();

void ManageInput(precision**& al, precision*& b, int& n, int& p);

void CalculateDecomposition(precision**& al, const int& n, const int& p);

void PrintVariables(precision** al, precision* b, const int& n, const int& p);

void SolveForX(precision** al, precision* &b, const int& n, const int& p);

void PrintDenseMatrix(precision** al, const int& n, const int& p);
