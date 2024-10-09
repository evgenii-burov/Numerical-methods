#pragma once

using precision = double;

// ќткрывает файл dimension,p,matrix_type, считывает из него размерность n, 
// параметр p и параметр, определ€ющим тип стро€щейс€ матрицы 
// и записывает файлы matrix_al и vector_b
int CreateMatrixFiles();


void ManageInput(precision**& al, precision*& b, int& n, int& p);

void CalculateDecomposition(precision**& al, const int& n, const int& p);

void PrintVariables(precision** al, precision* b, const int& n, const int& p);

void SolveForX(precision** al, precision* &b, const int& n, const int& p);

void PrintFullMatrix(precision** al, const int& n, const int& p);

void CalculateError(precision**& al, precision* b, const int& n, const int& p);