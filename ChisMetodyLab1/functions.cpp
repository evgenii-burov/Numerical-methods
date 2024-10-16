#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include "functions.h"

int CreateMatrixFiles()
{
	int n = -1;
	int p = -1;
	int matrix_type = -1;
	std::ifstream input;
	input.open("dimension,p,matrix_type.txt");
	if (!input.is_open())
	{
		std::cout << "Unable to read the matrix parameters file, exiting..";
		return -1;
	}
	input >> n >> p >> matrix_type;
	input.close();
	if (n < 1 || p < 0 || matrix_type == -1)
	{
		std::cout << "Invalid matrix parameters, exiting..";
		return -1;
	}
	std::ofstream output;
	output.open("matrix_al.txt");
	switch (matrix_type)
	{
	default:
		std::cout << "Invalid matrix type, exiting..";
		return -1;
	case 1:
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < p+1; j++)
				output << (p  - i > j ? 0 : (j + i + 1)) << '\t'; //(j + i + 1))
			//10*double(rand())/RAND_MAX)
			//(rand()/RAND_MAX + 3)*4
			output << '\n';
		}
		break;
	}
	output.close();

	output.open("vector_b.txt");
	for (int i = 0; i < n; i++)
		output << i*1.5 + 1 << '\t';
	output.close();
}

void ManageInput(precision**& al, precision*& b, int& n, int& p)
{
	std::ifstream input;
	input.open("dimension,p,matrix_type.txt");
	input >> n >> p;
	input.close();
	input.open("vector_b.txt");
	b = new precision[n];
	for (int i = 0; i < n; i++)
		input >> b[i];
	input.close();
	input.open("matrix_al.txt");
	al = new precision* [n];
	for (int i = 0; i < n; i++)
	{
		al[i] = new precision[p+1];
		for (int j = 0; j < p+1; j++)
			input >> al[i][j];
	}
	input.close();
}

void CalculateDecomposition(precision**& al, const int& n, const int& p)
{
	precision l_i_k = 0;
	precision d_k_k = 0;
	precision l_j_k = 0;
	precision sum_over_k = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0, lij=p-i; j < i; lij++, j++)
		{
			sum_over_k = 0;
			for (int k = j-1, lik = p - i+j-1, ljk = p - 1; k>=0 && lik >=0; k--, lik--, ljk--)
			{
				sum_over_k += al[i][lik] * al[k][p] * al[j][ljk];
			}
			al[i][lij] = (al[i][lij] - sum_over_k) / al[j][p];
		}
		sum_over_k = 0;
		for (int k = i - 1, lik = p - 1; k >= 0 && lik >= 0; k--, lik--)
		{
			sum_over_k += al[i][lik] * al[k][p] * al[i][lik];
		}
		al[i][p] = (al[i][p] - sum_over_k);
	}
}

void PrintVariables(precision** al, precision* b, const int& n,const int& p)
{
	std::cout << std::setprecision(4) << "n: " << n << "\tp: " << p;
	std::cout << "\nb: ";
	for (int i = 0; i < n; i++)
	{
		std::cout << b[i] << '\t';
	}
	std::cout << "\nal:\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j <p+1; j++)
			std::cout << al[i][j] << '\t';
		std::cout << "\n";
	}
	std::cout << "\n\n";
}

void SolveForX(precision** al, precision* &b, const int& n, const int& p)
{
	precision sum_over_k = 0;
	for (int i = 0; i < n; i++)
	{
		sum_over_k = 0;
		for (int lk = 0, k = i-p; lk < p; k++, lk++)
		{
			sum_over_k += al[i][lk] * b[k];
		}
		b[i] -= sum_over_k;
	}
	for (int i = n - 1; i >= 0; i--)
	{
		sum_over_k = 0;
		for (int k = i + 1, li = p - 1; k < n && k <=p+i; k++, li--)//int k = i + 1; k < n; k++
		{
			sum_over_k += al[k][li] * b[k];
		}
		b[i] = b[i] / al[i][p] - sum_over_k;
	}
	std::cout << "x: ";
	for (int i = 0; i < n; i++)
		std::cout << b[i] << ' ';
	std::ofstream output;
	output.open("vector_x.txt");
	output << std::setprecision(4);
	for (int i = 0; i < n; i++)
		output << b[i] << ' ';
	output.close();
}

void PrintDenseMatrix(precision** al, const int& n, const int& p)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			std::cout << (i-j<=p ? al[i][p - i + j] : 0)<<'\t';
		}
		std::cout << "\n\n";
	}
}
