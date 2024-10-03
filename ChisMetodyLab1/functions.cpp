#include <iostream>
#include <fstream>
#include <iomanip>
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
				output << (p - i > j ? 0 : (j + i + 1)) << '\t';
			output << '\n';
		}
		break;
	}
	output.close();

	output.open("vector_b.txt");
	for (int i = 0; i < n; i++)
		output << i + 1 << '\t';
	output.close();
}

void ManageInput(precision**& al, precision*& b, int& n, int& p)
{
	std::ifstream input;
	input.open("dimension,p,matrix_type.txt");
	input >> n >> p;
	//if (!input.is_open())
	//	std::cout << "\nNo file\n";
	input.close();
	input.open("vector_b.txt");
	b = new precision[n];
	for (int i = 0; i < n; i++)
		input >> b[i];
	input.close();
	input.open("matrix_al.txt");
	al = new precision * [n];
	for (int i = 0; i < n; i++)
	{
		al[i] = new precision[n];
		for (int j = 0; j < n; j++)
			input >> al[i][j];
	}
	input.close();
}

void CalculateDecomposition(precision**& al, precision*& b, const int& n, const int& p)
{
	precision l_i_k = 0;
	precision d_k_k = 0;
	precision l_j_k = 0;
	precision sum_over_k = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i + 1; j++)
		{
			if (i == j)
			{
				sum_over_k = 0;
				for (int k = 0; k < i; k++)
				{
					l_i_k = al[i][n - 1 - i + k];
					d_k_k = al[k][n - 1];
					sum_over_k += l_i_k * l_i_k * d_k_k;
				}
				al[i][n - 1] -= sum_over_k;
			}
			else
			{
				sum_over_k = 0;
				for (int k = 0; k < j; k++)
				{
					l_i_k = al[i][n - 1 - i + k];
					d_k_k = al[k][n - 1];
					l_j_k = al[j][n - 1 - j + k];
					sum_over_k += l_i_k * d_k_k * l_j_k;
				}
				al[i][n - 1 - i + j] = (al[i][n - 1 - i + j] - sum_over_k) / al[j][n - 1];
			}
		}
	}
}

void PrintVariables(precision** al, precision* b, const int& n)
{
	std::cout << std::setprecision(3) << "n: " << n;
	std::cout << "\nb: ";
	for (int i = 0; i < n; i++)
	{
		std::cout << b[i] << '\t';
	}
	std::cout << "\nal:\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			std::cout << al[i][j] << '\t';
		std::cout << "\n";
	}
	std::cout << "\n";
}

void DecompositionCheck(precision** al, const int& n)
{
	precision** a = new precision * [n];
	precision l_i_k = 0;
	precision d_k_j = 0;
	for (int i = 0; i < n; i++)
	{
		a[i] = new precision[n];
		for (int j = 0; j < n; j++)
		{
			a[i][j] = 0;
			for (int k = 0; k < n; k++)
			{
				l_i_k = i > k ? al[i][n - 1 - i + k] : (i < k ? 0 : 1);
				d_k_j = k == j ? al[k][n - 1] : 0;
				a[i][j] += l_i_k * d_k_j;
			}
		}
	}
	precision** b = new precision * [n];
	for (int i = 0; i < n; i++)
	{
		b[i] = new precision[n];
		for (int j = 0; j < n; j++)
		{
			b[i][j] = 0;
			for (int k = 0; k < n; k++)
			{
				l_i_k = k > j ? 0 : (k < j ? al[j][n - 1 - j + k] : 1);
				b[i][j] += a[i][k] * l_i_k;
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			std::cout << a[i][j] << '\t';
		std::cout << "\n";
	}
	std::cout << "\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			l_i_k = i > j ? 0 : (i < j ? al[j][n - 1 - j + i] : 1);
			std::cout << l_i_k << '\t';
		}
		std::cout << "\n";
	}
	std::cout << "\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			std::cout << b[i][j] << '\t';
		std::cout << "\n";
	}
	for (int i = 0; i < n; i++)
	{
		delete[] a[i];
		delete[] b[i];
	}
}

void SolveForX(precision** al, precision* b, precision*& x, const int& n)
{
	x = new precision[n];
	precision sum_over_k = 0;
	for (int i = 0; i < n; i++)
	{
		x[i] = 0;
		sum_over_k = 0;
		for (int k = 0; k < i; k++)
			sum_over_k += al[i][n - 1 - i + k] * x[k];
		x[i] = b[i] - sum_over_k;
	}

	for (int i = n - 1; i >= 0; i--)
	{
		sum_over_k = 0;
		for (int k = i + 1; k < n; k++)
			sum_over_k += al[k][n - 1 - k + i] * b[k];
		b[i] = x[i] / al[i][n - 1] - sum_over_k;
	}
	std::cout << "x: ";
	for (int i = 0; i < n; i++)
		std::cout << b[i] << ' ';
}
