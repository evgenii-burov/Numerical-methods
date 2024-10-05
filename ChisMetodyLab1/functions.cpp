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
		for (int j = 0; j < i+1; j++)
		{
			sum_over_k = 0;
			for (int k = 0; k < j; k++)
			{
				l_i_k = i-k <= p ? al[i][p - i + k] : 0;
				d_k_k = al[k][p];
				l_j_k = j - k <= p ? al[j][p - j + k] : 0;
				sum_over_k += l_i_k * d_k_k * l_j_k;
			}
			al[i][p - i + j] = (al[i][p - i + j] - sum_over_k) / (i!=j ? al[j][p] : 1);
		}
	}
}

void PrintVariables(precision** al, precision* b, const int& n,const int& p)
{
	std::cout << std::setprecision(3) << "n: " << n << "p: " << p;
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
	std::cout << "\n";
}

void SolveForX(precision** al, precision* b, const int& n, const int& p)
{
	precision sum_over_k = 0;
	precision l_i_k;
	for (int i = 0; i < n; i++)
	{
		sum_over_k = 0;
		for (int k = 0, jl = n-1-i; k < i; k++, jl++)
		{
			l_i_k = abs(i - k) <= p ? al[i][jl] : 0;
			sum_over_k += l_i_k * b[k];
		}
		b[i] = b[i] - sum_over_k;
	}

	std::cout << "x: ";
	for (int i = 0; i < n; i++)
		std::cout << b[i] << ' ';
	for (int i = n - 1; i >= 0; i--)
	{
		sum_over_k = 0;
		//pamagite
		for (int k = i + 1, jl = n-1; k < n; k++, jl--)
		{
			l_i_k = abs(i - k) <= p ? al[k][n - 1 - k + i] : 0;
			sum_over_k += l_i_k * b[k];
		}
		b[i] = b[i] / al[i][p] - sum_over_k;
	}
	std::cout << "x: ";
	for (int i = 0; i < n; i++)
		std::cout << b[i] << ' ';
}

void PrintFullMatrix(precision** al, const int& n, const int& p)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			std::cout << (i-j<=p ? al[i][p - i + j] : 0)<<'\t';
		}
		std::cout << '\n';
	}
}