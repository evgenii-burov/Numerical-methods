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
				output << (p  - i > j ? 0 : int(10 * double(rand()) / RAND_MAX)+1) << '\t'; //(j + i + 1))
			//10*double(rand())/RAND_MAX)
			//(rand()/RAND_MAX + 3)*4
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
				l_i_k = i - k <= p ? (i!=k? al[i][p - i + k]:1) : 0; //al[i][p - i + k]
				d_k_k = al[k][p];
				l_j_k = j - k <= p ? (j != k ? al[j][p - j + k] : 1) : 0;
				sum_over_k += l_i_k * d_k_k * l_j_k;
			}
			al[i][p - i + j] = (al[i][p - i + j] - sum_over_k) / (i!=j ? (al[j][p]==0? 1: al[j][p]) : 1);
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
	std::cout << "\n\n";
}

void SolveForX(precision** al, precision* &b, const int& n, const int& p)
{
	precision l_i_k;
	precision l_k_i;
	precision sum_over_k = 0;
	for (int i = 0; i < n; i++)
	{
		sum_over_k = 0;
		for (int k = 0; k < i; k++)
		{
			l_i_k = i - k <= p ? (i!=k ? al[i][p - i + k] : 1) : 0;
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
		for (int k = i + 1; k < n; k++) // was k < n
		{
			l_k_i = k - i <= p ? (i!=k? al[k][p - k + i] : 1) : 0; //al[k][p - k + i]
			std::cout << '\n' << k << ' ' << i << ' ' << l_k_i << '\n';
			sum_over_k += l_k_i * b[k];
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
		std::cout << "\n\n";
	}
}

void CalculateError(precision**& al, precision* b, const int& n, const int& p)
{
	std::ifstream input;
	input.open("matrix_al.txt");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p + 1; j++)
			input >> al[i][j];
	}
	input.close();

	precision l_i_j = 0;
	precision elem = 0;
	std::cout << '\n';
	for (int i = 0; i < n; i++)
	{
		elem = 0;
		for (int j = 0; j < n; j++)
		{
			l_i_j = abs(i - j) <= p ? (i>j?al[i][p - i + j] : al[j][p - j + i]) : 0;
			std::cout << al[i][p - i + j] << ' ';
			elem += l_i_j * b[i];
		}
		std::cout << '\n' << elem << ' ';
	}
}