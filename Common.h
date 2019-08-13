#pragma once

#include<iostream>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <bitset>

#include "Parameters.h"
#include "Stopwatch.hpp"

using namespace std;

#define		PII				3.1415
#define MIN(a,b,c)  (a<b?(a<c?a:c):(b<c?b:c))

int Rand_NDigits(int nDigits);
int Rand_Uniform_Int(int minVal, int maxVal);

// transfer string to bitset
template<size_t N>
bitset<N> String2Bit(const std::string &s)
//bitset<N> from_string(const unsigned long& s)
{
    return bitset<N>(s);
}

template<typename RealT>
RealT Rand_Uniform(RealT minVal, RealT maxVal)
{
	return ( minVal + ( maxVal - minVal ) * ( (RealT)rand() / (RealT)(RAND_MAX) ) );
}

template<typename RealT>
RealT Rand_Gaussian()
{
	RealT x1, x2, ret;
	do
	{
		x1 = Rand_Uniform<RealT>( 0.0 , 1.0 );
	} while( x1 == 0 );
	x2 = Rand_Uniform<RealT>( 0.0 , 1.0 );
	ret = sqrt( -2.0 * log( x1 ) ) * cos( 2.0 * PII * x2 );
	return ret;
}


template<typename RealT>
int Compute_Edit_Distance(RealT *a, RealT *b,int dim)
{
    int d[dim+1];
    int i, j, old, temp;

    for (j = 0; j <= dim; j++) 
    {
        d[j] = j;
    }

    for (i = 1; i <= dim; i++) 
    {
        old = i - 1;
        d[0] = i;
        for (j = 1; j <= dim; j++) 
        {
            temp = d[j];
            // 算法中 a, b 字符串下标从 1 开始，c 语言从 0 开始，所以 -1
            if (a[i-1] - b[j-1] <= 0.5) 
            {
                d[j] = old;
            } else 
            {
                d[j] = MIN(d[j] + 1, d[j-1] + 1, old + 1);
            }
            old = temp;
        }
    }

    return d[dim];
}


template<typename RealT>
RealT Compute_Distance_L2Sq(RealT *v0, RealT *v1, int dim, int start)
{
	RealT ret = 0.0;
	for(int i=0;i<dim;i++)
	{
		//std::cout << v0[i] << "-" << v1[start + i] << "=" << v0[i] - v1[start + i] << std::endl;
		ret += ( v0[i] - v1[start + i] ) * ( v0[i] - v1[start + i] );
	}
	
	return ret;
}

template<typename RealT>
RealT Compute_Distance_L2Sq(RealT *v0, RealT *v1, int dim)
{
	RealT ret = 0.0;
	for(int i=0;i<dim;i++)
	{
		ret += ( v0[i] - v1[i] ) * ( v0[i] - v1[i]);
	}
	return ret;
}

template<typename RealT>
void SetVector_Val(RealT *vec, int dim, RealT val)
{
	for(int i=0;i<dim;i++)
	{
		vec[i] = val;
	}
}

template<typename RealT>
void SetVector_Vec(RealT *dest, RealT *src, int dim)
{
	for(int i=0;i<dim;i++)
	{
		dest[i] = src[i];
	}
}

template<typename RealT>
void Scalar_Vector(RealT *A, RealT s, int dim)
{
	for(int i=0;i<dim;i++)
	{
		A[i] = s * A[i];
	}
}

template<typename RealT>
void Sub_Vector(RealT *A, RealT *B, RealT *ret, int dim)
{
	for(int i=0;i<dim;i++)
	{
		ret[i] = A[i] - B[i];
	}
}


template<typename RealT>
void Add_Vector(RealT *A, RealT *B, RealT *ret, int dim)
{
	for(int i=0;i<dim;i++)
	{
		ret[i] = A[i] + B[i];
	}
}

template<typename RealT>
void SetMatrix_Val(RealT **mat, int nY, int nX, RealT val)
{
	for(int i=0;i<nY;i++)
	{
		for(int j=0;j<nX;j++)
		{
			mat[i][j] = val;
		}
	}
}