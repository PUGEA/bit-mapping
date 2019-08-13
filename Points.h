#pragma once
#include <iostream>
#include <fstream>

#include "Parameter.hpp"
#include "Utils.hpp"
#include "MemoryMapped.h"
#include "Common.h"

class Points
{
public :
	int nP;
	int dim;
	REAL_TYPE **d;
	MemoryMapped all_reads;

	void A_Read(int i,REAL_TYPE *read)
	{
		for (int j = 0; j < dim; ++j)
		{
			read[j] = d[i][j];
		//	read[j] = ictoi_table[all_reads[i * (dim + 1) + j]];
		}
	}

	void Initialize_MemoryMapped(std::string filename )
	{
		
		MemoryMapped reads_file(filename);
	/*	for (int i = 0; i < nP; ++i)
		{
			for (int j = 0; j < dim; ++j)
			{
			//	std::cout << i * (dim + 1) + j << std::endl;
				d[i][j] = ictoi_table[reads_file[i * (dim + 1) + j]];
			}
		//	std::cout << std::endl;
		}
		reads_file.close();*/
	}

	void Initialize(int _nP, int _dim)
	{
		nP = _nP;
		dim = _dim;
		d = new REAL_TYPE * [ nP ];
		for(int i=0;i<nP;i++)
		{
			d[i] = new REAL_TYPE [ dim ];
		}
	}

	// this function is for read point set from file
	// format:
	// (num. points) (dimensionality)
	// v0 (floats, number of elements is equal to dimensionality)
	// v1
	// ...
	void Initialize_From_File(std::string filename, int num = 0)
	{
		std::ifstream srcfile(filename);
		std::string tmp;
		if (!srcfile.is_open())
		{
			perror(INPUT_REF_FILE_NAME);
			exit(EXIT_FAILURE);
		}
		int cond = 0;
		if (num != 0)
		{
			cond = num;
		}else
		{
			cond = nP;
		}
		for(int i=0;i < cond;i++)
		{
			std::getline(srcfile,tmp);
			for(int k=0;k<dim;k++)
			{
				d[i][k] = ictoi_table[tmp[k]];
			//	std::cout << d[i][k];
			}
		//	std::cout << std::endl;
		}
	}

	// computing center of points for zero centering
	void Compute_Center(REAL_TYPE *center)
	{
		double *tCenter = new double [dim];
		SetVector_Val<double>( tCenter , dim , 0.0 );
		for(int i=0;i<nP;i++)
		{
			for(int k=0;k<dim;k++)
			{
				tCenter[k] += d[i][k];
			}
		}
		for(int k=0;k<dim;k++)
		{
			tCenter[k] /= (double)(nP);
			center[k] = (REAL_TYPE)( tCenter[k] );
		}
		delete [] tCenter;
	}

	void ReleaseMem()
	{
		for (int i = 0; i < nP; ++i)
	    {
	        delete [] d[i];
	    }
	    delete [] d;
	    d = NULL;
		//all_reads.close();
	}

	
};
