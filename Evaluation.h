#include <iostream>

#pragma once

template<typename DistType>
class Result_Element
{
public :
	int index;
	DistType dist;
	bool isCorrect;
	
	bool operator < (const Result_Element<DistType> &T) const
	{
		if( this->dist < T.dist )
		{
			return true;
		}
		return false;
	}
};


// compute average precision for each query
template<typename DistType>
double Compute_AP(int *gt, Result_Element<DistType> *res, int nP, Points dps)
//double Compute_AP( Result_Element<DistType> *res, int nP, Points dps)
{
	double ret = 0.0;
#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<nP;i++)
	{
		res[i].isCorrect = false;
	}

#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif

/*	for(int i=0;i<K;i++)
	{
		res[ gt[i] ].isCorrect = true;
		//-----------------
		std::cout << "true data:";
		for (int j = 0; j < 50; ++j)
		{
			std::cout << dps.d[gt[i]][j];
		}
		std::cout << std::endl;
		//----------------
	}*/
	sort( &res[0] , &res[nP] );

	for(int i=0;i<KNN;i++)
	{
		if (res[i].dist == res[0].dist)
		{
			//-----------------
			std::cout << "hash data " << i + 1 << ": ";
			for (int j = 0; j < 50; ++j)
			{
				std::cout << dps.d[res[i].index][j];
			}
			std::cout << std::endl;
			//----------------
		}
	}
	return ret;
}