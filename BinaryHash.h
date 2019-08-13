#pragma once

#include "Common.h"
#include "Points.h"
#include <fstream>
#include <sstream>

// hamming distance function
__inline int Compute_HD(bitset<BCODE_LEN> &a, bitset<BCODE_LEN> &b)
{
	return ( ( a ^ b ).count() );
}

// spherical hamming distance function
__inline double Compute_SHD(bitset<BCODE_LEN> &a, bitset<BCODE_LEN> &b)
{
	return ( ( (double)( ( a ^ b ).count() ) ) / ( (double)( ( a & b ).count() ) + 0.1 ) );
}

class Index_Distance
{
public :
	bool flag;
	int index;
	REAL_TYPE dist, distSq;
	bool operator < (const Index_Distance &T) const
	{
		if( this->distSq < T.distSq )	{			return true;		}
		return false;
	}
};

class Sphere
{
public :
	REAL_TYPE *c, r, rSq;
	
	void Initialize(int _dim);

	// function to set radius to include desired portion of training set
	void Set_Radius(Points *ps, Index_Distance *ids);
};

class SphericalHashing
{
public :
	Points *ps;

	// training set
	Points tps;

	Sphere *s;

	int code_len;
	int dim;

	Index_Distance **ids;
	bitset<NUM_TRAIN_SAMPLES> *table;

	void Initialize(Points *_p,int code_len,int dim);
	void Compute_Table();
	void Compute_Num_Overlaps(int **overlaps);
	void Set_Spheres();

	void ReleaseMem();

	template<typename argType, typename bitType>
	__inline void compute(argType *x, bitType &y, int start, int idx, int len)
	{
	#ifdef USE_PARALLELIZATION
		#pragma omp parallel for
	#endif
		for(int i = 0;i < len;i++)
		{
			argType dis = Compute_Distance_L2Sq<REAL_TYPE>( s[i + idx].c , x , dim , start);
			if( dis - s[i].rSq > 0.000001)
			{
				y.set(i, 0);
			}
			else
			{
				y.set(i, 1);
			}
		//	std::cout << dis << " " << s[i].rSq << " " << y[i] << std::endl;
		//	start += dim;
			start += 1;
		}
		//std::cout << y << std::endl;
	}

	template<typename argType>
	__inline void Compute_BCode_Read(argType *x, bitset<BCODE_LEN> &y)
	{
		int start;
		
		start = Parameter::region_searching + Parameter::skip;	
	
		compute<argType, bitset<BCODE_LEN>>(x, y, start, 0, Parameter::bcode_len);
	}

	template<typename argType>
	__inline void Compute_BCode_Ref(argType *x, bitset<BCODE_64> &y, int idx)
	{
		int start;
	
		start = Parameter::region_searching + idx;
		
		compute<argType, bitset<BCODE_64>>(x, y, start, idx, BCODE_64);
	}

	void Save_Sphere_Info();

	void Load_Sphere_Info(int _code_len,int _seg_len);
};
