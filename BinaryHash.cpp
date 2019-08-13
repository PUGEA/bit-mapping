#include "BinaryHash.h"

#include <omp.h>
#include <iostream>

void Sphere::Set_Radius(Points *ps, Index_Distance *ids)
{
	for(int i=0;i<ps->nP;i++)
	{
		ids[i].index = i;
		ids[i].distSq = Compute_Distance_L2Sq<REAL_TYPE>( c , ps->d[i] , ps->dim );
	//	ids[i].distSq = Compute_Edit_Distance<REAL_TYPE>( c , ps->d[i] , ps->dim );
		ids[i].dist = sqrt( ids[i].distSq );
	}
	sort( &ids[0] , &ids[ ps->nP ] );
	int tIndex = (int)( (REAL_TYPE)( ps->nP ) * INCLUDING_RATIO );
	r = ids[tIndex-1].dist;
	rSq = r * r;
}

void Sphere::Initialize(int _dim)
{
	c = new REAL_TYPE [ _dim ];
	r = 0.0;		rSq = 0.0;
}

void SphericalHashing::Initialize(Points *_ps,int _code_len,int _seg_len)
{
	ps = _ps;
	code_len = _code_len;
	dim = _seg_len;

	tps.Initialize( NUM_TRAIN_SAMPLES , ps->dim );	
	bool *checkList = new bool [ ps->nP ];
	SetVector_Val<bool>( checkList , ps->nP , true );
	int index;

	// sampling training set
	for(int i=0;i<NUM_TRAIN_SAMPLES;i++)
	{
		while(true)
		{
			index = Rand_Uniform_Int( 0 , ps->nP - 1 );
			if( checkList[index] )
			{
				checkList[index] = false;
				break;
			}
		}
		SetVector_Vec<REAL_TYPE>( tps.d[i] , ps->d[index] , tps.dim );		
	}
	delete [] checkList;
	s = new Sphere [ code_len ];
	for(int i=0;i<code_len;i++)
	{
		s[i].Initialize( tps.dim );
	}
	table = new bitset<NUM_TRAIN_SAMPLES> [ code_len ];
	for(int i=0;i<code_len;i++)
	{
		for(int j=0;j<NUM_TRAIN_SAMPLES;j++)
		{
			table[i][j] = 0;
		}
	}
	ids = new Index_Distance * [ code_len ];
	for(int i=0;i<code_len;i++)
	{
		ids[i] = new Index_Distance [ NUM_TRAIN_SAMPLES ];
	}
}


void SphericalHashing::Save_Sphere_Info()
{
	std::ofstream file;
	file.open("bin/sphere_info.log");	
	for (int i = 0; i < code_len; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			file << s[i].c[j] << " " ;
		//	std::cout << s[i].c[j] << "\t" ;
		}
		file << s[i].rSq << std::endl;
	//	std::cout << s[i].rSq << std::endl;
	}
	file.close();
}

void SphericalHashing::Load_Sphere_Info(int _code_len,int _seg_len)
{
	code_len = _code_len;
	dim = _seg_len;

	std::ifstream file;
	file.open("bin/sphere_info.log");
	if (!file.is_open())
	{
		perror("bin/sphere_info.log");
		exit(0);
	}
	std::string tmp;

	s = new Sphere [ code_len ];
	for (int i = 0; i < code_len; ++i)
	{
		s[i].Initialize( dim );

		getline(file,tmp);
		std::stringstream ss(tmp);
		for (int j = 0; j < dim; ++j)
		{
			ss >> s[i].c[j];
		//	std::cout << s[i].c[j] << "\t";
		}
		ss >> s[i].rSq;
	//	std::cout << s[i].rSq << std::endl;
	}
	file.close();
}

void SphericalHashing::ReleaseMem()
{
	tps.ReleaseMem();
	delete [] table;
	for(int i=0;i<code_len;i++)
	{
		delete [] ids[i];
	}
	delete ids;
}

void SphericalHashing::Compute_Table()
{
	for(int i=0;i<code_len;i++)
	{
#ifdef USE_PARALLELIZATION
		#pragma omp parallel for
#endif
		for(int j=0;j<NUM_TRAIN_SAMPLES;j++)
		{
			table[i][j] = 0;
			if( Compute_Distance_L2Sq<REAL_TYPE>( s[i].c , tps.d[j] , tps.dim ) < s[i].rSq )
			//if( Compute_Edit_Distance<REAL_TYPE>( s[i].c , tps.d[j] , tps.dim ) < s[i].rSq )
			{
				table[i][j] = 1;
			}
		}
	}
}

void SphericalHashing::Compute_Num_Overlaps(int **overlaps)
{
	for(int i=0;i<code_len-1;i++)
	{
		overlaps[i][i] = table[i].count();
#ifdef USE_PARALLELIZATION
		#pragma omp parallel for
#endif
		for(int j=i+1;j<code_len;j++)
		{
			overlaps[i][j] = ( table[i] & table[j] ).count();
			overlaps[j][i] = overlaps[i][j];
		}
	}
}



void SphericalHashing::Set_Spheres()
{
	REAL_TYPE	marginT = 0.05;

	double allowedErrorMean, allowedErrorVar;
	allowedErrorMean = (double)tps.nP * OVERLAP_RATIO * EPSILON_MEAN;
	allowedErrorVar = (double)tps.nP * OVERLAP_RATIO * EPSILON_STDDEV;

	// initial pivots are determined as center of 10 randomly sampled points
	// for locating initial pivots near the data center
	int rIndex;
	for(int i=0;i<code_len;i++)
	{
		SetVector_Val<REAL_TYPE>( s[i].c , tps.dim , 0.0 );
		for(int j=0;j<10;j++)
		{
			rIndex = Rand_Uniform_Int( 0 , tps.nP - 1 );
			Add_Vector<REAL_TYPE>( s[i].c , tps.d[ rIndex ] , s[i].c , tps.dim );
		}
		Scalar_Vector<REAL_TYPE>( s[i].c , 0.1 , tps.dim );
	}

#ifdef USE_PARALLELIZATION
	#pragma omp parallel for
#endif
	for(int i=0;i<code_len;i++)
	{
		s[i].Set_Radius( &tps , &ids[i][0] );
	}
	
	int **overlaps = new int * [ code_len ];
	for(int i=0;i<code_len;i++)
	{
		overlaps[i] = new int [ code_len ];
	}
	REAL_TYPE **forces = new REAL_TYPE * [ code_len ];
	for(int i=0;i<code_len;i++)
	{
		forces[i] = new REAL_TYPE [ tps.dim ];
	}
	REAL_TYPE *force = new REAL_TYPE [ tps.dim ];
	REAL_TYPE tmpOverlap, alpha;

	for(int k=0;k<MAX_NUM_ITERATIONS;k++)
	{
		Compute_Table();
		Compute_Num_Overlaps( overlaps );

		double mean, variance, cnt;
		mean = 0.0;		cnt = 0.0;		variance = 0.0;
		for(int i=0;i<code_len-1;i++)
		{
			for(int j=i+1;j<code_len;j++)
			{
				mean += (double)( overlaps[i][j] );
				cnt += 1.0;
			}
		}
		mean /= cnt;
		for(int i=0;i<code_len-1;i++)
		{
			for(int j=i+1;j<code_len;j++)
			{
				variance += ( (double)( overlaps[i][j] ) - mean ) * ( (double)( overlaps[i][j] ) - mean );
			}
		}
		variance /= cnt;
		
		// iteration convergence condition
		if( fabs( mean - ( (double)tps.nP * OVERLAP_RATIO ) ) < allowedErrorMean && sqrt(variance) < allowedErrorVar )
		{
			break;
		}
		

		// force computation
		SetMatrix_Val<REAL_TYPE>( forces , code_len , tps.dim , 0.0 );
		for(int i=0;i<code_len-1;i++)
		{
			for(int j=i+1;j<code_len;j++)
			{
				tmpOverlap = (REAL_TYPE)overlaps[i][j] / (REAL_TYPE)tps.nP;
				alpha = ( tmpOverlap - OVERLAP_RATIO ) / OVERLAP_RATIO;
				alpha /= 2.0;

				Sub_Vector<REAL_TYPE>( s[j].c , s[i].c , force , tps.dim );
				Scalar_Vector<REAL_TYPE>( force , alpha , tps.dim);
				Add_Vector<REAL_TYPE>( forces[j] , force , forces[j] , tps.dim );
				Scalar_Vector<REAL_TYPE>( force , -1.0 , tps.dim );
				Add_Vector<REAL_TYPE>( forces[i] , force , forces[i] , tps.dim );
			}
		}

		// move pivots and adjust radius
#ifdef USE_PARALLELIZATION
		#pragma omp parallel for
#endif
		for(int i=0;i<code_len;i++)
		{
			Scalar_Vector<REAL_TYPE>( forces[i] , 1.0 / code_len , tps.dim );
			Add_Vector<REAL_TYPE>( s[i].c , forces[i] , s[i].c , tps.dim );
			s[i].Set_Radius( &tps , &ids[i][0] );
		}
	}

	for(int i=0;i<code_len;i++)
	{
		delete [] overlaps[i];
		delete [] forces[i];
	}
	delete [] overlaps;
	delete [] forces;
	delete [] force;
}