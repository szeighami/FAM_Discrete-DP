#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstring>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/random_device.hpp>

using namespace std;

ifstream fileRead;
double* maxInData;

int memoryUsage = 0;

#ifndef WIN32
#include <sys/resource.h>
#include <sys/times.h>
#endif

struct points{
    int index;
    points* next;
};


#ifndef WIN32
void calculateExecutionTime(struct rusage *myTimeStart, struct rusage *myTimeEnd, double *userTime, double *sysTime)
{
        (*userTime) =
                ((double) (myTimeEnd->ru_utime.tv_sec  - myTimeStart->ru_utime.tv_sec)) +
                ((double) (myTimeEnd->ru_utime.tv_usec - myTimeStart->ru_utime.tv_usec)) * 1e-6;
        (*sysTime) =
                ((double) (myTimeEnd->ru_stime.tv_sec  - myTimeStart->ru_stime.tv_sec)) +
                ((double) (myTimeEnd->ru_stime.tv_usec - myTimeStart->ru_stime.tv_usec)) * 1e-6;

}
#endif



double gain(double* utility, double* data, int d)
{
	double gain = 0;
	for (int i = 0; i < d; i++)
	{
		gain += utility[i]*data[i];
	}
	return gain;
}

double maxInSolution(double* utility, double** data, int d, int* selected, int k)
{
	double max = 0;
	double gained = 0;
	for (int i = 0; i < k; i++)
	{
		gained = gain(utility, data[selected[i]], d);
		if (max < gained)
			max = gained;
	}
	return max;
}


double arr(int N, int n, int k, int d, double** data, double** utility, int* selected)
{
	double sum = 0;
	double temp = 0;
	double temp2 = 0;
	double temp3 = 0;
	for (int i = 0; i < N; i++)
	{

		temp2 = maxInSolution(utility[i], data, d, selected, k);
		temp3 = maxInData[i];
		temp = temp2/temp3;
		sum += temp;
	}

	return sum;
}

double calcArr(int N, int n, int k, int d, double** data, double** utility, string& bitmaskResult)
{
    string bitmask(k, 1); // K leading 1's
    bitmask.resize(n, 0); // N-K trailing 0's

    // print integers and permute bitmask
    int* selected = new int[k];
	
    int selectedIndex = 0;
    double max = 0;
    double arrValue = 0;
    int counter = 0;
    do {
	counter++;
        for (int i = 0; i < n; ++i) // [0..N-1] integers
        {
            if (bitmask[i])
            {
            	selected[selectedIndex] = i;
            	selectedIndex++;
            }
        }
        selectedIndex = 0;
        arrValue = arr(N, n, k, d, data, utility, selected);
        if (max < arrValue)
	{
		bitmaskResult = bitmask;
        	max = arrValue;
	}
    } while (prev_permutation(bitmask.begin(), bitmask.end()));
	delete[] selected;
	return max;
}



double getOnePointARR(double** data, double** utility, double* maxInData,int i, int j, int N, int n, int pointIndex)
{
    double regretRatioSum = 0;
    for (int a = i; a <= j; a++)
    {
        double pointUtil = 0;
        for (int d = 0; d < 2; d++)
            pointUtil += data[pointIndex][d] * utility[a][d];

        regretRatioSum += 1 - (pointUtil / maxInData[a]);
    }
    return regretRatioSum / (j-i+1);
}

double findOnePointARR(double** data, double** utility, double* maxInData,int i, int j, int N, int n)
{
    double minARR = 1;
    for (int pointIndex = 0; pointIndex < n; pointIndex++)
    {
        double arr = getOnePointARR(data, utility, maxInData, i, j, N, n, pointIndex);
        if (arr < minARR)
            minARR = arr;
    }
    return minARR;
}

double DP(double** utility, double** data, int k, int N, int n, double* maxInData)
{
   double*** S = new double**[k];
   for (int r = 0; r < k; r++) 
   {
       S[r] = new double*[N];
       for (int i = 0; i < N; i++)
       {
           S[r][i] = new double[N];
           for (int j = i; j < N; j++)
           {
               if (r == 0)
               {
                    S[r][i][j] = findOnePointARR(data, utility, maxInData, i, j, N, n);       
                    continue;
               }

               if (i == j)
               {
                   S[r][i][j] = S[0][i][j];
                   continue;
               }

               double min = S[0][i][i] + (j-i)*S[r-1][i+1][j];
               for (int a = i + 1; a <= j; a++)
               {
                   double value = (a-i+1)*S[0][i][a] + (j-a)*S[r-1][(a+1>j)?j:(a+1)][j]; 
                   if (value < min)
                       min = value; 
               }
               S[r][i][j] = min / (j-i+1);
               
           }
       }
   }
 
   return S[k-1][0][N-1];

}

int comp (const void * elem1, const void * elem2) 
{
    double* f = *((double**)elem1);
    double* s = *((double**)elem2);
    if (f[1]/f[0] > s[1]/s[0]) return  1;
    if (f[1]/f[0] < s[1]/s[0]) return -1;
    return 0;
}

int main(int argc, char** argv)
{
#ifndef WIN32
	double  userTime, sysTime;
	struct rusage myTime_start, myTime_end;
	double  userTime2, sysTime2;
	struct rusage myTime_start2, myTime_end2;
    double  userTime3, sysTime3;
	struct rusage myTime_start3, myTime_end3;
#endif

    boost::mt19937 rng(time(0));
	boost::uniform_real<double> u(0, 1);
	boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > gen(rng, u);

	int N = 0, n = 10000, k = 0, d = 2;
	memoryUsage = 0;
    string filename = argv[1];
	k = atoi(argv[2]);
	n = atoi(argv[3]);
	N = atoi(argv[4]);
	ifstream inFile;
	inFile.open(filename.c_str());
	
	double** data = new double*[n];
	memoryUsage += sizeof(double) * n;
	double** utility = new double*[N];
	memoryUsage += sizeof(double) * N;
	maxInData = new double[N];
	memoryUsage += sizeof(double) * N;

#ifndef WIN32
	// start time
	getrusage(RUSAGE_SELF,&myTime_start);
#endif

    //Reading points
	for(int i = 0; i < n; i++)
	{
		data[i] = new double[d];
		memoryUsage += sizeof(double) * d;
		for(int j = 0; j < d; j++)
		{
			inFile >> data[i][j];
		}
	}
	inFile.close();

    //Samling Utility functions
	for (int i = 0; i < N; i++)
	{
		utility[i] = new double[d];
		memoryUsage += sizeof(double) * d;
		for(int j = 0; j < d; j++)
		{
			utility[i][j] = gen();
		}
		maxInData[i] = 0;				
		for (int j = 0; j < n; j++)
		{
			double util = 0;
			for(int z = 0; z < d; z++)
			{
				util += utility[i][z] * data[j][z];
			}
			if (util > maxInData[i])
			{
				maxInData[i] = util;
			}
		}
	}

    FILE * fp;
    fp = fopen ("result.txt", "a");
	fprintf(fp, "value of k = %d\n", k);

#ifndef WIN32
	// end time
	getrusage(RUSAGE_SELF,&myTime_end);
	// output execution time
	calculateExecutionTime(&myTime_start, &myTime_end, &userTime, &sysTime);
	fprintf(fp, "User time preprocessing: %f seconds\n", userTime);
	fprintf(fp, "System time preprocessing: %f seconds\n\n", sysTime);
	// start time
	getrusage(RUSAGE_SELF,&myTime_start2);
#endif



	cout << "k:" << k << endl;

    //Finding Best points
	string result = "";
    qsort (utility, N, sizeof(double*), comp);
    for (int i = 0; i < N; i++)
	{
		maxInData[i] = 0;				
		for (int j = 0; j < n; j++)
		{
			double util = 0;
			for(int z = 0; z < d; z++)
			{
				util += utility[i][z] * data[j][z];
			}
			if (util > maxInData[i])
			{
				maxInData[i] = util;
			}
		}
	}
    //Running DP
    cout << "running DP" << endl;
	double DPSol = DP(utility, data, k, N, n, maxInData); 

 #ifndef WIN32
	// end time
	getrusage(RUSAGE_SELF,&myTime_end2);
	// output execution time
	calculateExecutionTime(&myTime_start2, &myTime_end2, &userTime2, &sysTime2);
	fprintf(fp, "User time DP: %f seconds\n", userTime2);
	fprintf(fp, "System time DP: %f seconds\n\n", sysTime2);
	getrusage(RUSAGE_SELF,&myTime_start3);
#endif

    cout << "running brute force" << endl;
    double BruteForceSol = (N - calcArr(N, n, k, d, data, utility, result))/N;

#ifndef WIN32
	// end time
	getrusage(RUSAGE_SELF,&myTime_end3);
	// output execution time
	calculateExecutionTime(&myTime_start3, &myTime_end3, &userTime3, &sysTime3);
	fprintf(fp, "User time Brute-Force: %f seconds\n", userTime3);
	fprintf(fp, "System time Brute-Force: %f seconds\n\n", sysTime3);
#endif

	cout << DPSol << endl;
	cout << BruteForceSol << endl;

	fprintf(fp, "memory : %d\n", memoryUsage);
	fprintf(fp, "average regret ratio (DP): %f\n", DPSol);
	fprintf(fp, "average regret ratio (Brute force): %f\n", BruteForceSol);
	fprintf(fp, "-------------------------\n");


	cout << "end" << endl;

	return 0;
}


