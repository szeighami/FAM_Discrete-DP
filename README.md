# FAM_Discrete-DP
This package contains the source code for the discrete utility function in two dimension for solving the FAM probelm described in (the paper only contains the continuous case): S. Zeighami and R. C.-W. Wong, “Finding average regret ratio minimizing set in database,” arXiv preprint: https://arxiv.org/abs/1810.08047

Prerequisite
===========
We have used "Boost library" (www.boost.org) in our program for random sampling of utility functions from a uniform distribution. Thus, the user should download this library and place it under the folder "ARR" (with the folder name "boost"). If the the library is in any other directory, replace the directory in -I option below to the desired directory. (i.e. if the boost file is in /home/myFiles/boost, then your compilation command should include -I /home/myFiles)

Usage Step
===========
a. Compilation
	g++ -o run -I . DP.cpp
b. Execution
	./run PointsFile k n N
Where PointsFile contains the all the points in the dataset, k is the size of the solution returned, n is the size of the dataset, and N is the sample size. The algorithm assumes a uniform distribution of linear utility funcitons in a two dimensional database and samples N utility functions from that distribution.
In the output file, you can see the results of the algorithm. 

===========
About dataset:
The dataset should contain n points in 2 dimensions. Each point must be written in a separate line with its dimensions being separated by tabs.

===========
Example (The hotel example from the paper):
The dataset contains 50 points in 2 dimensions. Use
	./run examplePoints.txt 2 50 100

to select 2 points from the 50 points in the 2-dimensional dataset using a sample size of 100 utility functions. The program outputs the results to the file result.txt. 
