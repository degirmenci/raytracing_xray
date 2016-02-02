// ray_tracing_vs13.cpp : Defines the entry point for the console application.
//



#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <ctime>
#include <algorithm>
#include <omp.h>
#include <windows.h>
#include <cilk/cilk.h>

using namespace std;
#include <fstream>

#if defined(AOS_CILK_FOR) || defined(SOA_CILK_FOR)
#include <cilk/cilk.h>
#endif

int main() {
	std::string line;
	int ctr = 0;
	int trial_ctr = 0, trial_max = 20;
	// number of nonzero elements = 43864476
	// i_ptr_start: 160
	int *i_offset = new int[126464];
	int *j_list = new int[43864476];
	int i_start = 160;
	float *val_list = new float[43864476];
	float *fp_ones = new float[126464];
	float *img = new float[240 * 400];
	clock_t clock_temp;
	double time_temp;
	double time_total = 0.0f, avg_time;

	int N_x, N_y;
	float delta_x, delta_y;
	//int img_ind, img_ind_max = 0;
	//int x_max, y_max;
	//float img_val, img_val_max = 0.0f;
	//int i_ind;

	cout << "Reading j_list and value_list from the text file...";
	std::ifstream infile("psf_2d_j_list_val_list.txt");
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		if (!(iss >> j_list[ctr] >> val_list[ctr])) { break; } // error
		//cout << src_x[ctr] << " " << src_y[ctr] << " " << det_x[ctr] << " "  <<  det_y[ctr]  << endl;

		if (ctr % 1000000 == 0)
			cout << "ctr = " << ctr << endl;

		ctr++;
		// process pair (a,b)
	}

	cout << "Done!" << endl;

	cout << "Reading i_ptr_offset from the text file...";
	std::ifstream infile2("psf_2d_i_offset_list.txt");
	ctr = 0;
	while (std::getline(infile2, line))
	{
		std::istringstream iss(line);
		if (!(iss >> i_offset[ctr])) { break; } // error
		//cout << src_x[ctr] << " " << src_y[ctr] << " " << det_x[ctr] << " "  <<  det_y[ctr]  << endl;
		//cout << "ctr = " << ctr << endl;
		if (ctr % 1000000 == 0)
			cout << "ctr = " << ctr << endl;

		ctr++;
		// process pair (a,b)
	}
	cout << "Done!" << endl;


	// Run the whole thing 20 times and get the average time.
	for (trial_ctr = 0; trial_ctr < trial_max; trial_ctr++)
	{

		cout << "Initializing the image...";

		for (int i = 0; i < 240 * 400; i++)
			img[i] = 1.0f;
		cout << "Done!" << endl;

		cout << "Starting sparse vector multiplication..." << endl;



		unsigned __int64 cycle_temp = __rdtsc();
		double start = clock();


		// SPARSE MATRIX VECTOR MULTIPLICATION 
		// SOYSAL DEGIRMENCI
		// 04/20/15

		//for (i = 0; i<126464; i++)
		int i;
		cilk_for (i = 0; i<126464; i++)
		{
			fp_ones[i] = 0;
			//cout << "i = " << i << " i_offset[i] = " << i_offset[i] << " i_offset[i+1] = " << i_offset[i + 1] << endl;
			int ckey;
			for (ckey = i_offset[i]; ckey<i_offset[i + 1]; ckey++)
			{
				int j = j_list[ckey];
				//cout << "ckey = " << ckey << " img[j] = " << img[j] << endl;
				fp_ones[i] += val_list[ckey] * img[j];
			}
		}

		double end = clock();
		clock_temp = __rdtsc() - cycle_temp;
		printf_s("%I64d ticks\n", clock_temp);
		time_temp = (end - start) / CLOCKS_PER_SEC;
		cout << "Trial " << trial_ctr << " Time = " << time_temp << endl;
		float fp_ones_max = 0.0f;

		for (int i = 0; i<126464; i++)
		{
			if (fp_ones[i]>fp_ones_max)
				fp_ones_max = fp_ones[i];
		}
		cout << "fp_ones_max = " << fp_ones_max << endl;
		time_total += time_temp;
		//system("pause");
	}
	avg_time = time_total / trial_max;
	cout << "Ray tracing cpp | OMP v1 | forward projection of ones | Average time = " << avg_time << endl;
	system("pause");
	return 0;
}