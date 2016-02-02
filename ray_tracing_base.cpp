// ray_tracing_vs13.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <ctime>
#include <algorithm>
#include <cilk/cilk.h>
#include <windows.h>

using namespace std;
#include <fstream>
std::ifstream infile("src_det_pos.txt");

#if defined(AOS_CILK_FOR) || defined(SOA_CILK_FOR)
#include <cilk/cilk.h>
#endif

int main() {
	std::string line;
	int ctr = 0;
	int trial_ctr = 0, trial_max = 20;
	float *src_x = new float[126464];
	float *src_y = new float[126464];
	float *det_x = new float[126464];
	float *det_y = new float[126464];
	float *fp_ones = new float[126464];
	float *img = new float[241 * 400];
	clock_t clock_temp;
	double time_temp;
	double time_total = 0.0f, avg_time;

	int N_x, N_y;
	float delta_x, delta_y;
	//int img_ind, img_ind_max = 0;
	//int x_max, y_max;
	//float img_val, img_val_max = 0.0f;
	//int i_ind;

	cout << "Reading source-detector positions from the text file...";
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		if (!(iss >> src_x[ctr] >> src_y[ctr] >> det_x[ctr] >> det_y[ctr])) { break; } // error
		//cout << src_x[ctr] << " " << src_y[ctr] << " " << det_x[ctr] << " "  <<  det_y[ctr]  << endl;
		//cout << "ctr = " << ctr << endl;
		ctr++;
		// process pair (a,b)
	}
	cout << "Done!" << endl;




	// Run the whole thing 20 times and get the average time.
	for (trial_ctr = 0; trial_ctr < trial_max; trial_ctr++)
	{

		cout << "Initializing forward projection vector and the image...";
		for (int i = 0; i < 126464; i++)
			fp_ones[i] = 0.0f;
		for (int i = 0; i < 241 * 400; i++)
			img[i] = 1.0f;
		cout << "Done!" << endl;

		cout << "Starting ray_tracing..." << endl;

		N_x = 400;
		N_y = 240;
		delta_x = 400 / N_x;
		delta_y = 240 / N_y;


		unsigned __int64 cycle_temp = __rdtsc();
		double start = clock();
		cilk_for (int i = 0; i < 126464; i++)
		//for (int i = 0; i < 126464; i++)
		{
			float src_x_t, src_y_t, det_x_t, det_y_t;
			float L, l_total = 0.0f;
			float lambda_x_0, lambda_x_Nx, lambda_y_0, lambda_y_Ny;
			float lambda_min, lambda_max;
			//cout << "i = " << i << endl;
			src_x_t = src_x[i];
			src_y_t = src_y[i];
			det_x_t = det_x[i];
			det_y_t = det_y[i];
			// compute the length
			L = sqrt(pow(src_x_t - det_x_t, 2.0f) + pow(src_y_t - det_y_t, 2.0f));


			// find lambda values 
			lambda_x_0 = -L*src_x_t / (det_x_t - src_x_t);
			lambda_x_Nx = L*(N_x - src_x_t) / (det_x_t - src_x_t);
			lambda_y_0 = -L*src_y_t / (det_y_t - src_y_t);
			lambda_y_Ny = L*(N_y - src_y_t) / (det_y_t - src_y_t);
			lambda_min = max(max(0.0f, min(lambda_x_0, lambda_x_Nx)), min(lambda_y_0, lambda_y_Ny));
			lambda_max = min(min(L, max(lambda_x_0, lambda_x_Nx)), max(lambda_y_0, lambda_y_Ny));


			if ((lambda_min >= 0.0f) && (lambda_min <= L) && (lambda_max <= L) && (lambda_max >= 0) && (lambda_min <= lambda_max))
			{
				float x_min, y_min;
				int flag;
				float u_x, u_y;
				int ind_x, ind_y;
				int ind_x_old, ind_y_old;
				float lambda_x, lambda_y, lambda_x_next, lambda_y_next;
				float inc_x, inc_y;
				float l;

				//Compute the index of the first voxel    
				x_min = src_x_t + lambda_min*(det_x_t - src_x_t) / L;
				y_min = src_y_t + lambda_min*(det_y_t - src_y_t) / L;
				if (y_min <= 1e-13f)
					y_min = 0.0f;
				if (x_min <= 1e-13f)
					x_min = 0.0f;

				// Set the initial values of (lambda_x, lambda_y), the parameters where the
				// ray first intersects in the array with a plane perpendicular to the x, y,
				// z axis.
				if ((det_y_t >= src_y_t) && (det_x_t >= src_x_t))
				{
					flag = 1;
					u_x = ceil(x_min);
					u_y = ceil(y_min);
					ind_x = u_x;
					ind_y = u_y;
				}

				else if ((det_y_t < src_y_t) && (det_x_t >= src_x_t))
				{
					flag = 2;
					u_x = ceil(x_min);
					u_y = floor(y_min);
					ind_x = u_x;
					ind_y = u_y + 1;
				}

				else if ((det_y_t < src_y_t) && (det_x_t < src_x_t))
				{
					flag = 3;
					u_x = floor(x_min);
					u_y = floor(y_min);
					ind_x = u_x + 1;
					ind_y = u_y + 1;
				}

				else if ((det_y_t >= src_y_t) && (det_x_t < src_x_t))
				{
					flag = 4;
					u_x = floor(x_min);
					u_y = ceil(y_min);
					ind_x = u_x + 1;
					ind_y = u_y;
				}

				// 4 possible planes of entrance.
				if (u_x == 0) // Soldan giriyor
				{
					ind_x = 1;
					if (flag == 1) // ust saga gidiyor
						ind_y = u_y;
					else if (flag == 3) // alt saga gidiyor
						ind_y = u_y + 1;
				}

				else if (u_x == 400) // Sagdan giriyor
				{
					ind_x = 400;
					if (flag == 3) // alt sola gidiyor
						ind_y = u_y + 1;
					else if (flag == 4) // ust sola gidiyor
						ind_y = u_y;
				}

				if (u_y == 0) // Alttan giriyor
				{
					ind_y = 1;
					if (flag == 1) // ust saga gidiyor
						ind_x = u_x;
					else if (flag == 3) // ust sola gidiyor
						ind_x = u_x + 1;
				}

				else if (u_y == 240) // Ustten giriyor
				{
					ind_y = 240;
					if (flag == 1) // alt saga gidiyor
						ind_x = u_x;
					else if (flag == 3) // alt sola gidiyor
						ind_x = u_x + 1;
				}

				if ((u_x == 400) && (u_y == 0))
				{
					ind_x = 400;
					ind_y = 1;
				}

				lambda_x = L*(u_x - src_x_t) / (det_x_t - src_x_t);
				lambda_y = L*(u_y - src_y_t) / (det_y_t - src_y_t);




				inc_x = delta_x*abs(L / (det_x_t - src_x_t));
				inc_y = delta_y*abs(L / (det_y_t - src_y_t));

				while (lambda_x < lambda_min)
					lambda_x = lambda_x + inc_x;

				while (lambda_y < lambda_min)
					lambda_y = lambda_y + inc_y;

				lambda_x_next = lambda_x + inc_x;
				lambda_y_next = lambda_y + inc_y;

				if (det_x_t == src_x_t) // parallel in x plane
				{
					lambda_x = 1e23f;
					lambda_x_next = 1e24f;
				}

				if (det_y_t == src_y_t) // parallel in y plane
				{
					lambda_y = 1e23f;
					lambda_y_next = 1e24f;
				}

				//cout << "L = " << L << " u_x = " << u_x << " u_y = " << u_y <<  endl;
				//cout << " src_x = " << src_x_t << " src_y = " << src_y_t << " det_x = " << det_x_t << " det_y = " << det_y_t << endl;

				while ((lambda_x - lambda_max) < -1e-2f || (lambda_y - lambda_max) < -1e-2f)
				{
					//cout << "ind_x = " << ind_x << " ind_y = " << ind_y << endl;
					//cout << "lambda_x = " << lambda_x << "lambda_x_next = " << lambda_x_next <<  " lambda_y = " << lambda_y <<  " lambda_y_next = " << lambda_y_next << endl << "lambda_max = " << lambda_max << endl;
					//system("pause");
					ind_x_old = ind_x;
					ind_y_old = ind_y;
					if (lambda_x < lambda_y)
					{
						// Girisi lambda_x ten. Iki ihtimal - cikis ya lambda_y ya da
						// lambda_x_next
						// disp('Enter from x-plane');
						// x-plane, for example x = 0, parallel to y-axis
						if (lambda_x_next < lambda_y)
						{
							// disp('Exit from x-plane');
							// Cikis x'ten
							l = lambda_x_next - lambda_x;
							if (flag < 3)
								ind_x = ind_x + 1;
							else
								ind_x = ind_x - 1;
						}

						else
						{
							// disp('Exit from y-plane');
							// Cikis y'den
							l = lambda_y - lambda_x;
							if ((flag % 4) < 2)
								ind_y = ind_y + 1;
							else
								ind_y = ind_y - 1;
						}

						lambda_x = lambda_x_next;
						lambda_x_next = lambda_x + inc_x;
					}

					else
					{
						// disp('Enter from y-plane');
						// Girisi lambda_y ten. Iki ihtimal - cikis ya lambda_x ya da
						// lambda_y_next
						if (lambda_y_next < lambda_x)
						{
							// disp('Exit from y-plane');
							// Cikis y'den
							l = lambda_y_next - lambda_y;
							if ((flag % 4) < 2)
								ind_y = ind_y + 1;
							else
								ind_y = ind_y - 1;
						}
						else
						{
							// Cikis x'ten
							// disp('Exit from x-plane');
							l = lambda_x - lambda_y;
							if (flag < 3)
								ind_x = ind_x + 1;
							else
								ind_x = ind_x - 1;
						}

						lambda_y = lambda_y_next;
						lambda_y_next = lambda_y + inc_y;
					}

					//l_total += l*1.0f;
					//cout << "x_old = " << ind_x_old << " y_old = " << ind_y_old << endl;
					//cout << "lambda_x = " << lambda_x << " next = " << lambda_x_next << " lambda_y = " << lambda_y << " next = " << lambda_y_next << endl;
					l_total += l*img[(ind_x_old - 1) + (ind_y_old - 1) * 400];
				}
				fp_ones[i] = l_total;
				//cout << "Forward projection onto index " << i << " = " << l_total << endl;
				// system("pause");
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
	cout << "Ray tracing cpp | Cilk v1 | forward projection of ones | Average time = " << avg_time << endl;
	system("pause");
	return 0;
}