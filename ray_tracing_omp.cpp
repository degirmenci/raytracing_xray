// ray_tracing_vs13.cpp : Defines the entry point for the console application.
//
// do 3d version for this one.
//   /Qvec-report:1 

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
#include <omp.h>
#include <immintrin.h>
#include <vector>
#include <malloc.h>


#define CLOCKS_PER_SEC 1000

using namespace std;
#include <fstream>
std::ifstream infile("src_det_pos_2d_v1.txt");

#if defined(AOS_CILK_FOR) || defined(SOA_CILK_FOR)
#include <cilk/cilk.h>
#endif

static int DATA_SIZE = 62976;
static int IMG_SIZE = 96000;
static int NUM_SLICES = 160;
static int IMG_SIZE_2 = 160 * 240;

int main() {
	std::string line;
	int ctr = 0;
	int trial_ctr = 0, trial_max = 20;
	float *src_x = new float[DATA_SIZE];
	float *src_y = new float[DATA_SIZE];
	float *det_x = new float[DATA_SIZE];
	float *det_y = new float[DATA_SIZE];
	float *fp_ones = new float[DATA_SIZE*NUM_SLICES];
	//float *img = new float[IMG_SIZE*NUM_SLICES];
	//std::vector<float> img(IMG_SIZE*NUM_SLICES);
	//__declspec(align(32)) float img[96000*160];
	//float *img = (float *)malloc(96000 * 160);
	float *img = (float *)_aligned_malloc(96000 * 160 * sizeof(float), 32);
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
		for (int i = 0; i < DATA_SIZE*NUM_SLICES; i++)
			fp_ones[i] = 0.0f;
		for (int i = 0; i < IMG_SIZE*NUM_SLICES; i++)
			img[i] = 1.0f;
		cout << "Done!" << endl;

		cout << "Starting ray_tracing..." << endl;

		N_x = 400;
		N_y = 240;
		delta_x = 400 / N_x;
		delta_y = 240 / N_y;


		unsigned __int64 cycle_temp = __rdtsc();
		double start = clock();
		//omp_set_nested(1);
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < DATA_SIZE; i++)
//		int i;
//		cilk_for (i = 0; i < DATA_SIZE; i++)
		{
			//std::cout << "Data slice " << data_ind << " out of " << DATA_SLICES << endl;
			//int img_offset = data_ind * IMG_SIZE;
			//int data_offset = data_ind * DATA_SIZE;
			//#pragma omp for schedule(dynamic)
			//for (int data_ind = 0; data_ind < NUM_SLICES; data_ind++)

			//int data_offset = data_ind*DATA_SLICES;
			//int img_offset = data_ind * 2 * IMG_SIZE;
			float src_x_t, src_y_t, det_x_t, det_y_t;
			int img_z_offset_t;
			float L;
			//float *l_total_array = new float[NUM_SLICES];
			float *l_total_array = (float *)_aligned_malloc(NUM_SLICES * sizeof(float), 32);
			for (int ind_temp = 0; ind_temp < NUM_SLICES; ind_temp++)
				l_total_array[ind_temp] = 0.0f;

			float lambda_x_0, lambda_x_Nx, lambda_y_0, lambda_y_Ny;
			float lambda_min, lambda_max;
			//cout << "i = " << i << endl;
			src_x_t = src_x[i];
			src_y_t = src_y[i];
			det_x_t = det_x[i];
			det_y_t = det_y[i];

			// compute the length
			L = sqrt((src_x_t - det_x_t)*(src_x_t - det_x_t) + (src_y_t - det_y_t)*(src_y_t - det_y_t));

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

				if (flag == 1)
				{
					while ((lambda_x - lambda_max) < -1e-2f || (lambda_y - lambda_max) < -1e-2f)
					{
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
								ind_x++;
							}

							else
							{
								// disp('Exit from y-plane');
								// Cikis y'den
								l = lambda_y - lambda_x;
								ind_y++;

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
								ind_y++;
							}
							else
							{
								// Cikis x'ten
								// disp('Exit from x-plane');
								l = lambda_x - lambda_y;
								ind_x++;

							}

							lambda_y = lambda_y_next;
							lambda_y_next = lambda_y + inc_y;
						}

						//l_total += l*1.0f;
						//cout << "x_old = " << ind_x_old << " y_old = " << ind_y_old << endl;
						//cout << "lambda_x = " << lambda_x << " next = " << lambda_x_next << " lambda_y = " << lambda_y << " next = " << lambda_y_next << endl;
						int img_offset = (ind_y_old - 1)*IMG_SIZE_2 + +(ind_x_old - 1)*NUM_SLICES;
						int data_offset = i*NUM_SLICES;

						__m256 l_temp = _mm256_set1_ps(l);
						for (int data_ind_2 = 0; data_ind_2 < NUM_SLICES; data_ind_2 += 8)
						{
							__m256 img_array = _mm256_load_ps(&img[data_ind_2 + img_offset]);
							__m256 l_total_array_temp = _mm256_load_ps(&l_total_array[data_ind_2]);
							l_total_array_temp = _mm256_add_ps(l_total_array_temp, _mm256_mul_ps(l_temp, img_array));
							_mm256_store_ps(&l_total_array[data_ind_2], l_total_array_temp);
						}
					}
				}
				else if (flag == 2)
				{
					while ((lambda_x - lambda_max) < -1e-2f || (lambda_y - lambda_max) < -1e-2f)
					{
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
								ind_x++;
							}

							else
							{
								// disp('Exit from y-plane');
								// Cikis y'den
								l = lambda_y - lambda_x;
								ind_y--;
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
								ind_y--;
							}
							else
							{
								// Cikis x'ten
								// disp('Exit from x-plane');
								l = lambda_x - lambda_y;
								ind_x++;
							}

							lambda_y = lambda_y_next;
							lambda_y_next = lambda_y + inc_y;
						}

						//l_total += l*1.0f;
						//cout << "x_old = " << ind_x_old << " y_old = " << ind_y_old << endl;
						//cout << "lambda_x = " << lambda_x << " next = " << lambda_x_next << " lambda_y = " << lambda_y << " next = " << lambda_y_next << endl;
						int img_offset = (ind_y_old - 1)*IMG_SIZE_2 + +(ind_x_old - 1)*NUM_SLICES;
						int data_offset = i*NUM_SLICES;

						__m256 l_temp = _mm256_set1_ps(l);
						for (int data_ind_2 = 0; data_ind_2 < NUM_SLICES; data_ind_2 += 8)
						{
							__m256 img_array = _mm256_load_ps(&img[data_ind_2 + img_offset]);
							__m256 l_total_array_temp = _mm256_load_ps(&l_total_array[data_ind_2]);
							l_total_array_temp = _mm256_add_ps(l_total_array_temp, _mm256_mul_ps(l_temp, img_array));
							_mm256_store_ps(&l_total_array[data_ind_2], l_total_array_temp);
						}
					}
				}
				else if (flag == 3)
				{
					while ((lambda_x - lambda_max) < -1e-2f || (lambda_y - lambda_max) < -1e-2f)
					{
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
								ind_x--;
							}

							else
							{
								// disp('Exit from y-plane');
								// Cikis y'den
								l = lambda_y - lambda_x;
								ind_y--;
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
								ind_y--;
							}
							else
							{
								// Cikis x'ten
								// disp('Exit from x-plane');
								l = lambda_x - lambda_y;
								ind_x--;
							}

							lambda_y = lambda_y_next;
							lambda_y_next = lambda_y + inc_y;
						}

						//l_total += l*1.0f;
						//cout << "x_old = " << ind_x_old << " y_old = " << ind_y_old << endl;
						//cout << "lambda_x = " << lambda_x << " next = " << lambda_x_next << " lambda_y = " << lambda_y << " next = " << lambda_y_next << endl;
						int img_offset = (ind_y_old - 1)*IMG_SIZE_2 +(ind_x_old - 1)*NUM_SLICES;
						int data_offset = i*NUM_SLICES;

						__m256 l_temp = _mm256_set1_ps(l);
						for (int data_ind_2 = 0; data_ind_2 < NUM_SLICES; data_ind_2 += 8)
						{
							//cout << "data_ind_2 = " << data_ind_2 << " img_offset = " << img_offset << endl;
							__m256 img_array = _mm256_load_ps(&img[data_ind_2 + img_offset]);
							__m256 l_total_array_temp = _mm256_load_ps(&l_total_array[data_ind_2]);
							l_total_array_temp = _mm256_add_ps(l_total_array_temp, _mm256_mul_ps(l_temp, img_array));
							_mm256_store_ps(&l_total_array[data_ind_2], l_total_array_temp);
						}

						/*
						int img_offset = (ind_y_old - 1)*IMG_SIZE_2 + +(ind_x_old - 1)*NUM_SLICES;
						int data_offset = i*NUM_SLICES;

						for (int data_ind = 0; data_ind < NUM_SLICES; data_ind++)
						{
							// THIS IS V4P2, OLD MEMORY SEQUENCE
							//int img_offset = data_ind * IMG_SIZE;
							//int data_offset = data_ind * DATA_SIZE;
							//l_total_array[data_ind] += l * img[(ind_x_old - 1) + (ind_y_old - 1) * 400 + img_offset];
							// How fast the indices change: x > y > z
							// ind_x + ind_y*size_x + ind_z*size_x*size_y

							// THIS IS V4P3, NEW PROPOSED ONE
							//int img_offset = (ind_y_old - 1)*IMG_SIZE_2;
							//int data_offset = i*NUM_SLICES;
							//l_total_array[data_ind] += l * img[(data_ind) + (ind_x_old - 1)*NUM_SLICES + img_offset];
							// How fast the indices change: z > x > y
							// ind_z + ind_x*size_z + ind_y*size_z*size_y
							l_total_array[data_ind] += l * img[data_ind + img_offset];

						}
						*/

					}
				}
				else //flag = 4
				{
					while ((lambda_x - lambda_max) < -1e-2f || (lambda_y - lambda_max) < -1e-2f)
					{
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
								ind_x--;
							}

							else
							{
								// disp('Exit from y-plane');
								// Cikis y'den
								l = lambda_y - lambda_x;
								ind_y++;
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
								ind_y++;
							}
							else
							{
								// Cikis x'ten
								// disp('Exit from x-plane');
								l = lambda_x - lambda_y;
								ind_x--;
							}

							lambda_y = lambda_y_next;
							lambda_y_next = lambda_y + inc_y;
						}

						//l_total += l*1.0f;
						//cout << "x_old = " << ind_x_old << " y_old = " << ind_y_old << endl;
						//cout << "lambda_x = " << lambda_x << " next = " << lambda_x_next << " lambda_y = " << lambda_y << " next = " << lambda_y_next << endl;
						int img_offset = (ind_y_old - 1)*IMG_SIZE_2 + +(ind_x_old - 1)*NUM_SLICES;
						int data_offset = i*NUM_SLICES;

						__m256 l_temp = _mm256_set1_ps(l);
						for (int data_ind_2 = 0; data_ind_2 < NUM_SLICES; data_ind_2 += 8)
						{
							__m256 img_array = _mm256_load_ps(&img[data_ind_2 + img_offset]);
							__m256 l_total_array_temp = _mm256_load_ps(&l_total_array[data_ind_2]);
							l_total_array_temp = _mm256_add_ps(l_total_array_temp, _mm256_mul_ps(l_temp, img_array));
							_mm256_store_ps(&l_total_array[data_ind_2], l_total_array_temp);
						}
					}
				}
				for (int data_ind = 0; data_ind < NUM_SLICES; data_ind++)
				{
					int data_offset = data_ind * DATA_SIZE;
					fp_ones[i + data_offset] = l_total_array[data_ind];

				}

				//cout << "Forward projection onto index " << i << " = " << l_total << endl;
				// system("pause");
			}
		}
		double end = clock();
		clock_temp = __rdtsc() - cycle_temp;
		printf_s("%I64d ticks\n", clock_temp);
		time_temp = (end - start) / CLOCKS_PER_SEC;
		std::cout << "Trial " << trial_ctr << " Time = " << time_temp << endl;
		float fp_ones_max = 0.0f;

		for (int i = 0; i<126464; i++)
		{
			if (fp_ones[i]>fp_ones_max)
				fp_ones_max = fp_ones[i];
		}
		std::cout << "fp_ones_max = " << fp_ones_max << endl;
		time_total += time_temp;
		//system("pause");
	}
	avg_time = time_total / trial_max;
	std::cout << "Ray tracing cpp | SureScan | v4p4 - modified for loop / modified memory structure / AVX usage | forward projection of ones | Average time = " << avg_time << endl;
	//std::system("pause");
	return 0;
}