#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <thread>
#include <omp.h>
#include "lodepng.h"
#include "functions.h"

using namespace std;

/*
	This function basically calculates the ZNCC values of every raw gray-scaled image samples using an openMP approach.

	Parameters
	- vector<vector<unsigned char>> &sample1 , The first gray-scaled image sample.
	- vector<vector<unsigned char>> &sample2, The second gray-scaled image sample.
	- vector<vector<zncc_parameters>> &zncc1 , Data structure to hold ZNCC values for image 2 on image 1.
	- vector<vector<zncc_parameters>> &zncc2 , Data structure to hold ZNCC values for image 1 on image 2.
	- vector<vector<zncc_parameters>> &zncc3 , Data structure to hold disparity map values.

	Return
	- This function does not return anything.
*/
void ZNCC_parallel_openMP(vector<vector<unsigned char>> &sample1, vector<vector<unsigned char>> &sample2, 
	vector<vector<zncc_parameters>> &zncc1, vector<vector<zncc_parameters>> &zncc2, vector<vector<zncc_parameters>> &zncc3)
{
	cout << "ZNCC parallel evaluation stage has begun." << endl;

	// Declaring the variables that will be used in ZNCC algorithm.
	const int HEIGHT = sample1.size(), WIDTH = sample1[0].size(), MAX_DISP = 65, B = 9, DISP_DIFF = 60;
	double mean1, mean2;

	// Starting the algorithm.
	for (int y = 0; y < HEIGHT; y++)
	{
		// Temporary data structure to store calculations of a row of image data. 
		vector<zncc_parameters> zncc1_temp, zncc2_temp, zncc3_temp;

		// Creating template objects.
		for (int x = 0; x < WIDTH; x++)
		{
			// New ZNCC parameters object with the pixel (x,y) for image 2 on image 1
			zncc1_temp.push_back(zncc_parameters());
			zncc1_temp.back().x = x;
			zncc1_temp.back().y = y;
			zncc1_temp.back().zncc = -10;

			// New ZNCC parameters object with the pixel (x,y) for image 2 on image 1
			zncc2_temp.push_back(zncc_parameters());
			zncc2_temp.back().x = x;
			zncc2_temp.back().y = y;
			zncc2_temp.back().zncc = -10;

			// New ZNCC parameters with the pixel (x,y) for disparity values.
			zncc3_temp.push_back(zncc_parameters());
			zncc3_temp.back().x = x;
			zncc3_temp.back().y = y;
			zncc3_temp.back().zncc = -10;
		}

		// Pushing values into the storing data structures.
		zncc1.push_back(zncc1_temp);
		zncc2.push_back(zncc2_temp);
		zncc3.push_back(zncc3_temp);
	}

	// From now on, parallel processing starts...

	#pragma omp parallel for default(none) shared(HEIGHT, WIDTH, MAX_DISP, B, DISP_DIFF, sample1, sample2, zncc1, zncc2, zncc3)
	for (int y = 0; y < HEIGHT; y++)
	{
		/*	
		If you wanna make width also parallel, uncomment this section.
		#pragma omp parallel for  default(none) shared(HEIGHT, WIDTH, MAX_DISP, B, DISP_DIFF, sample1, sample2, zncc1, zncc2, zncc3) 
		*/
		for (int x = 0; x < WIDTH; x++)
		{
			// Calculating the mean for both images for pixel (x,y) using window size B.
			getMean(sample1, sample2, x, y, B, mean1, mean2, HEIGHT, WIDTH);

			// Calculate the ZNCC values.
			for (int d = 0; d < MAX_DISP; d++)
			{
				// Initialize the necessary variables.
				const int H = ((B - 1) / 2);
				double zncc1_sum1 = 0, zncc1_sum2 = 0, zncc1_sum3 = 0;
				double zncc2_sum1 = 0, zncc2_sum2 = 0, zncc2_sum3 = 0;

				// Only taking the pixels inside the window frame sized B and in the picture.
				for (int y_w = y - H; y_w <= (y + H) && (y + H) < WIDTH; y_w++)
				{
					for (int x_w = x - H; x_w <= x + H; x_w++)
					{
						if (x_w > -1 && y_w > -1 && y_w < HEIGHT)
						{
							// ZNCC sum equations for image 2 on image 1.
							if ((x_w - d) > -1 && (x_w - d) < WIDTH)
							{
								unsigned char s1 = sample1[y_w][x_w], s2 = sample2[y_w][x_w - d];

								zncc1_sum1 = zncc1_sum1 + (s1 - mean1)*(s2 - mean2);
								zncc1_sum2 = zncc1_sum2 + (s1 - mean1) *(s1 - mean1);
								zncc1_sum3 = zncc1_sum3 + (s2 - mean2) *(s2 - mean2);
							}

							// ZNCC sum equations for image 1 on image 2.
							if (x_w + d < WIDTH)
							{
								unsigned char s1 = sample1[y_w][x_w + d], s2 = sample2[y_w][x_w];

								zncc2_sum1 = zncc2_sum1 + (s2 - mean2)*(s1 - mean1);
								zncc2_sum2 = zncc2_sum2 + (s2 - mean2) *(s2 - mean2);
								zncc2_sum3 = zncc2_sum3 + (s1 - mean1) *(s1 - mean1);
							}
						}
					}
				}

				// Finding the final ZNCC value for the first case by evaluating the values described in the algorithm.
				zncc1_sum2 = sqrt(zncc1_sum2);
				zncc1_sum3 = sqrt(zncc1_sum3);
				double final1 = zncc1_sum1 / (zncc1_sum2 * zncc1_sum3);

				// Finding the final ZNCC value for the second case by evaluating the values described in the algorithm.
				zncc2_sum2 = sqrt(zncc2_sum2);
				zncc2_sum3 = sqrt(zncc2_sum3);
				double final2 = zncc2_sum1 / (zncc2_sum2 * zncc2_sum3);

				// Executing the algorithm's last step.
				if (zncc1[y][x].zncc <= final1)
				{
					zncc1[y][x].zncc = final1;
					zncc1[y][x].disparity = (unsigned char)(d * 255 / MAX_DISP);
				}

				if (zncc2[y][x].zncc <= final2)
				{
					zncc2[y][x].zncc = final2;
					zncc2[y][x].disparity = (unsigned char)(d * 255 / 65);
				}
			}

			// If the disparity difference is inside the threshold, taking it into the disparity map.
			if (abs(zncc1[y][x].disparity - zncc2[y][x].disparity) < DISP_DIFF)
				zncc3[y][x].disparity = zncc1[y][x].disparity > zncc2[y][x].disparity ? zncc2[y][x].disparity : zncc1[y][x].disparity;
			else
				zncc3[y][x].disparity = 0;
		}
	}
}

/*
	This function basically helps to calculate the ZNCC values of every raw gray-scaled image samples using threads approach as for the multiprogramming.

	Parameters
	- vector<vector<unsigned char>> &sample1 , The first gray-scaled image sample.
	- vector<vector<unsigned char>> &sample2, The second gray-scaled image sample.
	- vector<vector<zncc_parameters>> &zncc1 , Data structure to hold ZNCC values for image 2 on image 1.
	- vector<vector<zncc_parameters>> &zncc2 , Data structure to hold ZNCC values for image 1 on image 2.
	- vector<vector<zncc_parameters>> &zncc3 , Data structure to hold disparity map values.

	Return
	- This function does not return anything.
*/
void ZNCC_parallel(vector<vector<unsigned char>> &sample1, vector<vector<unsigned char>> &sample2, 
	vector<vector<zncc_parameters>> &zncc1, vector<vector<zncc_parameters>> &zncc2, vector<vector<zncc_parameters>> &zncc3)
{
	cout << "ZNCC parallel evaluation stage has begun." << endl;

	// Declaring the variables that will be used in ZNCC algorithm.
	const int HEIGHT = sample1.size(), WIDTH = sample1[0].size(), MAX_DISP = 65, B = 9, DISP_DIFF = 60;
	double mean1, mean2;

	// Starting the algorithm.
	for (int y = 0; y < HEIGHT; y++)
	{
		// Temporary data structure to store calculations of a row of image data. 
		vector<zncc_parameters> zncc1_temp, zncc2_temp, zncc3_temp;

		// Creating template objects.
		for (int x = 0; x < WIDTH; x++)
		{
			// New ZNCC parameters object with the pixel (x,y) for image 2 on image 1
			zncc1_temp.push_back(zncc_parameters());
			zncc1_temp.back().x = x;
			zncc1_temp.back().y = y;

			// New ZNCC parameters object with the pixel (x,y) for image 2 on image 1
			zncc2_temp.push_back(zncc_parameters());
			zncc2_temp.back().x = x;
			zncc2_temp.back().y = y;

			// New ZNCC parameters with the pixel (x,y) for disparity values.
			zncc3_temp.push_back(zncc_parameters());
			zncc3_temp.back().x = x;
			zncc3_temp.back().y = y;
		}

		// Pushing values into the storing data structures.
		zncc1.push_back(zncc1_temp);
		zncc2.push_back(zncc2_temp);
		zncc3.push_back(zncc3_temp);
	}

	/*
	If you wanna make 504 threads and run them as parallel, uncomment this section.
	That makes each row in a seperate thread.

	vector<thread> threads;
	for (int y = 0; y < HEIGHT; y++)
	{
		threads.push_back(thread(processing_per_line, ref(sample1), ref(sample2), ref(zncc1[y]), ref(zncc2[y]), ref(zncc3[y]), y));
	}

	for (auto& th : threads)
		th.join();
	*/

	/*
	If you wanna make HEIGHT / M (M is an integer) threads and run them as parallel, uncomment this section.

	vector<thread> threads;
	int M = 252;

	for (int y = 0; y < HEIGHT / M; y++)
	{
		for (int i = 0; i < M; i++)
			threads.push_back(thread(processing_per_line, ref(sample1), ref(sample2), ref(zncc1[int(y + i * HEIGHT / M)]), 
			ref(zncc2[int(y + i * HEIGHT / M)]), ref(zncc3[int(y + i * HEIGHT / M)]), int(y + i * HEIGHT / M)));
	}

	for (auto& th : threads)
		th.join();
	*/

	// Using one single thread looping around all lines
	for (int y = 0; y < HEIGHT; y++)
	{
		processing_per_line(sample1, sample2, zncc1[y], zncc2[y], zncc3[y], y);
	}
}

/*
	This function basically helps to calculate the ZNCC values of every raw gray-scaled image samples using threads approach as for the multiprogramming.
	It enables multithreading as line-based.

	Parameters
	- vector<vector<unsigned char>> &sample1 , The first gray-scaled image sample.
	- vector<vector<unsigned char>> &sample2, The second gray-scaled image sample.
	- vector<vector<zncc_parameters>> &zncc1 , Data structure to hold ZNCC values for image 2 on image 1.
	- vector<vector<zncc_parameters>> &zncc2 , Data structure to hold ZNCC values for image 1 on image 2.
	- vector<vector<zncc_parameters>> &zncc3 , Data structure to hold disparity map values.
	- int y , The row number in the image.

	Return
	- This function does not return anything.
*/
void processing_per_line(vector<vector<unsigned char>> &sample1, vector<vector<unsigned char>> &sample2, 
	vector<zncc_parameters> &zncc1, vector<zncc_parameters> &zncc2, vector<zncc_parameters> &zncc3, int y)
{
	// Declaring the variables that will be used in ZNCC algorithm.
	const int HEIGHT = sample1.size(), WIDTH = sample1[0].size(), MAX_DISP = 65, B = 9, DISP_DIFF = 50;

	/*
	If you want to use thread per pixel, uncomment this section.

	vector<thread> threads;
	for (int x = 0; x < WIDTH; x++)
	{
		threads.push_back(thread(processing_per_pixel, ref(sample1), ref(sample2), ref(zncc1[x]), ref(zncc2[x]), ref(zncc3[x]), x, y));
	}

	for (auto& th : threads)
		th.join();
	*/

	/*
	If you wanna make WIDTH / M (M is an integer) threads and run them as parallel, uncomment this section.

	vector<thread> threads;
	int M = 1;

	for (int x = 0; x < WIDTH/M; x++)
	{
		for(int i = 0; i < M; i++)
			threads.push_back(thread(processing_per_pixel, ref(sample1), ref(sample2), 
			ref(zncc1[x + i*WIDTH/M]), ref(zncc2[x + i * WIDTH / M]), ref(zncc3[x + i * WIDTH / M]), x + i * WIDTH / M, y));
	}

	for (auto& th : threads)
		th.join();
	*/

	// Using one thread per line...
	for (int x = 0; x < WIDTH; x++)
	{
		processing_per_pixel(sample1, sample2, zncc1[x], zncc2[x], zncc3[x], x, y);
	}
}

/*
	This function basically helps to calculate the ZNCC values of every raw gray-scaled image samples using threads approach as for the multiprogramming.
	It enables multithreading as pixel-based.

	Parameters
	- vector<vector<unsigned char>> &sample1 , The first gray-scaled image sample.
	- vector<vector<unsigned char>> &sample2, The second gray-scaled image sample.
	- vector<vector<zncc_parameters>> &zncc1 , Data structure to hold ZNCC values for image 2 on image 1.
	- vector<vector<zncc_parameters>> &zncc2 , Data structure to hold ZNCC values for image 1 on image 2.
	- vector<vector<zncc_parameters>> &zncc3 , Data structure to hold disparity map values.
	- int y , The row number in the image for pixel (x,y).
	- int y , The column number in the image for pixel (x,y).

	Return
	- This function does not return anything.
*/
void processing_per_pixel(vector<vector<unsigned char>> &sample1, vector<vector<unsigned char>> &sample2, zncc_parameters &zncc1, zncc_parameters &zncc2, zncc_parameters &zncc3, int x, int y)
{
	// Declaring the variables that will be used in ZNCC algorithm.
	const int HEIGHT = sample1.size(), WIDTH = sample1[0].size(), MAX_DISP = 65, B = 9, DISP_DIFF = 50;
	double mean1, mean2;

	// New ZNCC parameters object with the pixel (x,y) for image 2 on image 1
	zncc1.x = x;
	zncc1.y = y;

	// New ZNCC parameters object with the pixel (x,y) for image 1 on image 2
	zncc2.x = x;
	zncc2.y = y;

	// New ZNCC parameters with the pixel (x,y) for disparity values.
	zncc3.x = x;
	zncc3.y = y;

	// Calculating the mean for both images for pixel (x,y) using window size B.
	getMean(sample1, sample2, x, y, B, mean1, mean2, HEIGHT, WIDTH);

	// Calculate the ZNCC values.
	for (int d = 0; d < MAX_DISP; d++)
	{
		// Initialize the necessary variables.
		const int H = ((B - 1) / 2);
		double zncc1_sum1 = 0, zncc1_sum2 = 0, zncc1_sum3 = 0;
		double zncc2_sum1 = 0, zncc2_sum2 = 0, zncc2_sum3 = 0;

		// Only taking the pixels inside the window frame sized B and in the picture.
		for (int y_w = y - H; y_w <= (y + H) && (y + H) < WIDTH; y_w++)
		{
			for (int x_w = x - H; x_w <= x + H; x_w++)
			{
				if (x_w > -1 && y_w > -1 && y_w < HEIGHT)
				{
					// ZNCC sum equations for image 2 on image 1.
					if ((x_w - d) > -1 && (x_w - d) < WIDTH)
					{
						unsigned char s1 = sample1[y_w][x_w], s2 = sample2[y_w][x_w - d];

						zncc1_sum1 = zncc1_sum1 + (s1 - mean1)*(s2 - mean2);
						zncc1_sum2 = zncc1_sum2 + (s1 - mean1) *(s1 - mean1);
						zncc1_sum3 = zncc1_sum3 + (s2 - mean2) *(s2 - mean2);
					}

					// ZNCC sum equations for image 1 on image 2.
					if (x_w + d < WIDTH)
					{
						unsigned char s1 = sample1[y_w][x_w + d], s2 = sample2[y_w][x_w];

						zncc2_sum1 = zncc2_sum1 + (s2 - mean2)*(s1 - mean1);
						zncc2_sum2 = zncc2_sum2 + (s2 - mean2) *(s2 - mean2);
						zncc2_sum3 = zncc2_sum3 + (s1 - mean1) *(s1 - mean1);
					}
				}
			}
		}

		// Finding the final ZNCC value for the first case by evaluating the values described in the algorithm.
		zncc1_sum2 = sqrt(zncc1_sum2);
		zncc1_sum3 = sqrt(zncc1_sum3);
		double temp = zncc1_sum1 / (zncc1_sum2*zncc1_sum3);

		if (zncc1.zncc <= temp)
		{
			zncc1.zncc = temp;
			zncc1.disparity = (unsigned char)(d * 255 / MAX_DISP);
		}

		// Finding the final ZNCC value for the second case by evaluating the values described in the algorithm.
		zncc2_sum2 = sqrt(zncc2_sum2);
		zncc2_sum3 = sqrt(zncc2_sum3);
		temp = zncc2_sum1 / (zncc2_sum2*zncc2_sum3);

		if (zncc2.zncc <= temp)
		{
			zncc2.zncc = temp;
			zncc2.disparity = (unsigned char)(d * 255 / 65);
		}
	}

	// If the disparity difference is inside the threshold, taking it into the disparity map.
	if (abs(zncc1.disparity - zncc2.disparity) < DISP_DIFF)
		zncc3.disparity = zncc1.disparity > zncc2.disparity ? zncc2.disparity : zncc1.disparity;
	else
		zncc3.disparity = 0;
}