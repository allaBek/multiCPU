#include <iostream>
#include "lodepng.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>  // for string streams
#include "functions.h"
#include <thread>
#include <omp.h>
using namespace std;

void ZNCC_parallel_openMP(vector< vector <unsigned char> >&sample1, vector< vector <unsigned char> >&sample2, vector<vector<zncc_parameters>> &zncc1, vector<vector<zncc_parameters>> &zncc2, vector<vector<zncc_parameters>> &zncc3)
{
	cout << "zncc" << endl;
	//declare ZNCC algorithm parameters
	const int HEIGHT = sample1.size(), WIDTH = sample1[0].size(), MAX_DISP = 65, B = 9, DISP_DIFF = 60;
	double mean1, mean2;
	//create template vectors
	for (int y = 0; y < HEIGHT; y++)
	{
		vector<zncc_parameters> zncc1_temp, zncc2_temp, zncc3_temp;
		for (int x = 0; x < WIDTH; x++)
		{
			//zncc for image 2 on image one
			zncc1_temp.push_back(zncc_parameters());
			zncc1_temp.back().x = x;
			zncc1_temp.back().y = y;
			//zncc for image 1 on image 2
			zncc1_temp.back().zncc = -10;
			zncc2_temp.push_back(zncc_parameters());
			zncc2_temp.back().x = x;
			zncc2_temp.back().y = y;
			zncc2_temp.back().zncc = -10;
			//disparity map values
			zncc3_temp.push_back(zncc_parameters());
			zncc3_temp.back().x = x;
			zncc3_temp.back().y = y;
			zncc3_temp.back().zncc = -10;
		}
		zncc1.push_back(zncc1_temp);
		zncc2.push_back(zncc2_temp);
		zncc3.push_back(zncc3_temp);
	}

#pragma omp parallel for  default(none) shared(HEIGHT, WIDTH, MAX_DISP, B, DISP_DIFF, sample1, sample2, zncc1, zncc2, zncc3)
	for (int y = 0; y < HEIGHT; y++)
	{
		//	#pragma omp parallel for  default(none) shared(HEIGHT, WIDTH, MAX_DISP, B, DISP_DIFF, sample1, sample2, zncc1, zncc2, zncc3)
		for (int x = 0; x < WIDTH; x++)
		{
			//calculate the mean
			getMean(sample1, sample2, x, y, B, mean1, mean2, HEIGHT, WIDTH);
			for (int d = 0; d < MAX_DISP; d++)
			{
				//calculate the ZNCC
				const int H = ((B - 1) / 2);
				double zncc1_sum1 = 0, zncc1_sum2 = 0, zncc1_sum3 = 0;
				double zncc2_sum1 = 0, zncc2_sum2 = 0, zncc2_sum3 = 0;
				for (int y_w = y - H; y_w <= (y + H) && (y + H) < WIDTH; y_w++)
				{
					for (int x_w = x - H; x_w <= x + H; x_w++)
					{
						if (x_w > -1 && y_w > -1 && y_w < HEIGHT)
						{
							//cout << (x_w) << "		" << d << "		" << (x_w - d)<< endl;
							if ((x_w - d) > -1 && (x_w - d) < WIDTH)
							{
								//zncc sum equations
								//img2 on img1
								unsigned char s1 = sample1[y_w][x_w], s2 = sample2[y_w][x_w - d];
								zncc1_sum1 = zncc1_sum1 + (s1 - mean1)*(s2 - mean2);
								zncc1_sum2 = zncc1_sum2 + (s1 - mean1)*(s1 - mean1);
								zncc1_sum3 = zncc1_sum3 + (s2 - mean2)*(s2 - mean2);

							}
							if (x_w + d < WIDTH)
							{
								unsigned char s1, s2;
								s2 = sample2[y_w][x_w]; s1 = sample1[y_w][x_w + d];
								//img1 on img2
								zncc2_sum1 = zncc2_sum1 + (s2 - mean2)*(s1 - mean1);
								zncc2_sum2 = zncc2_sum2 + (s2 - mean2) *(s2 - mean2);
								zncc2_sum3 = zncc2_sum3 + (s1 - mean1) *(s1 - mean1);
							}

						}
					}
				}
				//finding final zncc value for case 1
				zncc1_sum2 = sqrt(zncc1_sum2);
				zncc1_sum3 = sqrt(zncc1_sum3);
				double temp = zncc1_sum1 / (zncc1_sum2*zncc1_sum3);
				if (zncc1[y][x].zncc <= temp)
				{
					zncc1[y][x].zncc = temp;
					zncc1[y][x].disparity = (unsigned char)(d * 255 / MAX_DISP);
				}
				//finding final zncc value for case 2
				zncc2_sum2 = sqrt(zncc2_sum2);
				zncc2_sum3 = sqrt(zncc2_sum3);
				temp = zncc2_sum1 / (zncc2_sum2*zncc2_sum3);
				if (zncc2[y][x].zncc <= temp)
				{
					zncc2[y][x].zncc = temp;
					zncc2[y][x].disparity = (unsigned char)(d * 255 / 65);
				}



			}
			if (abs(zncc1[y][x].disparity - zncc2[y][x].disparity) < DISP_DIFF)
				zncc3[y][x].disparity = zncc1[y][x].disparity > zncc2[y][x].disparity ? zncc2[y][x].disparity : zncc1[y][x].disparity;
			else
				zncc3[y][x].disparity = 0;

		}
	}

}

void ZNCC_parallel(vector< vector <unsigned char> >&sample1, vector< vector <unsigned char> >&sample2, vector<vector<zncc_parameters>> &zncc1, vector<vector<zncc_parameters>> &zncc2, vector<vector<zncc_parameters>> &zncc3)
{
	cout << "zncc" << endl;
	//declare ZNCC algorithm parameters
	const int HEIGHT = sample1.size(), WIDTH = sample1[0].size(), MAX_DISP = 65, B = 9, DISP_DIFF = 50;
	double mean1, mean2;
	//this loop is used to create a template vector, to just place later values on it
	for (int y = 0; y < HEIGHT; y++)
	{
		vector<zncc_parameters> zncc1_temp, zncc2_temp, zncc3_temp;
		for (int x = 0; x < WIDTH; x++)
		{
			//zncc for image 2 on image one
			zncc1_temp.push_back(zncc_parameters());
			zncc1_temp.back().x = x;
			zncc1_temp.back().y = y;
			//zncc for image 1 on image 2
			zncc2_temp.push_back(zncc_parameters());
			zncc2_temp.back().x = x;
			zncc2_temp.back().y = y;
			//disparity map values
			zncc3_temp.push_back(zncc_parameters());
			zncc3_temp.back().x = x;
			zncc3_temp.back().y = y;
		}
		zncc1.push_back(zncc1_temp);
		zncc2.push_back(zncc2_temp);
		zncc3.push_back(zncc3_temp);
	}
	//uncomment if you want to make 504 threads, each line in a separate thread
	/*
	vector<thread> threads;

	for (int y = 0; y < HEIGHT; y++)
	{
		threads.push_back(thread(processing_per_line, ref(sample1), ref(sample2), ref(zncc1[y]), ref(zncc2[y]), ref(zncc3[y]), y));

	}
	for (auto& th : threads)
		th.join();
	*/
	//uncomment of you want to make H/M threads

	vector<thread> threads;
	int M = 252;
	for (int y = 0; y < HEIGHT / M; y++)
	{
		for (int i = 0; i < M; i++)
			threads.push_back(thread(processing_per_line, ref(sample1), ref(sample2), ref(zncc1[int(y + i * HEIGHT / M)]), ref(zncc2[int(y + i * HEIGHT / M)]), ref(zncc3[int(y + i * HEIGHT / M)]), int(y + i * HEIGHT / M)));
	}
	for (auto& th : threads)
		th.join();


	//uncomment if you want to make a single thread looping aroundall lines
	/*
	for (int y = 0; y < HEIGHT; y++)
	{
		processing_per_line(sample1,sample2,zncc1[y],zncc2[y],zncc3[y], y);

	}
	*/
}
void processing_per_line(vector< vector <unsigned char> >&sample1, vector< vector <unsigned char> >&sample2, vector<zncc_parameters> &zncc1, vector<zncc_parameters> &zncc2, vector<zncc_parameters> &zncc3, int y)
{

	const int HEIGHT = sample1.size(), WIDTH = sample1[0].size(), MAX_DISP = 65, B = 9, DISP_DIFF = 50;
	//uncomment this when using thread per pixel
	/*
	vector<thread> threads;
	for (int x = 0; x < WIDTH; x++)
	{
		threads.push_back(thread(processing_per_pixel, ref(sample1), ref(sample2), ref(zncc1[x]), ref(zncc2[x]), ref(zncc3[x]), x, y));
	}
	for (auto& th : threads)
		th.join();
	*/
	//uncomment if you want to make W/M threads
	/*
	vector<thread> threads;
	int M = 1;
	for (int x = 0; x < WIDTH/M; x++)
	{
		for(int i = 0; i < M; i++)
			threads.push_back(thread(processing_per_pixel, ref(sample1), ref(sample2), ref(zncc1[x + i*WIDTH/M]), ref(zncc2[x + i * WIDTH / M]), ref(zncc3[x + i * WIDTH / M]), x + i * WIDTH / M, y));
	}
	for (auto& th : threads)
		th.join();
	*/
	//uncomment this when using one thread per line

	for (int x = 0; x < WIDTH; x++)
	{
		processing_per_pixel(sample1, sample2, zncc1[x], zncc2[x], zncc3[x], x, y);
	}


}

void processing_per_pixel(vector< vector <unsigned char> >&sample1, vector< vector <unsigned char> >&sample2, zncc_parameters &zncc1, zncc_parameters &zncc2, zncc_parameters &zncc3, int x, int y)
{

	const int HEIGHT = sample1.size(), WIDTH = sample1[0].size(), MAX_DISP = 65, B = 9, DISP_DIFF = 50;
	double mean1, mean2;
	//zncc for image 2 on image one
	zncc1.x = x;
	zncc1.y = y;
	//zncc for image 1 on image 2
	zncc2.x = x;
	zncc2.y = y;
	//disparity map values
	zncc3.x = x;
	zncc3.y = y;
	//calculate the mean
	getMean(sample1, sample2, x, y, B, mean1, mean2, HEIGHT, WIDTH);
	for (int d = 0; d < MAX_DISP; d++)
	{
		//calculate the ZNCC
		const int H = ((B - 1) / 2);
		double zncc1_sum1 = 0, zncc1_sum2 = 0, zncc1_sum3 = 0;
		double zncc2_sum1 = 0, zncc2_sum2 = 0, zncc2_sum3 = 0;
		for (int y_w = y - H; y_w <= (y + H) && (y + H) < WIDTH; y_w++)
		{
			for (int x_w = x - H; x_w <= x + H; x_w++)
			{
				if (x_w > -1 && y_w > -1 && y_w < HEIGHT)
				{
					//cout << (x_w) << "		" << d << "		" << (x_w - d)<< endl;
					if ((x_w - d) > -1 && (x_w - d) < WIDTH)
					{
						//zncc sum equations
						//img2 on img1
						unsigned char s1 = sample1[y_w][x_w], s2 = sample2[y_w][x_w - d];
						zncc1_sum1 = zncc1_sum1 + (s1 - mean1)*(s2 - mean2);
						zncc1_sum2 = zncc1_sum2 + (s1 - mean1) *(s1 - mean1);
						zncc1_sum3 = zncc1_sum3 + (s2 - mean2) *(s2 - mean2);

					}
					if (x_w + d < WIDTH)
					{
						unsigned char s1, s2;
						s2 = sample2[y_w][x_w]; s1 = sample1[y_w][x_w + d];
						//img1 on img2
						zncc2_sum1 = zncc2_sum1 + (s2 - mean2)*(s1 - mean1);
						zncc2_sum2 = zncc2_sum2 + (s2 - mean2) *(s2 - mean2);
						zncc2_sum3 = zncc2_sum3 + (s1 - mean1) *(s1 - mean1);
					}

				}
			}
		}
		//finding final zncc value for case 1
		zncc1_sum2 = sqrt(zncc1_sum2);
		zncc1_sum3 = sqrt(zncc1_sum3);
		double temp = zncc1_sum1 / (zncc1_sum2*zncc1_sum3);
		if (zncc1.zncc <= temp)
		{
			zncc1.zncc = temp;
			zncc1.disparity = (unsigned char)(d * 255 / MAX_DISP);
		}
		//finding final zncc value for case 2
		zncc2_sum2 = sqrt(zncc2_sum2);
		zncc2_sum3 = sqrt(zncc2_sum3);
		temp = zncc2_sum1 / (zncc2_sum2*zncc2_sum3);
		if (zncc2.zncc <= temp)
		{
			zncc2.zncc = temp;
			zncc2.disparity = (unsigned char)(d * 255 / 65);
		}
	}
	if (abs(zncc1.disparity - zncc2.disparity) < DISP_DIFF)
		zncc3.disparity = zncc1.disparity > zncc2.disparity ? zncc2.disparity : zncc1.disparity;
	else
		zncc3.disparity = 0;

}