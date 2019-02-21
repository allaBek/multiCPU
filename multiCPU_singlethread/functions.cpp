#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include "functions.h"
#include "lodepng.h"

using namespace std;

/*
	This function basically calculates the ZNCC values of every raw gray-scaled image samples.

	Parameters
	- vector<vector<unsigned char>> &sample1 , The first gray-scaled image sample.
	- vector<vector<unsigned char>> &sample2, The second gray-scaled image sample.
	- vector<vector<zncc_parameters>> &zncc1 , Data structure to hold ZNCC values for image 2 on image 1.
	- vector<vector<zncc_parameters>> &zncc2 , Data structure to hold ZNCC values for image 1 on image 2.
	- vector<vector<zncc_parameters>> &zncc3 , Data structure to hold disparity map values.

	Return
	- This function does not return anything.
*/
void ZNCC(vector<vector<unsigned char>> &sample1, vector<vector<unsigned char>> &sample2, 
	vector<vector<zncc_parameters>> &zncc1, vector<vector<zncc_parameters>> &zncc2, vector<vector<zncc_parameters>> &zncc3)
{
	cout << "ZNCC evaluation stage has begun." << endl;

	// Declaring the variables that will be used in ZNCC algorithm.
	const int HEIGHT = sample1.size(), WIDTH = sample1[0].size(), MAX_DISP = 65, B = 9, DISP_DIFF = 60;
	double mean1, mean2;

	// Starting the algorithm.
	for (int y = 0; y < HEIGHT; y++)
	{
		// Temporary data structure to store calculations of a row of image data. 
		vector<zncc_parameters> zncc1_temp, zncc2_temp, zncc3_temp;

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
				for (int y_w = y - H; y_w <= (y + H) && (y+H) < WIDTH; y_w++)
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
				if (zncc1_temp.back().zncc <= final1)
				{
					zncc1_temp.back().zncc = final1;
					zncc1_temp.back().disparity = (unsigned char) (d * 255 / MAX_DISP);
				}

				if (zncc2_temp.back().zncc <= final2)
				{
					zncc2_temp.back().zncc = final2;
					zncc2_temp.back().disparity = (unsigned char) (d * 255 / 65);
				}
			}

			// If the disparity difference is inside the threshold, taking it into the disparity map.
			if (abs(zncc1_temp.back().disparity - zncc2_temp.back().disparity) < DISP_DIFF)
				zncc3_temp.back().disparity = zncc1_temp.back().disparity>zncc2_temp.back().disparity ? zncc2_temp.back().disparity : zncc1_temp.back().disparity;
			else
				zncc3_temp.back().disparity = 0;
		}

		// Pushing values into the storing data structures.
		zncc1.push_back(zncc1_temp);
		zncc2.push_back(zncc2_temp);
		zncc3.push_back(zncc3_temp);
	}
}

/*
	This function is basically used in calculating mean values in a defined window frame of pixels.

	Parameters
	- vector<vector<unsigned char>>&img1 , The first gray-scaled image sample.
	- vector<vector<unsigned char>> &img2, The second gray-scaled image sample.
	- unsigned int x_bias , x-axis pixel coordinate.
	- unsigned int y_bias , y-axis pixel coordinate.
	- const unsigned int &B , Window frame size.
	- double &mean1 , The first mean value to return evaluated using img1.
	- double &mean2 , The second mean value to return evaluated using img2.
	- const unsigned int HEIGHT , The hight of the images.
	- const unsigned int WIDTH , The width of the images.

	Return
	- This function does not return anything.
*/
void getMean(vector<vector<unsigned char>>&img1, vector<vector<unsigned char>> &img2,
	unsigned int x_bias, unsigned int y_bias, const unsigned int &B, double &mean1, double &mean2, 
	const unsigned int HEIGHT, const unsigned int WIDTH)
{
	// Initializing mean values with default value zero in case of it wouldn't be in frame.
	mean1 = 0; mean2 = 0;

	// Initializing axis limits and counter with default value zero.
	int lim_x = 0, lim_y = 0, counter = 0;

	// Evaluating the limits for the beginning of the axis
	lim_y = (y_bias > ((B - 1) / 2)) ? ((B - 1) / 2) : y_bias;
	lim_x = (x_bias > ((B - 1) / 2)) ? ((B - 1) / 2) : x_bias;

	// Calculating the both means inside the window frame.
	for (int x = x_bias - lim_x; (x <= (x_bias + lim_x)) && (x < HEIGHT); x++)
	{
		for (int y = y_bias - lim_y; (y <= (y_bias + lim_y)) && (y < WIDTH); y++)
		{
			mean1 = mean1 + (unsigned char)img1[x][y];
			mean2 = mean2 + (unsigned char)img2[x][y];

			// Counter value will be used in determining how many pixels are inside of the frame.
			counter++;
		}
	}

	// Returning the average.
	if (counter != 0) {
		mean1 /= counter;
		mean2 /= counter;
	}
}

void occlusion_filling_x(vector<vector<zncc_parameters>> &zncc)
{
	cout << "Occlusion filling through x"<<endl;
	for (int i = 0; i < zncc.size(); i++)
	{
		for (int j = 0; j < zncc[0].size(); j++)
		{
			if (zncc[i][j].disparity == 0)
			{
				for (int k = 0; true && (k < zncc[0].size()); k++)
				{
					if ((j + k) < zncc[0].size())
					{
						if (zncc[i][j + k].disparity != 0)
						{
							zncc[i][j].disparity = zncc[i][j + k].disparity;
							break;
						}

					}
					else if ((j - k) > 0)
					{
						if (zncc[i][j - k].disparity != 0)
						{
							zncc[i][j].disparity = zncc[i][j - k].disparity;
							break;
						}
					}


				}
			}

		}
	}
}


void occlusion_filling_y(vector<vector<zncc_parameters>> &zncc)
{
	cout << "Occlusion filling through y"<<endl;
	for (int i = 0; i < zncc.size(); i++)
	{
		for (int j = 0; j < zncc[0].size(); j++)
		{
			if (zncc[i][j].disparity == 0)
			{
				for (int k = 0; true && (k < zncc.size()); k++)
				{
					if ((i + k) < zncc.size())
					{
						if (zncc[i+k][j].disparity != 0)
						{
							zncc[i][j].disparity = zncc[i+k][j].disparity;
							break;
						}

					}
					else if ((i - k) > 0)
					{
						if (zncc[i-k][j].disparity != 0)
						{
							zncc[i][j].disparity = zncc[i-k][j ].disparity;
							break;
						}
					}


				}
			}

		}
	}
}


void two_maps_to_one(vector<vector<zncc_parameters>> &znccX, vector<vector<zncc_parameters>> &znccY)
{
	cout << "Merging two images";
	for (int i = 0; i < znccX.size(); i++)
	{
		for (int j = 0; j < znccX[0].size(); j++)
		{
			int temp1 = znccX[i][j].disparity < znccY[i][j].disparity ? znccY[i][j].disparity : znccX[i][j].disparity;
			int temp2 = znccX[i][j].disparity > znccY[i][j].disparity ? znccY[i][j].disparity : znccX[i][j].disparity;
			znccX[i][j].disparity = int(temp2 / 3 + temp1 * 2 / 3);
		}
	}
}



//Decode from disk to raw pixels with a single function call
void acquireImage(const char* filename, unsigned &width, unsigned &height, vector<unsigned char> &image)
{
	//the raw pixels
	//decode
	unsigned error = lodepng::decode(image, width, height, filename);
	//if there's an error, display it
}

void grayDownSampled(vector <unsigned char> &image, vector <unsigned char> &grayImage)
{
	//downsampling images by a factor of 1/16
	int j = 0;
	int k = 0;
	for (int j = 0; j < 2016; j += 4)
	{
		for (int i = 11760 * j; i < image.size(); i = i + 16)
		{
			if (i == 11760 * (j + 1))
			{
				break;
			}
			unsigned char gray = (unsigned char)(0.2126*image[i] + 0.7152*image[i + 1] + 0.0722*image[i + 2]);
			grayImage.push_back(gray);
		}
	}
}
void toDoubleDimension(vector <unsigned char> &oneDimension, vector< vector <unsigned char> >&twoDimensions)
{
	vector <unsigned char> temp;

	for (int i = 0; i < 504; i++)
	{
		for (int j = 0; j < 735; j++)
		{
			if (j + i * 735 < oneDimension.size())
				temp.push_back(oneDimension[j + i * 735]);

		}
		twoDimensions.push_back(temp);
		temp.erase(temp.begin(), temp.end());
	}

}
void zncc_to_one_dimension_gray(vector<vector<zncc_parameters>> &two_D, vector<unsigned char> &one_D)
{
	const int I = two_D.size(), J = two_D[0].size();
	for (int i = 0; i < I; i++)
	{
		for (int j = 0; j < J; j++)
		{
			for (int k = 0; k < 3; k++)
				one_D.push_back(two_D[i][j].disparity);
			one_D.push_back((unsigned char) 255);
		}
	}
}

void copyVector(vector<vector<zncc_parameters>> &zncc1, vector<vector<zncc_parameters>> &zncc2)
{
	const int I = zncc1.size(), J = zncc1[0].size();
	for (int i = 0; i < I; i++)
	{
		vector<zncc_parameters> temp;
		for (int j = 0; j < J; j++)
		{
			temp.push_back(zncc1[i][j]);
		}
		zncc2.push_back(temp);
	}
}