/*
Authors : Mohamed Bekoucha & Berke Esmer
For the course, 521288S Multiprogramming Project

University of Oulu, Spring 2019
*/

#include <iostream>
#include <fstream>
#include <sstream> 
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include "lodepng.h"
#include "functions.h"

using namespace std;

// Main function is the only runnable part of the application.
int main()
{
	cout << "Starting the code here" << endl; 

	// Starting the clock time in order to measure statistics.
	auto start = chrono::high_resolution_clock::now();

	// Data structures to store raw and processed image data.
	vector<unsigned char> img1, img2, grayImg1, grayImg2, out1, out2, out3, out4, out5;
	vector<vector<unsigned char>> sample1, sample2;

    // Declaring parameters to acquire size of images.
    unsigned width, height;

	// File path for the two images and the outputs.
    const char* fileName1 = ".\\images\\im0.png";
    const char* fileName2 = ".\\images\\im1.png";
	const char* zncc1_path = ".\\images\\zncc1.png";
	const char* zncc2_path = ".\\images\\zncc2.png";
	const char* zncc3_path = ".\\images\\zncc3.png";

	// Obtaining the raw image data as vector.
    acquireImage(fileName1, width, height, img1);
    acquireImage(fileName2, width, height, img2);

	// Down-sampling and getting the gray-scale image by processing the raw image data.
	grayDownSampled(img1, grayImg1, height, width);
    grayDownSampled(img2, grayImg2, height, width);

	// Clearing unnecessary vectors to free up memory.
	img1.clear(); img1.shrink_to_fit();
	img2.clear(); img2.shrink_to_fit();

	// Converting images from one dimensional vector to 2D vector.
	toDoubleDimension(grayImg1, sample1, height / 4, width / 4);
    toDoubleDimension(grayImg2, sample2, height / 4, width / 4);

	// Clearing unnecessary vectors to free up memory.
	grayImg1.clear(); grayImg1.shrink_to_fit();
	grayImg2.clear(); grayImg2.shrink_to_fit();

	// Decleration of data structure to process ZNCC algorithm.
	vector<vector<zncc_parameters>> zncc1, zncc2, zncc3, zncc3_x_filling, zncc3_y_filling, zncc3_x_filling_copy;

	// Processing the ZNCC algorithm.
    ZNCC(sample1, sample2, zncc1, zncc2, zncc3);
	
	// Creating copies of vector zncc3, just for size template.
	copyVector(zncc3, zncc3_x_filling);
	copyVector(zncc3, zncc3_y_filling);
	copyVector(zncc3, zncc3_x_filling_copy);

	// Clearing unnecessary vectors to free up memory.
	sample1.clear(); sample1.shrink_to_fit();
	sample2.clear(); sample2.shrink_to_fit();

	// Post-processing on data using two different methods.
	occlusion_filling_x(zncc3_x_filling);
	occlusion_filling_y(zncc3_y_filling);

	// Joining two results into one to reach the final conclusion.
	two_maps_to_one(zncc3_x_filling_copy, zncc3_y_filling);

	// Conversion back from ZNCC structure to one-dimensional vector form.
	zncc_to_one_dimension_gray(zncc1, out1);
	zncc_to_one_dimension_gray(zncc2, out2);
	zncc_to_one_dimension_gray(zncc3_x_filling, out3);
	zncc_to_one_dimension_gray(zncc3_y_filling, out4);
	zncc_to_one_dimension_gray(zncc3_x_filling_copy, out5);

	// Clearing unnecessary vectors to free up memory.
	zncc1.clear(); zncc1.shrink_to_fit();
	zncc2.clear(); zncc2.shrink_to_fit();
	zncc3.clear(); zncc3.shrink_to_fit();
	zncc3_x_filling.clear(); zncc3_x_filling.shrink_to_fit();
	zncc3_y_filling.clear(); zncc3_y_filling.shrink_to_fit();
	zncc3_x_filling_copy.clear(); zncc3_x_filling_copy.shrink_to_fit();

	// Obtaining the .png images using the 1D vectors.
	lodepng::encode(".\\images\\first_disparity_map.png", out1, width / 4, height/4);
	lodepng::encode(".\\images\\second_disparity_map.png", out2, width /4, height/4);
	lodepng::encode(".\\images\\X_occlusion_filling.png", out3, width /4, height/4);
	lodepng::encode(".\\images\\Y_occlusion_filling.png", out4, width / 4, height / 4);
	lodepng::encode(".\\images\\X_Y_occlusion_filling.png", out5, width / 4, height / 4);

	// Clearing unnecessary vectors to free up memory.
	out1.clear(); out1.shrink_to_fit();
	out2.clear(); out2.shrink_to_fit();
	out3.clear(); out3.shrink_to_fit();

	// Stopping the clock time and retrieving elapsed time.
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;

	// Printing statistics.
	cout << "Elapsed time: " << elapsed.count() << " s\n";

	// Program exit.
	system("pause");
	return 0;
}
