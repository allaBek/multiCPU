#include <iostream>
#include "lodepng.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream> 
#include <string>
#include <chrono>
#include "functions.h"

using namespace std;
//declaring data types



int main()
{
	cout << "Starting the code here" << endl; 
	auto start = chrono::high_resolution_clock::now();
    //declaring parameters to acquire gray down-sampled image
    unsigned width, height;
    //file path for the two images
    const char* fileName1 = "images\\im0.png";
    const char* fileName2 = "images\\im1.png";
    vector<unsigned char> img1, img2, grayImg1, grayImg2;
    vector< vector <unsigned char> > sample1, sample2;
    //obtaining image
    acquireImage(fileName1, width, height, img1);
    acquireImage(fileName2, width, height, img2);
    //down-sampling and getting the gray scale image
    grayDownSampled(img1, grayImg1);
    grayDownSampled(img2, grayImg2);
	//clear unneeded vectors to free up memory
	img1.clear();
	img1.shrink_to_fit();
	img2.clear();
	img2.shrink_to_fit();
	//converting images from vector to 2D vector
    toDoubleDimension(grayImg1, sample1);
    toDoubleDimension(grayImg2, sample2);

	//clear vectors
	grayImg1.clear();
	grayImg1.shrink_to_fit();
	grayImg2.clear();
	grayImg2.shrink_to_fit();
	//start the ZNCC algorithmm
    int d = 0;
    vector<vector<zncc_parameters>> zncc1, zncc2, zncc3;
	vector<vector<zncc_parameters>> zncc3_copy;
    ZNCC(sample1, sample2, zncc1, zncc2, zncc3);
	copyVector(zncc3, zncc3_copy);
	occlusion_filling_x(zncc3);
	occlusion_filling_y(zncc3_copy);
	two_maps_to_one(zncc3, zncc3_copy);
	auto finish = chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Elapsed time: " << elapsed.count() << " s\n";
	//creating output
	vector<unsigned char> d1, d2, d3;
	const char* zncc1_path = "images\\zncc1.png";
	const char* zncc2_path = "images\\zncc2.png";
	const char* zncc3_path = "images\\zncc3.png";
	//convert
	zncc_to_one_dimension_gray(zncc1, d1);
	zncc_to_one_dimension_gray(zncc2, d2);
	zncc_to_one_dimension_gray(zncc3, d3);
	//output png
	lodepng::encode(zncc1_path, d1, sample1[0].size(), sample1.size());
	lodepng::encode(zncc2_path, d2, sample1[0].size(), sample1.size());
	lodepng::encode(zncc3_path, d3, sample1[0].size(), sample1.size());
	system("pause");
	return 0;


	
}
