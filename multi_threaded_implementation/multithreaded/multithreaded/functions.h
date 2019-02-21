#pragma once
 // for string streams
using namespace std;
//declare structure to hold data
struct zncc_parameters {
	int x = 0;
	int y = 0;
	int disparity = 0;
	double zncc = -10;
};

// functions signatures
void     acquireImage(const char* filename, unsigned &width, unsigned &height, vector<unsigned char> &image);
void grayDownSampled(vector <unsigned char> &image, vector <unsigned char> &grayImage);
void toDoubleDimension(vector <unsigned char> &oneDimension, vector< vector <unsigned char> >&twoDimensions);
void getMean(vector< vector <unsigned char> >&img1, vector< vector <unsigned char> >&img2,
	unsigned int x_bias, unsigned int y_bias, const unsigned int &B, double &mean1, double &mean2);
void ZNCC(vector< vector <unsigned char> >&sample1, vector< vector <unsigned char> >&sample2, vector<vector<zncc_parameters>> &zncc1, vector<vector<zncc_parameters>> &zncc2, vector<vector<zncc_parameters>> &zncc3);
void zncc_to_one_dimension_gray(vector<vector<zncc_parameters>> &two_D, vector<unsigned char> &one_D);
void occlusion_filling_x(vector<vector<zncc_parameters>> &zncc);
void occlusion_filling_y(vector<vector<zncc_parameters>> &zncc);
void two_maps_to_one(vector<vector<zncc_parameters>> &znccX, vector<vector<zncc_parameters>> &znccY);
void processing_per_pixel(vector< vector <unsigned char> >&sample1, vector< vector <unsigned char> >&sample2, zncc_parameters &zncc1, zncc_parameters &zncc2, zncc_parameters &zncc3, int x, int y);
void copyVector(vector<vector<zncc_parameters>> &zncc1, vector<vector<zncc_parameters>> &zncc2);
void processing_per_line(vector< vector <unsigned char> >&sample1, vector< vector <unsigned char> >&sample2, vector<zncc_parameters> &zncc1, vector<zncc_parameters> &zncc2, vector<zncc_parameters> &zncc3, int y);