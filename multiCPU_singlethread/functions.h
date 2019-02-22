using namespace std;

// Declaring a structure to hold data for ZNCC processing.
struct zncc_parameters {
	int x = 0;	// The x-axis of pixel, initialized as zero.
	int y = 0; // The y-axis of pixel, initialized as zero.
	int disparity = 0; // The disparity value, initialized as zero.
	double zncc = -10; // The ZNCC value, initialized as -10.
};

// Function signatures
void acquireImage(const char* filename, unsigned &width, unsigned &height, vector<unsigned char> &image);
void grayDownSampled(vector<unsigned char> &image, vector<unsigned char> &grayImage, unsigned int HEIGHT, unsigned int WIDTH);
void toDoubleDimension(vector<unsigned char> &oneDimension, vector<vector<unsigned char>> &twoDimensions, unsigned int HEIGHT, unsigned int WIDTH);
void getMean(vector<vector<unsigned char>>&img1, vector<vector<unsigned char>> &img2,
	unsigned int x_bias, unsigned int y_bias, const unsigned int &B, double &mean1, double &mean2, 
	const unsigned int HEIGHT, const unsigned int WIDTH);
void ZNCC(vector<vector<unsigned char>>&sample1, vector<vector<unsigned char>>&sample2, vector<vector<zncc_parameters>> &zncc1, 
	vector<vector<zncc_parameters>> &zncc2, vector<vector<zncc_parameters>> &zncc3);
void zncc_to_one_dimension_gray(vector<vector<zncc_parameters>> &two_D, vector<unsigned char> &one_D);
void occlusion_filling_x(vector<vector<zncc_parameters>> &zncc);
void occlusion_filling_y(vector<vector<zncc_parameters>> &zncc);
void two_maps_to_one(vector<vector<zncc_parameters>> &znccX, vector<vector<zncc_parameters>> &znccY);
void copyVector(vector<vector<zncc_parameters>> &zncc1, vector<vector<zncc_parameters>> &zncc2);