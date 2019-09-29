#ifndef NOISE_H
#define NOISE_H

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "opencv\cv.h"
#include "opencv\highgui.h"

#include <opencv2\core\mat.hpp> 
#include <opencv2\core\core.hpp>
using namespace cv;

#include <boost/filesystem.hpp>
#include <vector>
#include <numeric>  
#include <algorithm>
#include "RC_Pair.h"
#include "NoiseExtractFromImage.h"

class NoiseNode{
private:
	int rows;
	int cols;
	int type;
	//int depth;
	//int channels;
	Mat prnu_ns;
	
public:
	NoiseNode(): rows(0), cols(0), type(1), prnu_ns(1,1,1){}
	//Noise(): rows(0), cols(0), depth(0), channels(0), prnu_ns(1,1,CV_8U){}
	NoiseNode(int rows, int cols, int type);

	//-------------------------setters/getters----------------------
	int get_rows()			{return rows;		}
	int get_cols()			{return cols;		}
	int get_type()			{return type;		}

	//int get_depth()			{return depth;		}
	//int get_channels()		{return channels;	}
	Mat get_prnu()			{return prnu_ns;	}

	void set_noise_node(Mat mat);

	void set_rows(int r)		{rows = r;		}
	void set_cols(int c)		{cols = c;		}
	void set_type(int t)		{type = t;		}
	void set_prnu(Mat mat)		{prnu_ns = mat;	}

	//-----------------methods--------------------------------------
	void writeNoise(string ns_file, Mat ns);

	void readNoise(string ns_file, Mat &ns);
	static void NoiseNode::readNoise(string ns_file, NoiseNode& ns);
	static void writeNoise(string ns_file, NoiseNode nn);
	void printXbyY(Mat a, int x, int y);

	static RC_Pair get_rc_from_ns_file(string ns_file);

	//static NoiseNode rgb2NoiseNode(string imageName) {
	//	NoiseNode nn = NoiseNode();
	//	if (!boost::filesystem::is_regular_file(imageName))
	//		return nn;

	//	Mat noise;
	//	float sigma = 2;
	//	noise = noiseExtractionFromImage(imageName, sigma, false, false);
	//	Scalar mean = 0, stdDev = 0;

	//	meanStdDev(noise, mean, stdDev);
	//	Mat noisex;

	//	wienerInDFT(noisex, noise, (float)stdDev.val[0]);
	//	noise.~Mat();
	//	nn.set_noise_node(noisex);
	//	return nn;
	//}
	bool operator< (const NoiseNode& other) const {
		if(rows*cols == other.rows*other.cols)
			return rows < other.rows;
		else
			return rows * cols < other.rows * other.cols;
	}
};




#endif