#ifndef DIGEST_H
#define DIGEST_H

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
//#include "Hybrid.h"

class DigestNode {
private:
	int rows; // number of rows in digest
	int cols; // number of cols in digest
	int o_rows; // resolution of source fingerprint
	int o_cols; 
	int type;

	Mat prnu_ns;
	Mat indices;
public:
	DigestNode() : rows(0), cols(0), o_rows(0), o_cols(0), type(1), prnu_ns(1, 1, 1), indices(1, 1, 1) {}
	//Noise(): rows(0), cols(0), depth(0), channels(0), prnu_ns(1,1,CV_8U){}

	//-------------------------setters/getters----------------------
	int get_rows() { return rows; }
	int get_cols() { return cols; }
	int get_orows() { return o_rows; }
	int get_ocols() { return o_cols; }
	int get_type() { return type; }
	Mat get_prnu() { return prnu_ns; }
	Mat get_indices() { return indices; }

	//void set_noise(Mat mat);
	void set_rows(int r) { rows = r; }
	void set_cols(int c) { cols = c; }
	void set_orows(int or) { o_rows = or; }
	void set_ocols(int oc) { o_cols = oc; }
	void set_type(int t) { type = t; }

	//void set_depth(int d)		{depth = d;		}
	//void set_channels(int ch)	{channels = ch;	}
	void set_prnu(Mat mat) { prnu_ns = mat; }
	void set_indices(Mat mat) { indices = mat; }

	void set_digest(DigestNode dig, const Mat ns, const Mat digest, const Mat indices);

	//-----------------methods--------------------------------------
	static vector<int> get_header(DigestNode dig);

	static void writeDigestVec(string vec_file, vector<DigestNode> vec);
	static void writeDigest(string ns_file, DigestNode dig);

	static void readDigest(const string dig_file, DigestNode& dig);
	static void readDigestVec(string vec_file, vector<DigestNode> &vec);

	static void noise2digest(DigestNode &dig,const Mat ns, const double PERC);
};




#endif