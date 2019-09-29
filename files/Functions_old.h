#pragma once
#ifndef FUNCTIONS_H_old
#define FUNCTIONS_H_old


#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>    
#include <stdlib.h>
#include <stdio.h>
#include <direct.h>
#include <string> 
#include <sstream> 
#include <math.h>

#include <boost/filesystem.hpp>
#define bfs boost::filesystem 


using namespace std;

#include "opencv\cv.h"
#include "opencv\highgui.h"

#include <opencv2\core\mat.hpp> 
#include <opencv2\core\core.hpp>
using namespace cv;

constexpr auto HOUR = 3600000;
constexpr auto MINUTE = 60000;
constexpr auto SEC = 1000;
constexpr auto KILO = 1024;

#define mymax(A,B) (A > B ? A : B)
#define mymin(A,B) (A < B ? A : B)
#define even(x)  ((x & 1) ? 0 : 1)
#define isint(x) ((x - floor(x)) > 0.0 ? 0 : 1)


struct PCE_node {
	double pce_val;
	double peak_value;
	Point peakLocation;
}typedef;

void	rgb2gray(cv::Mat img, cv::Mat& output);
void	threshold(cv::Mat& coef_mat, cv::Mat sub_mat, float noiseVar);
void	waveNoise(cv::Mat initial_mat, cv::Mat& sub_mat, float noiseVar);
void	zeroMeanTotal(cv::Mat& output, cv::Mat input);
void	zeroMean(cv::Mat& output, int x, int y, string type);
//double	log2(double val);
//string	int2str(int val);
string	doub2str(double dbl);

cv::Mat		crossCorrelation(const cv::Mat mat1_real, const cv::Mat mat2_real);
double	pce_op_wo_image(cv::Mat mat1, cv::Mat mat2);

void	intenScale(cv::Mat& out, cv::Mat in);
cv::Mat		shiftRow(cv::Mat in, int numRight);
cv::Mat		shiftCol(cv::Mat in, int numDown);
cv::Mat		fliplr(cv::Mat in);
cv::Mat		flipud(cv::Mat in);
void	mat_arr2gray(const cv::Mat input[], cv::Mat& output);
void	dft_thread(cv::Mat planes[]);
cv::Mat		flip_cross(cv::Mat in);
void	subtract_mean(cv::Mat real_mat, cv::Mat& out_mat);

cv::Mat		removeNeighbourhood(cv::Mat cross_corr, cv::Point peak_loc, int squareSize);
PCE_node get_pce(cv::Mat cross_corr, int col_shift, int row_shift, int square_size);
cv::Mat		crossCorrelationFULL(const cv::Mat mat1_real, const cv::Mat mat2_real);

PCE_node get_pce(cv::Mat cross_corr, int col_shift=0, int row_shift=0, int square_size=0) {
	PCE_node output;

	int cols = cross_corr.cols;
	int rows = cross_corr.rows;
	if (square_size < 1)
		square_size = 11;

	if (row_shift > rows)
		row_shift = rows - 1;
	else if (row_shift < 0)
		row_shift = 0;

	if (col_shift > cols)
		col_shift = cols - 1;
	else if (col_shift < 0)
		col_shift = 0;

	if (cross_corr.at<double>(0, 0) == 0) {			//to make it more efficient
		if (countNonZero(cross_corr) == 0) {			//if the energy is of cross_correlation is zero
			//out icin bi typedef yapmak lazim 
			output.pce_val = 0;
			output.peak_value = 1;
			output.peakLocation = Point(0, 0);
			return output;
		}
	}

	cv::Mat c_inrange = cross_corr(cv::Rect(cols - col_shift - 1, rows - row_shift - 1, col_shift + 1, row_shift + 1));

	cv::MatIterator_<double> fw = max_element(c_inrange.begin<double>(), c_inrange.end<double>());
	double max_cc = *fw;			//peakheight

	double minVal = -10000000;
	double maxVal = 10000000000;
	cv::Point maxLoc;				//xpeak, ypeak burasi
	cv::Point minLoc;
	minMaxLoc(c_inrange, &minVal, &maxVal, &minLoc, &maxLoc);
	//	int index_max = maxLoc.x * c_inrange.cols + maxLoc.y+1;

	double peak_height = max_cc;
	output.peakLocation = cv::Point(row_shift - maxLoc.x, col_shift - maxLoc.y);

	cv::Mat cor_peakless = removeNeighbourhood(cross_corr, maxLoc, square_size);
	//	cout << output.peakLocation << endl;

	double correl = cross_corr.at<double>(rows - 1, cols - 1);
	pow(cor_peakless, 2.0, cor_peakless);

	double pce_energy;	cv::Scalar scl1, scl2;
	meanStdDev(cor_peakless, scl1, scl2);
	pce_energy = scl1[0];
	//	cout << pce_energy << "pce_energ\n";

	output.peak_value = pow(peak_height, 2.0) / pce_energy * (peak_height / abs(peak_height));
	//	cout << output.peak_value << endl;

	double tempval = peak_height / sqrt(pce_energy) / sqrt(2.0);
	//	cout << tempval << endl;
	//	output.pce_val = 1/2;
		//Out.pvalue = 1/2*erfc(peakheight/sqrt(PCE_energy)/sqrt(2));     % under simplifying assumption that C are samples from Gaussian pdf 

	//  [Out.P_FA,Out.log10P_FA] = FAfromPCE(Out.PCE,prod(shift_range+1));
	return output;
}

cv::Mat removeNeighbourhood(cv::Mat cross_corr, cv::Point peak_loc, int square_size) {
	int nrows = cross_corr.rows;
	int ncols = cross_corr.cols;

	int radius = (square_size - 1) / 2;

	cross_corr = shiftCol(cross_corr, radius - peak_loc.y);
	cross_corr = shiftRow(cross_corr, radius - peak_loc.x);

	cv::Mat Y(nrows - square_size, square_size, CV_64F);
	cross_corr(cv::Rect(0, square_size, square_size, nrows - square_size)).copyTo(Y(cv::Rect(0, 0, square_size, nrows - square_size)));

	transpose(Y, Y);
	Y.resize(1, 1);

	cv::Mat transpose_ccorr = cross_corr.t();

	transpose_ccorr = transpose_ccorr.reshape(1, 1);
	cv::Mat maymat = transpose_ccorr(cv::Range(0, 1), cv::Range(square_size * nrows + 1, transpose_ccorr.size().area()));

	hconcat(Y, maymat, Y);

	return Y;
}

void FAfromPCE(int pce, cv::Mat mymat, double& p_fa, double& log10p_fa) {}

void rot90(cv::Mat& matImage, int rotflag) {
	//1=CW, 2=CCW, 3=180
	if (rotflag == 1) {
		transpose(matImage, matImage);
		flip(matImage, matImage, 1); //transpose+flip(1)=CW
	}
	else if (rotflag == 2) {
		transpose(matImage, matImage);
		flip(matImage, matImage, 0); //transpose+flip(0)=CCW     
	}
	else if (rotflag == 3) {
		flip(matImage, matImage, -1);    //flip(-1)=180          
	}
	else if (rotflag != 0) { //if not 0,1,2,3:
		cout << "Unknown rotation flag(" << rotflag << ")" << endl;
	}
}

void rotate_CW(cv::Mat& img) {
	cv::transpose(img, img);
	cv::flip(img, img, 1);
}

void rotate_CCW(cv::Mat& img) {
	cv::transpose(img, img);
	cv::flip(img, img, 0);
}


//double log2(double val) {
//	return log10(val) / log10(2.0);
//}

//string int2str(int val){
//	return static_cast<ostringstream*>( &(ostringstream() << val) )->str();
//}

string doub2str(double dbl) {
	std::ostringstream strs;
	strs << dbl;
	return strs.str();
}

void bgr2gray(cv::Mat input, cv::Mat& output) {
	if (input.size().height == 0 && input.size().width == 0)
		cout << "hi";//here there is no image
	cvtColor(input, output, CV_BGR2GRAY);
}

void rgb2gray(cv::Mat input, cv::Mat& output) {
	if (input.size().height == 0 && input.size().width == 0)
		cout << "hi";//here there is no image
	cvtColor(input, output, CV_RGB2GRAY);
}

//void intenScale(cv::Mat& out, cv::Mat in) {
//	int thresh = 252;
//	int V = 6;
//
//	cv::Mat mask = (out > thresh) / 255;			//12-13ms
//	cv::Mat temp;
//
//	subtract(out, thresh, temp, mask);
//	pow(temp, 2.0, temp);
//	multiply(temp, (-1.0 / 6), temp);
//	exp(temp, temp);
//
//	multiply(temp, mask, temp);
//
//	subtract(1, mask, mask);
//	//252den kucuk olanlar icin
//	multiply(out, mask, out);
//	divide(out, thresh, out);
//	add(out, temp, out);
//}


void intenScale(cv::Mat& out, cv::Mat in) {
	int thresh = 252;
	int V = 6;
	if (in.type() != 5)
		in.convertTo(out, CV_32F);

	cv::Mat mask = (out > thresh) / 255;			//12-13ms
	cv::Mat temp;

	cv::subtract(out, thresh, temp, mask);
	cv::pow(temp, 2.0, temp);
	cv::multiply(temp, (-1.0 / 6), temp);
	cv::exp(temp, temp);

	if (mask.type() != 5)
		mask.convertTo(mask, CV_32F);			//5 ms
	cv::multiply(temp, mask, temp);

	cv::subtract(1, mask, mask);
	//252den kucuk olanlar icin
	cv::multiply(out, mask, out);
	cv::divide(out, thresh, out);
	cv::add(out, temp, out);
}


cv::Mat saturation(cv::Mat& in, String isGray) {
	int rows = in.rows;
	int cols = in.cols;
	cv::Mat saturMat; cv::Mat temp_in;
	if (in.type() != 5)
		in.convertTo(temp_in, CV_32F);

	cv::MatIterator_<uchar> fw = max_element(in.begin<uchar>(), in.end<uchar>());
	if ((int)* fw < 250) {
		saturMat = cv::Mat::ones(in.size(), CV_32F);
		return saturMat;
	}
	cv::Mat Xh; cv::Mat Xv;
	Xh = shiftCol(temp_in, 1);
	cv::subtract(temp_in, Xh, Xh);
	Xv = shiftRow(temp_in, 1);
	cv::subtract(temp_in, Xv, Xv);

	cv::Mat Xa;
	cv::multiply(Xh, Xv, Xa);
	cv::multiply(Xa, shiftCol(Xh, -1), Xa);
	cv::multiply(Xa, shiftRow(Xv, -1), Xa);

	cv::Mat mask = (cv::abs(Xa) > 0) / 255;

	if (in.channels() == 3) {
		cout << "channels 3 in saturation which is impossible in saturation" << endl;
		/*		//matlab code is the following
		for j=1:3
			maxX(j) = max(max(X(:,:,j)));
			if maxX(j)>250;
				SaturMap(:,:,j) = ~((X(:,:,j)==maxX(j)) & ~SaturMap(:,:,j));
			end
		end
		*/
	}
	else if (in.channels() == 1) {
		cv::MatIterator_<uchar> max_elem = max_element(in.begin<uchar>(), in.end<uchar>());
		cv::Mat newmask = (in == *fw) / 255;
		cv::subtract(1, mask, mask);
		cv::multiply(mask, newmask, mask);
		cv::subtract(1, mask, mask);
	}
	else { cerr << "invalid dimension " << endl; }

	if (isGray == "gray") {
		cout << " isGray is never gray in saturation" << endl;
		/*matlab code
		if size(X,3)==3,
			SaturMap = SaturMap(:,:,1)+SaturMap(:,:,3)+SaturMap(:,:,3);
			SaturMap(SaturMap>1) = 1;
		*/
	}
	if (mask.type() != 5)
		mask.convertTo(mask, CV_32F);
	return mask;
}


//numRight = number of shifts to right
cv::Mat shiftCol(cv::Mat in, int numRight) {
	int ncols = in.cols;
	int nrows = in.rows;
	cv::Mat out;
	out = cv::Mat::zeros(in.size(), in.type());

	numRight = numRight % ncols;
	if (numRight < 0)
		numRight = ncols + numRight;

	in(cv::Rect(ncols - numRight, 0, numRight, nrows)).copyTo(out(cv::Rect(0, 0, numRight, nrows)));
	in(cv::Rect(0, 0, ncols - numRight, nrows)).copyTo(out(cv::Rect(numRight, 0, ncols - numRight, nrows)));

	return out;
}
//shift down
cv::Mat shiftRow(cv::Mat in, int numDown) {
	int ncols = in.cols;
	int nrows = in.rows;
	cv::Mat out = cv::Mat::zeros(in.size(), in.type());

	numDown = numDown % nrows;
	if (numDown < 0)
		numDown = nrows + numDown;

	in(cv::Rect(0, 0, ncols, nrows - numDown)).copyTo(out(cv::Rect(0, numDown, ncols, nrows - numDown)));
	in(cv::Rect(0, nrows - numDown, ncols, numDown)).copyTo(out(cv::Rect(0, 0, ncols, numDown)));
	return out;
}

void mat_arr2gray(const cv::Mat input[], cv::Mat& output) {
	if (input[0].size().height == 0 && input[0].size().width == 0)
		cout << "hi";//here there is no input matrice
	cv::Mat temp(input[0].size(), CV_32F);// temp.create(input[0].size(),input[0].type());

	//25 ms
	cv::multiply(input[0], 0.299, output);
	cv::multiply(input[1], 0.587, temp);
	cv::add(temp, output, output);
	cv::multiply(input[2], 0.114, temp);
	cv::add(temp, output, output);
	temp.~Mat();
}

void threshold(cv::Mat& coef_mat, cv::Mat sub_mat, float noiseVar) {
	cv::subtract(sub_mat, noiseVar, sub_mat);
	cv::threshold(sub_mat, coef_mat, 0.0f, 1000000.0f, THRESH_TOZERO);
}



double pce_op_wo_image(cv::Mat mat1, cv::Mat mat2) {
	if (mat1.size().area() < 1 || mat2.size().area() < 1)
		return -1;

	clock_t s = clock();
	cv::Mat c_corra = crossCorrelation(mat1, mat2);
	//cv::Mat c_corra2 = crossCorrelationFAST(mat1, mat2);
	if (c_corra.type() != CV_64F)
		c_corra.convertTo(c_corra, CV_64F);

	return get_pce(c_corra, 0, 0, 0).peak_value;
}


cv::Mat crossCorrelationFULL(const cv::Mat mat1_real, const cv::Mat mat2_real) {
	cv::Mat mat1, mat2;
	cv::Mat result(mat1_real.size(), CV_32F);
	cv::Mat complexI(mat1_real.size(), CV_32FC2);
	mat1_real.copyTo(mat1);
	mat2_real.copyTo(mat2);

	cv::Scalar mean1 = mean(mat1_real);
	cv::Scalar mean2 = mean(mat2_real);

	cv::subtract(mat1, (float)mean1[0], mat1);
	cv::subtract(mat2, (float)mean2[0], mat2);

	cv::Mat tilted_mat2 = flip_cross(mat2);

	mat1.convertTo(mat1, CV_64F);
	cv::Mat planes[] = { cv::Mat_<float>(mat1), cv::Mat::zeros(mat1.size(), CV_32F) };
	cv::Mat planes2[] = { cv::Mat_<float>(tilted_mat2), cv::Mat::zeros(mat1.size(), CV_32F) };

	cv::merge(planes, 2, complexI);						//4		// Add to the expanded another plane with zeros
	cv::dft(complexI, complexI);						//72	// this way the result may fit in the source matrix
	cv::split(complexI, planes);						//1		// planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))

	cv::merge(planes2, 2, complexI);						//4		// Add to the expanded another plane with zeros
	cv::dft(complexI, complexI);						//72	// this way the result may fit in the source matrix
	cv::split(complexI, planes2);						//1		// planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))

	cv::multiply(planes[0], planes2[0], mat1);
	cv::multiply(planes[1], planes2[1], mat2);
	cv::subtract(mat1, mat2, mat1);

	cv::multiply(planes[1], planes2[0], mat2);
	cv::multiply(planes[0], planes2[1], tilted_mat2);
	cv::add(tilted_mat2, mat2, planes[1]);

	planes[0] = mat1;

	cv::merge(planes, 2, complexI);

	cv::idft(complexI, result, DFT_SCALE | DFT_REAL_OUTPUT);

	planes[0].~Mat();
	planes[1].~Mat();
	planes2[0].~Mat();
	planes2[1].~Mat();
	tilted_mat2.~Mat();
	mat1.~Mat();
	mat2.~Mat();
	return result;
}


cv::Mat crossCorrelation(const cv::Mat mat1_real, const cv::Mat mat2_real) {
	cv::Mat mat1, mat2;
	cv::Mat result(mat1_real.size(), CV_32F);
	cv::Mat complexI(mat1_real.size(), CV_32FC2);
	mat1_real.copyTo(mat1);
	mat2_real.copyTo(mat2);

	cv::Scalar mean1 = mean(mat1_real);
	cv::Scalar mean2 = mean(mat2_real);

	cv::subtract(mat1, (float)mean1[0], mat1);
	cv::subtract(mat2, (float)mean2[0], mat2);

	cv::Mat tilted_mat2 = flip_cross(mat2);

	mat1.convertTo(mat1, CV_64F);
	cv::Mat planes[] = { cv::Mat_<float>(mat1), cv::Mat::zeros(mat1.size(), CV_32F) };
	cv::Mat planes2[] = { cv::Mat_<float>(tilted_mat2), cv::Mat::zeros(mat1.size(), CV_32F) };

	cv::merge(planes, 2, complexI);						//4		// Add to the expanded another plane with zeros
	cv::dft(complexI, complexI);						//72	// this way the result may fit in the source matrix
	cv::split(complexI, planes);						//1		// planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))

	cv::merge(planes2, 2, complexI);						//4		// Add to the expanded another plane with zeros
	cv::dft(complexI, complexI);						//72	// this way the result may fit in the source matrix
	cv::split(complexI, planes2);						//1		// planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))

	cv::multiply(planes[0], planes2[0], mat1);
	cv::multiply(planes[1], planes2[1], mat2);
	cv::subtract(mat1, mat2, mat1);

	cv::multiply(planes[1], planes2[0], mat2);
	cv::multiply(planes[0], planes2[1], tilted_mat2); 
	cv::add(tilted_mat2, mat2, planes[1]);

	planes[0] = mat1;

	cv::merge(planes, 2, complexI);

	cv::idft(complexI, result, DFT_SCALE | DFT_REAL_OUTPUT);

	planes[0].~Mat();
	planes[1].~Mat();
	planes2[0].~Mat();
	planes2[1].~Mat();
	tilted_mat2.~Mat();
	mat1.~Mat();
	mat2.~Mat();
	return result;
}


void dft_thread(cv::Mat planes[]) {
	cv::Mat complexI;
	cv::merge(planes, 2, complexI);						//4		// Add to the expanded another plane with zeros
	cv::dft(complexI, complexI);						//72	// this way the result may fit in the source matrix
	cv::split(complexI, planes);						//1		// planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
}

void subtract_mean(cv::Mat mat_real, cv::Mat& out_mat) {
	cv::Scalar mean1 = cv::mean(mat_real);
	cv::subtract(mat_real, (float)mean1[0], out_mat);
}


cv::Mat fliplr(cv::Mat in) {
	int ncols = in.cols;
	int nrows = in.rows;
	cv::Mat out;
	out = cv::Mat::zeros(in.size(), in.type());

	for (int i = 0; i < ncols; i++)
		in(cv::Rect(i, 0, 1, nrows)).copyTo(out(cv::Rect(ncols - 1 - i, 0, 1, nrows)));
	return out;
}

cv::Mat flipud(cv::Mat in) {
	int ncols = in.cols;
	int nrows = in.rows;
	cv::Mat out;
	out = cv::Mat::zeros(in.size(), in.type());

	for (int i = 0; i < nrows; i++)
		in(cv::Rect(0, i, ncols, 1)).copyTo(out(cv::Rect(0, nrows - 1 - i, ncols, 1)));

	return out;
}

cv::Mat flip_cross(cv::Mat in) {
	int nrow = in.rows;
	int ncol = in.cols;
	cv::Mat out(in.rows, in.cols, in.type());

	if (in.type() == 5) {
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++)
				out.at<float>(i, j) = in.at<float>(nrow - 1 - i, ncol - 1 - j);
		}
	}
	else if (in.type() == 6) {
		for (int i = 0; i < nrow; i++) {
			for (int j = 0; j < ncol; j++)
				out.at<float>(i, j) = in.at<float>(nrow - 1 - i, ncol - 1 - j);
		}
	}

	return out;
}

void waveNoise(cv::Mat initial_mat, cv::Mat& out_mat, float noiseVar) {
	int h_size = initial_mat.size().width;
	int v_size = initial_mat.size().height;

	cv::Mat sub_mat;

	pow(initial_mat, 2.0, sub_mat);						//6

	cv::Mat newmat = cv::Mat::ones(3, 3, CV_32F);
	cv::divide(newmat, -3 * 3, newmat, -1);					//1

	cv::Mat coef_mat;// coef_mat.create(Size(v_size,h_size),sub_mat.type());
	filter2D(sub_mat, coef_mat, sub_mat.depth(), newmat, Point(-1, -1), 0.0, BORDER_ISOLATED); //21   //504-376 iken 6 sec

	threshold(coef_mat, coef_mat, noiseVar);

	cv::Mat est_mat;								//0

	clock_t begin = clock();

	for (int w = 5; w <= 9; w += 2) {				//1480 ms 1000-750 icin
		newmat = cv::Mat::ones(w, w, CV_32F);
		cv::divide(newmat, -w * w, newmat, -1);

		clock_t b1 = clock();
		filter2D(sub_mat, est_mat, sub_mat.depth(), newmat, Point(-1, -1), 0.0, BORDER_ISOLATED);
		clock_t b2 = clock();
		//		cout << "filter 2d :" << b2 - b1 << endl;
		threshold(est_mat, est_mat, noiseVar);
		clock_t b3 = clock();
		//		cout << "thresh :" << b3 - b2 << endl;
		cv::min(coef_mat, est_mat, coef_mat);
	}

	clock_t end1 = clock();

	cv::multiply(initial_mat, noiseVar, sub_mat);			//15 ms alttaki ucu
	cv::add(coef_mat, noiseVar, coef_mat);
	cv::divide(sub_mat, coef_mat, out_mat);					//bu divide'in outputu out_mat olursa son isleme gerek kalmaz

	newmat.~Mat();
	est_mat.~Mat();
	sub_mat.~Mat();
	coef_mat.~Mat();
}

void zeroMeanTotal(cv::Mat& output, cv::Mat input) {
	int nrow = input.rows;
	int ncol = input.cols;

	cv::Mat out = cv::Mat::zeros(input.size(), input.type());

	int in_channels = input.channels();

	if (in_channels == 1) { //gray
		zeroMean(input, 0, 0, "both");
		zeroMean(input, 0, 1, "both");
		zeroMean(input, 1, 0, "both");
		zeroMean(input, 1, 1, "both");
	}

	else {
		for (int i = 0; i < in_channels; i++) {
			//do zero mean total
		}
	}
}

void zeroMean(cv::Mat& input, int x, int y, string type) {
	int nrow = input.rows;
	int ncol = input.cols;

	if (type != "col") {			//row oldu bu
		for (int i = x; i < nrow; i += 2) {
			float sum = 0;
			for (int j = y; j < ncol; j += 2) {
				sum += input.at<float>(i, j);
			}
			sum /= (ncol + 1 - x) / 2;
			for (int j = y; j < ncol; j += 2) {
				input.at<float>(i, j) -= sum;
			}
		}
	}

	if (type != "row") {			//col oldu bu
		for (int i = y; i < ncol; i += 2) {
			float sum = 0;
			for (int j = x; j < nrow; j += 2) {
				sum += input.at<float>(j, i);
			}
			sum /= ((nrow + 1 - y) / 2);
			for (int j = x; j < nrow; j += 2) {
				input.at<float>(j, i) -= sum;
			}
		}
	}

}

std::vector<std::string>& split(const std::string& s, char delim, std::vector<std::string>& elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

std::vector<std::string> split(const std::string& s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

void trim(std::string& str) {
	str.erase(0, str.find_first_not_of(' '));       //prefixing spaces
	str.erase(str.find_last_not_of(' ') + 1);         //surfixing spaces
}

void replaceStringInPlace(std::string& subject, const std::string& search,
	const std::string& replace) {
	size_t pos = 0;
	while ((pos = subject.find(search, pos)) != std::string::npos) {
		subject.replace(pos, search.length(), replace);
		pos += replace.length();
	}
}

string sanitize_dir(string in_dir) {
	replaceStringInPlace(in_dir, "/", "\\");
	if (in_dir.at(in_dir.length() - 1) != '\\')
		in_dir += "\\";
	return in_dir;
}

void printXbyY(const cv::Mat matrix, int x, int y) {
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++)
			cout << matrix.at<float>(i, j) << " ";
		cout << endl;
	}
}

bool isprnuOK(const cv::Mat ns) {
	if (ns.size().area() <= 10)
		return false;
	return true;
}

string printable_fingerprint(const string fp_name) {
	vector<string> fp_info = split(fp_name, ';');
	string str = fp_name;

	if (fp_info.size() >= 3)
		str = fp_info.at(0) + "\n" + fp_info.at(1) + "\n" + fp_info.at(2) + "\n";

	return str;
}

string get_name(string img_file) {
	auto img_name_start = img_file.find_last_of('\\');
	auto img_name_end = img_file.find_last_of('.');

	string img_name = img_file.substr(img_name_start + 1, (img_name_end - img_name_start - 1));
	return img_name;
}


void list_images(string txt_file, vector<string>& image_list) {
	ifstream myFile;
	myFile.open(txt_file);
	string line;
	if (myFile.is_open())
		while (getline(myFile, line)) {
			image_list.push_back(get_name(line));
		}
	myFile.close();
}


#endif