#ifndef TEST_CORR_H
#define TEST_CORR_H

#include "Functions.h"
//#include "test_functions.h"


// This class will contain the following code pieces:
// 1) corr2
// 2) pce
// 3) crosscorr?
// 4) NCC
// 5) sd_search / sd_corr
#include <opencv/cv.h>
#include <opencv2/core/core.hpp>
using namespace cv;
#include "Digest.h"
#include "Hybrid.h"

using namespace std;

class Correlation {
public:
	static double corr2(Mat mat_x, Mat mat_y) {
		Scalar mean_x, mean_y, std_x, std_y;
		//mean_x and mean_y determined

		meanStdDev(mat_x, mean_x, std_x);
		meanStdDev(mat_y, mean_y, std_y);

		cv::Mat temp_x, temp_y;
		subtract(mat_x, (float)mean_x[0], temp_x);
		subtract(mat_y, (float)mean_y[0], temp_y);

		meanStdDev(temp_x, mean_x, std_x);
		meanStdDev(temp_y, mean_y, std_y);

		//dot product of x&y
		multiply(temp_x, temp_y, temp_x);

		//sum of all element in the final mat
		Scalar tot_t = sum(temp_x);

		return (double)tot_t[0] / ((double)std_y[0] * std_x[0]) / ((double)mat_x.rows * mat_x.cols);
	}

	static float find_crop2(Mat big, Mat small) {
		Functions fn;
		int r1 = big.rows;		int c1 = big.cols;
		int r2 = small.rows;	int c2 = small.cols;
		if (r1 < r2 || c1 < c2) {
			cout << "small is not small\n";
			return -1;
		}
		clock_t t1 = 0;

		map<string, cv::Mat> map;
		//if (r1 >= 2 * r2 && c1 >= 2 * c2) {
		//	small.copyTo(small_tmp);
		//	copyMakeBorder(small, small_tmp, 0, r2, 0, c2, 0);

		//	map[to_string(small_tmp.rows) + "," + to_string(small_tmp.cols)] = small_tmp;
		//}

		
		for (auto i = 0; i <= r1-r2; i += r2) {
			for (auto j = 0; j <= c1-c2; j += c2) {

				int p = std::min(r1-r2-i, r2);
				int q = std::min(c1-c2-j, c2);
				string key = to_string(r2 + p) + "," + to_string(c2 + q);

				//cout << key << endl;

				Mat small_tmp;
				if (map.find(key)== map.end()) {
					//cout << "inside\n";
					cv::Mat tmp;
					small.copyTo(tmp);
					copyMakeBorder(small, tmp, 0, p, 0, q, 0);
					//small_used = small_tmp;
					//cout << tmp.rows << " " << tmp.cols << endl;
					map[key] = tmp;
				}
				auto t2 = clock();
				small_tmp = map.find(key)->second;


				//cout << "p,q:" << p << " " << q << endl; 
				//cout << " " << "rows: " << big(cv::Rect(j, i, c2 + q, r2 + p)).rows << " " << small_tmp.rows <<
				//	"cols: "<< big(cv::Rect(j, i, c2 + q, r2 + p)).cols << " " << small_tmp.cols << endl;

				Mat ccc= fn.crossCorrelation(big(cv::Rect(j, i, c2 + q, r2 + p)), small_tmp);
				//cout << " *********************\n\n";

				//cout << i << "," << j << ": " << small_tmp.rows << " " << small_tmp.cols << "  " << clock()-t2 << endl;
			}
		}
		cout << clock() - t1 << "sec\n";
		cout << "\n";
	}

	static double template_match_pce(Mat mat1, Mat mat2, Point &maxloc) {
		Functions fn;
		cv::Mat result;

		if ((mat1.rows > mat2.rows && mat1.cols < mat2.cols) ||
			(mat1.rows < mat2.rows && mat1.cols > mat2.cols)) {
			return -1;
		}

		matchTemplate(mat1, mat2, result, cv::TM_CCORR_NORMED);

		// best match
		double maxval;
		cv::minMaxLoc(result, NULL, &maxval, NULL, &maxloc);

		if (mat1.size().area() > mat2.size().area())
			return fn.pce_op_wo_image(mat1(cv::Rect(maxloc.x, maxloc.y, mat2.cols, mat2.rows)), mat2);
		else
			return fn.pce_op_wo_image(mat2(cv::Rect(maxloc.x, maxloc.y, mat1.cols, mat1.rows)), mat1);
	}


	//speed this up
	static double sd_corr(NoiseNode q_ns, DigestNode dig, int x_beg=0, int y_beg=0) {
		int nrow = q_ns.get_rows(), ncol = q_ns.get_cols();
		cout << nrow << " " << ncol << endl;

		int d_orow = dig.get_orows(), d_ocol = dig.get_ocols();
		int dx_sz = dig.get_rows(), dy_sz = dig.get_cols();

		int x_end = x_beg + nrow, y_end = y_beg + ncol;
		if (x_end > d_orow || y_end > d_ocol) return -2;
		
		//clock_t t1 = clock();

		// find alingments
		vector<float> query_px, dig_px;
		Mat x_ind, y_ind;  
		cv::divide(dig.get_indices(), d_ocol, x_ind);
		cv::subtract(dig.get_indices(), x_ind * d_ocol, y_ind);


		clock_t t1 = clock();
		for (int i = 0; i < x_ind.rows; i++) {
			for (int j = 0; j < x_ind.cols; j++) {
				int p = x_ind.at<int>(i, j);
				int q = y_ind.at<int>(i, j);
				if(inRange(p, q, x_beg, y_beg, x_end, y_end)){
					dig_px.push_back(dig.get_prnu().at<float>(i, j));
					query_px.push_back(q_ns.get_prnu().at<float>(p - x_beg, q - y_beg));
				}
			}
		}
		cout << "time is : " << clock() - t1 << endl;;

		cout << dig_px.size() << "  " << query_px.size() << endl;

		cv::Mat q_mat(1, (int)query_px.size(), CV_32F);
		memcpy(q_mat.data, &query_px[0], query_px.size());

		cv::Mat d_mat(1, (int)dig_px.size(), CV_32F);
		memcpy(d_mat.data, &dig_px[0], dig_px.size());
		//cout << "time is : " << clock() - t1 << " corr: "<< corr2(q_mat, d_mat)  <<endl;;

		cout << q_mat.at<float>(0, 19) << " " << d_mat.at<float>(0, 19) << endl;

		return corr2(q_mat, d_mat);
	}
	static bool inRange(int p, int q, int x_beg, int y_beg, int x_end, int y_end) {
		return p >= x_beg && q >= y_beg && p <= x_end && q <= y_end;
	}


	//cv::Mat createSmallFP(cv::Mat FP, cv::Mat sd_index) {
	//	int sz = sd_index.cols;
	//	Mat sdFP(sz, sz, CV_32F);
	//	for (int i = 0; i < sz; i++)
	//		for (int j = 0; j < sz; j++)
	//			sdFP.at<float>(i, j) = FP.at<float>(((int)sd_index.at<float>(i, j) % KILO),
	//			(int)(sd_index.at<float>(i, j) / KILO));
	//	return sdFP;
	//}

	//double small_corr2(Mat FP, Mat FE, Mat FE_index) {
	//	Mat smallFP = createSmallFP(FP, FE_index);

	//	return corr2(smallFP, FE);
	//}

	//int matching = 0;
	//void sd_corr(cv::Mat FP, vector<Mat> sd_vec, vector<Mat> sd_vec_index) {
	//	cout << " " << sd_vec.size() << endl;
	//	for (int i = 0; i < sd_vec.size(); i++) {
	//		double x = small_corr2(FP, sd_vec.at(i), sd_vec_index.at(i));

	//		if (x > 0.053) {
	//			matching++;
	//			//cout << images.at(i) << " " << x << endl;
	//			cout << i << "th fp is matching " << x << endl;
	//			return;
	//		}
	//	}
	//}

	//int sd_corr_pos(Mat FP, vector<Mat> sd_vec, vector<Mat> sd_vec_index) {
	//	cout << " " << sd_vec.size() << endl;
	//	double max = 0.050;
	//	//double max = 0.053;
	//	int max_pos = -1;
	//	for (int i = 0; i < sd_vec.size(); i++) {
	//		double x = small_corr2(FP, sd_vec.at(i), sd_vec_index.at(i));
	//		if (x > max && x <= 1) {
	//			//cout << images.at(i) << " " << x << endl;
	//			cout << i << "th fp is matching " << x << endl;
	//			max_pos = i;
	//			//return i;
	//		}
	//	}
	//	return max_pos;
	//}

};

#endif