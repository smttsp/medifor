#ifndef WAVELET_H
#define WAVELET_H

#include "Functions.h"
#include "MDWT.h"
#include "MIDWT.h"


class Wavelet_cls {
public:
	void imgToFloat(float* x, Mat y) {
		int nrow = y.rows;
		int ncol = y.cols;

		for (int i = 0; i < ncol; i++) {
			for (int j = 0; j < nrow; j++) {
				x[nrow * i + j] = (float)y.at<uchar>(j, i);
			}
		}
	}

	void d_matToFloat(float* x, Mat y) {
		int nrow = y.rows;
		int ncol = y.cols;

		for (int i = 0; i < ncol; i++) {
			for (int j = 0; j < nrow; j++) {
				x[nrow * i + j] = y.at<float>(j, i);
			}
		}
	}

	void mdwt_init(Mat& wave_trans, Mat image, Mat qmf, int wave_level) {

		int im_width = image.size().width;
		int im_height = image.size().height;

		int qmf_height, qmf_width, lh;
		qmf_width = qmf.size().width;
		qmf_height = qmf.size().height;

		if (qmf_width > qmf_height)
			lh = qmf_width;
		else
			lh = qmf_height;

		if (wave_level < 0)
			cerr << "The number of levels, L, must be a non-negative integer\n";
		if (wave_level == -1) { //estimate L(wavelevel)
		}

		if (im_height > 1) {
			float mtest = (float)(im_height / pow(2.0, (double)wave_level));
			if (!isint(mtest))
				cerr << "The matrix row dimension must be of size m*2^(L)";
		}
		if (im_width > 1) {
			float ntest = (float)(im_width / pow(2.0, (double)wave_level));
			if (!isint(ntest))
				cerr << "The matrix col dimension must be of size m*2^(L)";
		}
		//noise.rows*noise.cols
		int tsize = im_height * im_width;
		float* y; y = new float[tsize];
		float* x; x = new float[tsize];


		clock_t beg = clock();

		imgToFloat(x, image);

		clock_t end = clock();
		//	cout << " img to float " << end - beg << endl;
		fill_n(y, tsize, 0);

		float* d_qmf; d_qmf = new float[qmf_width * qmf_height];

		d_matToFloat(d_qmf, qmf);

		MDWT_cls::MDWT(x, im_height, im_width, d_qmf, lh, wave_level, y);

		clock_t beginn = clock();

		Mat wave_trans2(im_width, im_height, CV_32F, y);
		wave_trans = wave_trans2.t();

		delete[] x;
		delete[] y;
		delete[] d_qmf;
		clock_t endd = clock();
		//	cout << "copy : " << endd - beginn << endl;

	}

	void midwt_init(Mat& image_noise, Mat wave, Mat qmf, int wave_level) {

		int im_width = wave.size().width;
		int im_height = wave.size().height;

		int qmf_height, qmf_width, lh;
		qmf_width = qmf.size().width;
		qmf_height = qmf.size().height;

		if (qmf_width > qmf_height)
			lh = qmf_width;
		else
			lh = qmf_height;

		if (wave_level < 0) { //estimate L(wavelevel)
			int i = im_width; int j = 0;
			while (even(i)) {
				i = (i >> 1);
				j++;
			}
			wave_level = im_height; i = 0;
			while (even(wave_level)) {
				wave_level = (wave_level >> 1);
				i++;
			}
			if (mymin(im_height, im_width) == 1)
				wave_level = mymax(i, j);
			else
				wave_level = mymin(i, j);
			if (wave_level == 0) {
				cerr << "Maximum number of levels is zero; no decomposition can be performed!";
				return;
			}
			//		cout << wave_level << " yeni wavelevel\n"; 
		}

		if (im_height > 1) {
			float mtest = (float)(im_height / pow(2.0, (double)wave_level));
			if (!isint(mtest))
				cerr << "The matrix row dimension must be of size m*2^(L)";
		}
		if (im_width > 1) {
			float ntest = (float)(im_width / pow(2.0, (double)wave_level));
			if (!isint(ntest))
				cerr << "The matrix col dimension must be of size m*2^(L)";
		}
		int tsize = im_height * im_width;
		float* y; y = new float[tsize];
		float* x; x = new float[tsize];

		clock_t beg = clock();


		d_matToFloat(y, wave);

		//	cout << y[500]<< " y here "; 
		//	cout << y[1000];
		clock_t end = clock();

		//	cout << "mat to float " << end - beg << endl;

		fill_n(x, tsize, 0);

		float* d_qmf; d_qmf = new float[qmf_width * qmf_height];
		d_matToFloat(d_qmf, qmf);

		MIDWT_cls::MIDWT(y, im_height, im_width, d_qmf, lh, wave_level, x);

		int rem = 0, div = 0;
		for (int i = 0; i < tsize; i++) {
			rem = i % im_height;
			div = i / im_height;
			image_noise.at<float>(rem, div) = x[i];
		}

		wave.~Mat();
		delete[] x;
		delete[] y;
		delete[] d_qmf;
	}

};
#endif