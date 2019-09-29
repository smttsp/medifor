/* Extraction of image noise residual from a color or grayscale image.
*/

#ifndef NOISEEXTRACTION_H
#define NOISEEXTRACTION_H

//#include "Functions.h"
#include "NoiseExtract.h"
#include "Functions.h"
#include "MakeONFilter.h"
#include "Wiener.h"
#include <sys\stat.h>

// *********** burdaki image nin type i = CV_8UC3 

/*
   function Noise = NoiseExtractFromImage(image,sigma,color,noZM) estimates PRNU from one image
  INPUT:
   image      test image filename
   sigma      std of noise to be used for identicication (recomended value between 2 and 3)
   color      if color then output Noise will be a colored noise
   noZM       if noZM then the output noise residual is not processed by zero-mean filter
  OUTPUT:
   Noise      extracted noise from the input image, a rough estimate of PRNU fingerprint
*/

//void readNoiseFloat(char* fpath,Mat& data);
//void saveNoiseFloat(string imageName, string fingerprintname);
//
//Mat getNoise(string imageName);
//Mat noiseExtractionFromImage(String imageName, float sigma, bool color, bool noZM);

class NoiseExtractionFromImage_cls {
public:
	int file_exist(const char* filename) {
		struct stat   buffer;
		return (stat(filename, &buffer) == 0);
	}

	Mat getNoise(string imageName) {
		Wiener_cls wie;
		if (!file_exist(imageName.c_str()))
			return Mat::zeros(1, 1, CV_32F);
		Mat noise;
		float sigma = 2;

		noise = noiseExtractionFromImage(imageName, sigma, false, false);
		Scalar mean = 0, stdDev = 0;

		meanStdDev(noise, mean, stdDev);
		Mat noisex;

		wie.wienerInDFT(noisex, noise, (float)stdDev.val[0]);
		noise.~Mat();			//*****dealloc*****	
		return noisex;
	}

	Mat noiseExtractionFromImage(String imageName, float sigma, bool color, bool noZM) {
		NoiseExtract_cls ne;
		MakeONFilter_cls mak;
		Functions fn;
		Mat qmf;
		Mat noise;
		Mat im_rgb;


		int waveletLevel = 4;		//L		// number of wavelet decomposition levels (between 2-5 as well)
		int img_height;				//M0	// height of the image  
		int img_width;				//N0	// width of the image
		int img_channels;			//three    //?????????????????? im not sure about that

		//temporary element
		Mat noise1C[3];
		try {
			im_rgb = imread(imageName, 1);

			//if(im_rgb.rows < KILO || im_rgb.cols < KILO){
			//	cout << "small image " << imageName << endl;
			//	Mat m(1,1,CV_32F);
			//	return m;
			//}

			if (im_rgb.rows < im_rgb.cols) {
				//cout << "rotated " << endl;
				fn.rotate_CW(im_rgb);
			}

			if (!im_rgb.data) {
				printf("!!! No data !!!");
				//			return;
			}
		}
		catch (exception& e) { e.what(); }

		//image exist here
		img_width = im_rgb.size().width;
		img_height = im_rgb.size().height;
		img_channels = im_rgb.channels();

		noise.create(im_rgb.size(), CV_32F);

		qmf = mak.MakeONFilter("Daubechies", 8);
		if (img_channels != 3) {
			ne.noiseExtract(noise, im_rgb, qmf, sigma, waveletLevel);
		}
		else {
			cv::Size sz = im_rgb.size();
			Mat im_splitted[3];
			cv::split(im_rgb, im_splitted);
			im_rgb.~Mat();
			//convert rgb to bgr
			for (int i = 0; i < img_channels; i++) {
				ne.noiseExtract(noise1C[2 - i], im_splitted[i], qmf, sigma, waveletLevel);
				im_splitted[i].~Mat();
			}

			if (!color) {
				noise.create(sz, CV_32F);
				fn.mat_arr2gray(noise1C, noise);
			}
			else {
				noise.create(sz, CV_32FC3);
				merge(noise1C, img_channels, noise);
			}
		}

		for (int i = 0; i < 3; i++)
			noise1C[i].~Mat();
		if (noZM)
			cout << "not removing the linear pattern";
		else
			fn.zeroMeanTotal(noise, noise);

		if (noise.type() != 5)
			noise.convertTo(noise, CV_32F);

		return noise;
	}
};


#endif 