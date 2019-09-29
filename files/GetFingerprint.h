#ifndef GETFINGERPRINT_H
#define GETFINGERPRINT_H

/* Estimation of camera fingerprint consisting of Photo-Response Non-Uniformity (PRNU) 
*  noise and traces of dark current. Color images from one source camera are processed 
*  to output the camera fingerprint. The images filenames must be listed in Images 
*  structure. An example is included in the header comments of getFingerprint.m.
*/
/*
 INPUT:
    Images      list of color images to process (obtained through dir command)
                these images have to be from the same camera and the same size
    sigma       local std of extracted noise
 OUTPUT:
    RP          reference pattern - estimate of PRNU (in the output file)
    LP          linear pattern structured data
*/

#include "Functions.h"
#include "NoiseExtract.h"
#include "MakeONFilter.h"
#include "Wiener.h"
#include <iostream>
#include <fstream>
using namespace std;
int max_fp_elem = 100;

cv::Mat getFingerprint(vector<string> images, float sigma = 0, int outRP = 0, int outLP = 0, int imagesInRP = 0){
	Functions fn; NoiseExtract_cls ne; MakeONFilter_cls mak; Wiener_cls wie;
	auto database_size = images.size();
	if (database_size == 0) {
		cv::Mat FP2(1,1,CV_32F);
		cout << "no images" << endl;
		return FP2;
	}
	if(sigma < 1 || sigma > 5)
		sigma = 3;

	int waveletLevel = 4;  // this is by default 4
	int height, width, dim;

	cv::Mat qmf = mak.MakeONFilter("Daubechies",8);

	cv::Mat RPsum[3]; cv::Mat NN[3];
	//here we are going to do noise extraction for each element in the database
	int t = 0;
	for(auto i = 0; i < database_size; i++){
		cout << "*** image: " << i << " ***\r";
		// seeProgress(i)
		cv::Mat img = imread(images.at(i));
		if(img.rows < img.cols){
			//cout << "rotated " << endl;
			fn.rotate_CW(img);
		}
		
//		cout << img.rows << "  " << img.cols << "  ";
		if(img.rows == 0 || img.cols == 0)
			cout << "image size is 0 " << endl;

		if(t==0){//first image whose noise is gonna be extracted		
			height = img.rows;
			width = img.cols;
			dim = img.channels();
			if(dim != 3){
				cout << "first one is gray\n";
				continue;       //gray images are ignored
			}
			//initialize RPsum & NN
			for(int i = 0; i < 3; i++){
				RPsum[i] = cv::Mat::zeros(img.size(),CV_32F);
				NN[i] = cv::Mat::zeros(img.size(),CV_32F);
			}
			t++;
		}
		else{
			if(img.rows != height || img.cols != width){
				cout << " size not matched in " << images.at(i) << endl;
				continue;
			}
			if(dim != 3){
				cout << "gray images are ignored in " << images.at(i) << endl;
				continue;
			}
		}
		
		for(int i = 0; i < 3; i++){//1:3 olarak degistir
			cv::Mat noise1C;
			for(int ii = 0; ii < 3; ii++)
				noise1C.create(img.size(),CV_32F);
			cv::Mat im_splitted[3];
			split(img, im_splitted);

			//convert rgb to bgr
			ne.noiseExtract(noise1C, im_splitted[2-i], qmf, sigma, waveletLevel);
			cv::Mat inten(img.size(), CV_32F);

			fn.intenScale(inten,im_splitted[2-i]);
			cv::Mat satur = fn.saturation(im_splitted[2-i],"");

			multiply(inten, satur, inten);
			multiply(noise1C,inten,noise1C);
			add(RPsum[i],noise1C,RPsum[i]);

			pow(inten,2.0,inten);
			add(NN[i],inten,NN[i]);
		}
	}
	cout << endl;
	for(int i = 0; i < 3; i++){
		cv::add(NN[i],1,NN[i]);
		cv::divide(RPsum[i],NN[i],RPsum[i]);
		
		fn.zeroMeanTotal(RPsum[i],RPsum[i]);
		NN[i].~Mat();
//		cout << "RP[0][0] " << RPsum[i].at<float>(0,0) << endl;
	}
	
	cv::Mat RP(RPsum[0].size(), CV_32F);

	fn.mat_arr2gray(RPsum,RP);

	for(int i = 0; i < 3; i++)
		RPsum[i].~Mat();
	cv::Mat FP;
	cv::Scalar mean = 0, stdDev = 0;

	meanStdDev(RP, mean, stdDev);
	wie.wienerInDFT(FP,RP,(float)stdDev[0]);
	RP.~Mat();
	return FP;
}

#endif