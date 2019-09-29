#ifndef WEINER_H
#define WEINER_H

#include "Functions.h"

class Wiener_cls {
public:
	void wienerInDFT(Mat& noiseClean, Mat init_noise, double sigma) {
		Functions fn;
		Mat noise;  init_noise.copyTo(noise);
		int height = noise.rows;
		int width = noise.cols;

		float noiseVar = (float)pow(sigma, 2.0);

		Mat planes[] = { Mat_<float>(noise), Mat::zeros(noise.size(), CV_32F) };

		Mat complexI;
		merge(planes, 2, complexI);						//4		// Add to the expanded another plane with zeros
		dft(complexI, complexI);						//72	// this way the result may fit in the source matrix
		split(complexI, planes);						//1		// planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))

		Mat plane0_double; Mat plane1_double;
		planes[0].convertTo(plane0_double, CV_32F);     //4
		planes[1].convertTo(plane1_double, CV_32F);		//4

		double mxn = sqrt((double)noise.rows * noise.cols);
		divide(plane0_double, mxn, plane0_double);		//7
		divide(plane1_double, mxn, plane1_double);		//7

		Mat noiseMag;	Mat noiseMag1;
		magnitude(plane0_double, plane1_double, noiseMag);   //6
		//get rid of plane0_double //plane1_double
		//clock_t beg = clock(); 
		fn.waveNoise(noiseMag, noiseMag1, noiseVar);			//1688 --> 1482 -->342 ye dustu
		//clock_t end = clock();

		divide(noiseMag1, noiseMag, noiseMag); //merge'e kadar 26 ms
		planes[0].convertTo(plane0_double, CV_32F); //return to initial planes[0]
		planes[1].convertTo(plane1_double, CV_32F); //return to initial planes[1]

		plane0_double = plane0_double.mul(noiseMag);
		plane1_double = plane1_double.mul(noiseMag);

		plane0_double.convertTo(planes[0], CV_32F); //return to initial planes[0]
		plane1_double.convertTo(planes[1], CV_32F); //return to initial planes[1]

		merge(planes, 2, complexI);

		cv::idft(complexI, noiseClean, DFT_SCALE | DFT_REAL_OUTPUT);	//39 ms

		noiseClean.convertTo(noiseClean, CV_32F);

		complexI.~Mat();
		plane0_double.~Mat();
		plane1_double.~Mat();
		planes[0].~Mat();
		planes[1].~Mat();
		//	delete[] planes;

	}
};
#endif