#ifndef NOISEEXTRACT_H
#define NOISEEXTRACT_H

#include "Functions.h"
#include "Wavelet_Transforms.h"

class NoiseExtract_cls {
public:
	void noiseExtract(Mat& noise, Mat input, Mat qmf, double sigma, int& waveletLevel) {
		Wavelet_cls wc; 
		int dataType = input.type();
		Functions fn;
		//	if(dataType is not double)
		//		convert it to double		

		int height = input.size().height; //H=M 
		int width = input.size().width;   //W=N

		int m = (int)pow(2, (double)waveletLevel);
		// use padding with mirrored image content 
		int minpad = 2;

		int nr = (int)ceil((double)(height + minpad) / m) * m;
		int nc = (int)ceil((double)(width + minpad) / m) * m;				// dimensions of the padded image (always pad 8 pixels or more)
		int pr = (int)ceil((double)(nr - height) / 2);					// number of padded rows on the top
		int prd = (int)floor((double)(nr - height) / 2);				    // number of padded rows at the bottom
		int pc = (int)ceil((double)(nc - width) / 2);					// number of padded columns on the left
		int pcr = (int)floor((double)(nc - width) / 2);					// number of padded columns on the right

		Mat padded;
		padded.create(Size(nc, nr), input.type());

		//1 ms
		//pad upside
		for (int i = 0; i < pr; i++)
			input(cv::Rect(0, pr - 1 - i, width, 1)).copyTo(padded(cv::Rect(pc, i, width, 1)));
		//pad downside
		for (int i = 0; i < prd; i++)
			input(cv::Rect(0, height - 1 - i, width, 1)).copyTo(padded(cv::Rect(pc, pr + height + i, width, 1)));
		//copy middle part
		input.copyTo(padded(cv::Rect(pc, pr, width, height)));
		//copy left side
		for (int i = 0; i < pc; i++)
			padded(cv::Rect(pc + i, 0, 1, nr)).copyTo(padded(cv::Rect(pc - 1 - i, 0, 1, nr)));
		//copy right side
		for (int i = 0; i < pcr; i++)
			padded(cv::Rect(pc + width - i - 1, 0, 1, nr)).copyTo(padded(cv::Rect(pc + width + i, 0, 1, nr)));


		// Precompute noise variance and initialize the output 
		float noiseVar = (float)(sigma * sigma);


		Mat wave_trans = Mat::zeros(Size(nc, nr), CV_32FC1);

		clock_t beg_n = clock();
		wc.mdwt_init(wave_trans, padded, qmf, waveletLevel);
		padded.~Mat();

		clock_t end_n = clock();
		//	cout << " mdwt init " << end_n - beg_n << endl;

		clock_t begin = clock();

		for (int i = 0; i < waveletLevel; i++) {
			int HHigh = nc / 2;
			int HLow = 0;
			int VLow = 0;
			int VHigh = nr / 2;

			fn.waveNoise(wave_trans(cv::Rect(HHigh, VLow, nc / 2, nr / 2)), wave_trans(cv::Rect(HHigh, VLow, nc / 2, nr / 2)), noiseVar);
			fn.waveNoise(wave_trans(cv::Rect(HLow, VHigh, nc / 2, nr / 2)), wave_trans(cv::Rect(HLow, VHigh, nc / 2, nr / 2)), noiseVar);
			fn.waveNoise(wave_trans(cv::Rect(HHigh, VHigh, nc / 2, nr / 2)), wave_trans(cv::Rect(HHigh, VHigh, nc / 2, nr / 2)), noiseVar);

			nr /= 2;
			nc /= 2;
		}

		clock_t end = clock();
		//	cout << "time for wavelet " << end - begin << endl;
		for (int i = 0; i < nr; i++) {
			for (int j = 0; j < nc; j++)
				wave_trans.at<float>(i, j) = 0;
		}

		int nr1 = wave_trans.size().height;
		int nc1 = wave_trans.size().width;

		Mat image_noise = Mat::zeros(Size(nc1, nr1), CV_32FC1);

		wc.midwt_init(image_noise, wave_trans, qmf, waveletLevel);

		noise = image_noise(cv::Rect(pc, pr, width, height));

		wave_trans.~Mat();
		qmf.~Mat();
		image_noise.~Mat();
	}
};
#endif