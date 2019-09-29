// Peak-to-Correlation Energy ratio

/* Like normalized correlation, Peak-to-Correlation Energy ratio (PCE) is a measure 
*  of similarity for two discrete signals. PCE is especially suitable for 2-dimensional 
*  camera fingerprints because presence of hidden periodic patterns (a latent source 
*  of false identification) decreases PCE. This implementation involves a fast 
*  cross-correlation computation using FFT.
*/

#ifndef PCE_H
#define PCE_H
#include "Functions.h"

#include <math.h>
//
//Mat removeNeighbourhood(Mat cross_corr, cv::Point peak_loc, int squareSize);
//PCE_node pce(Mat cross_corr, int col_shift, int row_shift, int square_size);
//
//PCE_node pce(Mat cross_corr, int col_shift, int row_shift, int square_size){
//	PCE_node output;
//
//	int cols = cross_corr.cols;
//	int rows = cross_corr.rows;
//	if(square_size < 1)
//		square_size = 11;
//
//	if(row_shift>rows)
//		row_shift = rows - 1; 
//	else if(row_shift<0)
//		row_shift=0;
//
//	if(col_shift>cols)
//		col_shift = cols- 1; 
//	else if(col_shift<0)
//		col_shift=0;
//
//	if(cross_corr.at<double>(0,0) == 0){			//to make it more efficient
//		if(countNonZero(cross_corr) == 0){			//if the energy is of cross_correlation is zero
//			//out icin bi typedef yapmak lazim 
//			output.pce_val = 0;
//			output.peak_value = 1;
//			output.peakLocation = Point(0,0);
//			return output;
//		}
//	}
//
//	Mat c_inrange = cross_corr(cv::Rect( cols-col_shift-1, rows-row_shift-1, col_shift+1,  row_shift+1));
//
//	MatIterator_<double> fw = max_element(c_inrange.begin<double>(),c_inrange.end<double>());
//	double max_cc = *fw;			//peakheight
//	
//	double minVal = -10000000;
//	double maxVal = 10000000000;
//	cv::Point maxLoc;				//xpeak, ypeak burasi
//	cv::Point minLoc;
//	minMaxLoc(c_inrange,&minVal,&maxVal,&minLoc,&maxLoc);
////	int index_max = maxLoc.x * c_inrange.cols + maxLoc.y+1;
//
//	double peak_height = max_cc;
//	output.peakLocation = cv::Point(row_shift - maxLoc.x, col_shift - maxLoc.y);
//
//	Mat cor_peakless = removeNeighbourhood(cross_corr, maxLoc, square_size);
////	cout << output.peakLocation << endl;
//
//	double correl = cross_corr.at<double>(rows-1, cols-1);
//	pow(cor_peakless,2.0,cor_peakless);
//
//	double pce_energy ;	Scalar scl1, scl2;
//	meanStdDev(cor_peakless,scl1,scl2);
//	pce_energy = scl1[0];
////	cout << pce_energy << "pce_energ\n";
//	
//	output.peak_value = pow(peak_height,2.0)/pce_energy*(peak_height/abs(peak_height));
////	cout << output.peak_value << endl;
//
//	double tempval = peak_height/sqrt(pce_energy)/sqrt(2.0);
////	cout << tempval << endl;
////	output.pce_val = 1/2;
//	//Out.pvalue = 1/2*erfc(peakheight/sqrt(PCE_energy)/sqrt(2));     % under simplifying assumption that C are samples from Gaussian pdf 
//
////  [Out.P_FA,Out.log10P_FA] = FAfromPCE(Out.PCE,prod(shift_range+1));
//	return output;
//}
//
//Mat removeNeighbourhood(Mat cross_corr, cv::Point peak_loc, int square_size){
//	int nrows = cross_corr.rows;
//	int ncols = cross_corr.cols;
//	
//	int radius = (square_size-1)/2;
//	
//	cross_corr = shiftCol(cross_corr,radius-peak_loc.y);
//	cross_corr = shiftRow(cross_corr,radius-peak_loc.x);
//
//	Mat Y(nrows-square_size, square_size,CV_32F);
//	cross_corr(cv::Rect(0,square_size, square_size,nrows-square_size)).copyTo(Y(cv::Rect(0,0 ,square_size, nrows-square_size)) );
//	
//	transpose(Y,Y);
//	Y.resize(1,1);
//
//	cv::Mat transpose_ccorr = cross_corr.t();
//
//	transpose_ccorr = transpose_ccorr.reshape(1,1);
//	Mat maymat = transpose_ccorr(cv::Range(0,1),cv::Range(8250,transpose_ccorr.size().area()));
//
//	hconcat(Y,maymat,Y);
//
//	return Y;
//}
//
//void FAfromPCE(int pce,Mat mymat, double& p_fa, double& log10p_fa){
//
//
//}

#endif
