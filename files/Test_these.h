#pragma once

#include "Digest.h"
#include "Hybrid.h"
#include "Correlation.h"
//#include "Functions.h"
#include "Functions2.h"
#include "NoiseExtractFromImage.h"

class Test_these {

public:
	void test_this() {
		const std::string folder = "C:/Users/Samet/Desktop/main/fps/";
		const std::string ext = "ndat";
		std::vector<std::string> vec = Func::boost_files_ina_directory(folder, ext);
		cout << vec.size() << " noises" << endl;

		vector<NoiseNode> noises;

		for (auto file : vec) {
			NoiseNode nn;
			NoiseNode::readNoise(file, nn);
			noises.push_back(nn);
		}

		for (auto nn : noises) {
			cout << nn.get_rows() << " " << nn.get_cols() << endl;
		}

		//HybridNode::compute_HybridNodes(noises, vec, 3);

		cout << "---------\nAfter\n";
		for (auto nn : noises) {
			cout << nn.get_rows() << " " << nn.get_cols() << endl;
		}
	}

	void test_this2() {
		string image = "C:/Users/Samet/Desktop/images/139999246.jpg";
		NoiseExtractionFromImage_cls nef;
		Mat ns = nef.getNoise(image);
		//cout << ns.size().area() << endl;

		DigestNode dg = DigestNode();
		DigestNode dg_b = DigestNode();
		DigestNode dg_c = DigestNode();


		DigestNode::noise2digest(dg, ns, 1.9);
		DigestNode::noise2digest(dg_b, ns, 2.6);
		DigestNode::noise2digest(dg_c, ns, 3.7);

		vector<DigestNode> vec1 = { dg, dg_b, dg_c };

		clock_t t1 = clock();
		DigestNode::writeDigestVec("digest_vecc.ndgv", vec1);
		cout << "writing " << clock() - t1 << endl;
		//DigestNode::writeDigest("digest.ndg", dg);

		clock_t t2 = clock();
		vector <DigestNode> vec2;
		DigestNode::readDigestVec("digest_vecc.ndgv", vec2);
		cout << "reading " << clock() - t2 << endl;
		//DigestNode dg2 = DigestNode();
		//DigestNode::readDigest("digest.ndg", dg2);

		//cout << dg.get_cols() << " " << dg2.get_cols() << endl;
		//cout << dg.get_rows() << " " << dg2.get_rows() << endl;
		//cout << dg.get_ocols() << " " << dg2.get_ocols() << endl;
		//cout << dg.get_orows() << " " << dg2.get_orows() << endl;
		//cout << dg.get_type() << " " << dg2.get_type() << endl;
		//cout << vec2.size() << "vec2 \n";
		//for(auto i = 0; i < 3; i++)
		//	cout << matIsEqual(vec1[i].get_prnu(), vec2[i].get_prnu()) << endl;

	}

	void test_sd_corr() {
		HybridNode hyb = HybridNode();
		HybridNode::readHybrid("hybrid0.nhyb", hyb);

		NoiseNode q = hyb.get_parent();
		//q		//this part for testing only
		q.set_noise_node(q.get_prnu()(cv::Rect(1, 1, 1000, 990)));

		Correlation::sd_corr(q, hyb.get_children()[4], 1, 1);

	}

	void test_NCC() {
		Mat m1(4500, 4500, CV_32F);
		randn(m1, Scalar(0), Scalar(1));

		Mat m2(1000, 1000, CV_32F);
		randn(m2, Scalar(0), Scalar(1));

		Correlation::find_crop2(m1, m2);
		//auto t1 = clock();
		//Mat m2; m1(cv::Rect(10, 10, 1000, 1000)).copyTo(m2);
		//cout << clock() - t1 << " ms\n";

		//m1.at<float>(0) = 1;
		//cout << m1.at<float>(0) << "  " << m2.at<float>(0) << endl;
		//Correlation::NCC_crop(m1, m2);

	}

	void template_matching_cv() {
		NoiseExtractionFromImage_cls nef;
		Functions fn;
		Mat img = nef.getNoise("D:/datasets/medifor/MFC19/i3.jpg");
		Mat templ = nef.getNoise("D:/datasets/medifor/MFC19/i2.jpg");

		//Mat1b img = imread("D:/datasets/medifor/MFC19/i1.jpg", cv::IMREAD_GRAYSCALE);
		//Mat1b templ = imread("D:/datasets/medifor/MFC19/i2.jpg", cv::IMREAD_GRAYSCALE);
		// compute match     
		auto t1 = clock();
		cv::Mat result;

		if ((img.rows > templ.rows && img.cols < templ.cols) ||
			(img.rows < templ.rows && img.cols > templ.cols)) {
			return;
		}
	
		matchTemplate(img, templ, result, cv::TM_CCORR_NORMED);
		cout << clock() - t1 << "time" << endl;
		// best match

		cout << "heyy";
		if (result.size().area() < 1) {
			return ;
		}
		Point maxloc;
		double maxval;
		cv::minMaxLoc(result, NULL, &maxval, NULL, &maxloc);
		// display result 

		cout << "result matrix: " << result.rows << " vs " << img.rows << endl;
		cout << "val: " << maxval << " ";
		cout << "match_location: " << maxloc.x << "  " << maxloc.y << endl;
		//cout << "pce: " << get_pce(result(cv::Rect(maxloc.x, maxloc.y, templ.cols, templ.rows))).peak_value << endl;
		//cout << img(cv::Rect(maxloc.x, maxloc.y, templ.cols, templ.rows)).rows << " " << templ.rows << endl;
		cout << "pce: " << fn.pce_op_wo_image(img(cv::Rect(maxloc.x, maxloc.y, templ.cols, templ.rows)), templ) << endl;;
		cout << clock() - t1 << "sec\n";

		//Mat3b res;
		//cvtColor(img, res, COLOR_GRAY2BGR);
		//rectangle(res, cv::Rect(maxloc.x, maxloc.y, templ.cols, templ.rows),cv::Scalar(0, 255, 0));
		//imshow("match", res);
		//waitKey();
	}

	//double template_match_pce(Mat mat1, Mat mat2) {

	//	cv::Mat result;
	//	matchTemplate(mat1, mat2, result, cv::TM_CCORR_NORMED);

	//	if ((mat1.rows > mat2.rows && mat1.cols < mat2.cols) ||
	//		(mat1.rows < mat2.rows && mat1.cols > mat2.cols)) {
	//		return -1;
	//	}

	//	// best match
	//	Point maxloc;
	//	double maxval;
	//	cv::minMaxLoc(result, NULL, &maxval, NULL, &maxloc);

	//	if (mat1.size().area() > mat2.size().area()) 
	//		return pce_op_wo_image(mat1(cv::Rect(maxloc.x, maxloc.y, mat2.cols, mat2.rows)), mat2);
	//	else 
	//		return pce_op_wo_image(mat2(cv::Rect(maxloc.x, maxloc.y, mat1.cols, mat1.rows)), mat1);
	//}

};