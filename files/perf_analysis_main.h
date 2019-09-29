#pragma once
#ifndef PERF_ANALYSIS_MAIN_H
#define PERF_ANALYSIS_MAIN_H


#include <opencv2/core/core.hpp>
#include "Test_these.h"
#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <string>
#include <map>
#include <set>
#include "AlphaNumSort.h"

#include <mutex>

std::mutex mut2;
class Perf_main {
public:
	//void search_ns_ori_in_comp(map<string, string> dir_params, vector<HybridNode> hyb_db, const string ns_file,
	//	map<string, string> img_srcs, map<string, int> srcs);

	void test_composite(map<string, string>& dir_params, map<string, string>& img_srcs, map<string, int>& srcs) {
		NoiseNode nn;
		string img_fold = dir_params["img_fold"];
		string hyb_fold = dir_params["hyb_fold"];

		vector<string> ns_files = Func::boost_files_ina_directory(img_fold, "ndat");
		vector<string> hyb_files = Func::boost_files_ina_directory(hyb_fold, "nhyb");
		std::sort(hyb_files.begin(), hyb_files.end(), doj::alphanum_less<std::string>());

		vector<HybridNode> hyb_db;
		for (auto h : hyb_files) {
			HybridNode hyb;
			HybridNode::readHybrid(h, hyb);
			hyb_db.push_back(hyb);
		}

		//size_t number_of_threads = boost::thread::hardware_concurrency() - 2;
		//cout << "there are " << number_of_threads << " threads\n";

		auto cnt = 0;
		//boost::asio::io_service io_service;
		//boost::thread_group threads; {
		//	boost::scoped_ptr< boost::asio::io_service::work > work(new boost::asio::io_service::work(io_service));
		//	for (size_t t = 0; t < number_of_threads; t++) {
		//		threads.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
		//	}
			for (auto ns_file : ns_files) {
				//Perf_main pm;
				//cout << cnt++ << "-th: " << ns_file << endl;
				//void (Perf_main:: * fn)(int) = &Perf_main::search_ns_ori_in_comp;
				//io_service.post(boost::bind(pm.search_ns_ori_in_comp, dir_params, hyb_db, ns_file, img_srcs, srcs));
				search_ns_ori_in_comp(dir_params, hyb_db, ns_file, img_srcs, srcs);

			}
		//}
	//	io_service.stop();
	//	threads.join_all();
	}


	void search_ns_ori_in_comp(map<string, string> dir_params, vector<HybridNode> hyb_db, const string ns_file,
		map<string, string> img_srcs, map<string, int> srcs) {

		NoiseNode nn; Functions fn;
		vector<vector<double>> match;
		vector<vector<double>> mismatch;

		string img_file = bfs::path(ns_file).filename().string();
		img_file = img_file.substr(0, img_file.size() - 5);


		if (img_srcs.count(img_file) && srcs.count(img_srcs[img_file])) {
			int pos = srcs[img_srcs[img_file]];
			//cout << img_file << " " << pos;
			string fp_file = dir_params["fp_fold"] + img_srcs[img_file];
			Mat ns; nn.readNoise(ns_file, ns);
			Mat fp; nn.readNoise(fp_file, fp);
			Mat hyb = hyb_db[pos].get_parent().get_prnu();

			//int r1 = rand() % (ns.rows - DROW), c1 = rand() % (ns.cols - DCOL);
			cv::Rect roi(0, 0, ns.cols, ns.rows);

			Mat hyb_roi = hyb_db[pos].get_parent().get_prnu()(roi);

			double pce1 = fn.pce_op_wo_image(ns, fp);
			double pce2 = fn.pce_op_wo_image(ns, hyb_roi);

			//cout << "pces: " << pce1 << " vs. " << pce2 << endl;

			mut2.lock();
			ofstream file(dir_params["res_fold"] + "match_comp50_norm.txt", ios::app);
			match.push_back(vector<double> {pce1, pce2});
			file << img_file << "\t" << pos << "\t" << hyb.rows << "\t" << hyb.cols << "\t";
			for (double a : match.at(match.size() - 1)) {
				cout << a << " ";
				file << a << "\t";
			}
			cout << endl;
			file << "\n";
			file.flush();
			file.close();
			mut2.unlock();

			//mismatch.push_back(vector<double>());
			//for (int p = 0; p < hyb_db.size(); p++) {
			//	if (p == pos) continue;
			//	NoiseNode par = hyb_db[p].get_parent();
			//	if (par.get_rows() < ns.rows || par.get_cols() < ns.cols) {
			//		mismatch.at(mismatch.size() - 1).push_back(-100);
			//	}
			//	else {
			//		double pce2 = fn.pce_op_wo_image(ns, par.get_prnu()(roi));
			//		mismatch.at(mismatch.size() - 1).push_back(pce2);
			//	}
			//}

			//for (double a : mismatch.at(mismatch.size() - 1))
			//	cout << a << " ";
			//cout << endl;
		}
	}

	void search_ns_in_comp(map<string, string> dir_params, vector<HybridNode> hyb_db, const string ns_file,
		map<string, string>& img_srcs, map<string, int>& srcs) {

		NoiseNode nn; Functions fn;
		vector<int> match;

		int DROW = 1024, DCOL = 1024;

		string img_file = bfs::path(ns_file).filename().string();
		img_file = img_file.substr(0, img_file.size() - 5);

		if (img_srcs.count(img_file) && srcs.count(img_srcs[img_file])) {
			int pos = srcs[img_srcs[img_file]];
			cout << img_file << " " << pos;
			string fp_file = dir_params["fp_fold"] + img_srcs[img_file];
			Mat ns; nn.readNoise(ns_file, ns);
			Mat fp; nn.readNoise(fp_file, fp);
			Mat hyb = hyb_db[pos].get_parent().get_prnu();

			int r1 = rand() % (ns.rows - DROW), c1 = rand() % (ns.cols - DCOL);
			cv::Rect roi(c1, r1, DCOL, DROW);
			Mat ns_roi = ns(roi);
			Mat fp_roi = fp(roi);

			double minVal;
			double maxVal;
			Point minLoc;
			Point maxLoc;

			minMaxLoc(fp, &minVal, &maxVal, &minLoc, &maxLoc);


			cout << cv::mean(fp)[0] << " " << minVal << " " << maxVal << " " << endl;

			Mat hyb_roi = hyb_db[pos].get_parent().get_prnu()(roi);

			//double p1 = fn.pce_op_wo_image(ns_roi, fp_roi);
			//double p2 = fn.pce_op_wo_image(ns_roi, hyb_roi);
			double p1 = fn.pce_op_wo_image(ns, fp);
			double p2 = fn.pce_op_wo_image(ns, hyb(cv::Rect(0, 0, ns.cols, ns.rows)));


			cout << "pces: " << p1 << " vs. " << p2 << endl;
		}
		//for (auto hyb : hyb_db) {
		//	Point maxloc;
		//	if (hyb.get_parent().get_rows() < ns.rows || hyb.get_parent().get_cols() < ns.cols) continue;

		//	double pce = Correlation::template_match_pce(hyb.get_parent().get_prnu(), ns, maxloc);

		//	if (pce > 45) {
		//		cout << "pce is: " << pce << " its position: " << maxloc.y << ", " << maxloc.x << endl;
		//	}
		//	if (pce > 70)
		//		break;

		//}
	}
};

#endif