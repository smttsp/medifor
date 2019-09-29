#ifndef DB_CREATE_MAIN_H
#define DB_CREATE_MAIN_H

#include <boost/asio.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/chrono.hpp>
#include <boost/thread.hpp>

#include "NoiseExtractFromImage.h"
#include "MakeONFilter.h"
//#include <windows.h>

#include "GetFingerprint.h"
#include "PCE.h"
#include "Functions.h"
#include "Functions2.h"
#include <conio.h>
#include <time.h>
#include "Noise.h"

#include <opencv2/core/core.hpp>
#include <boost/range/iterator_range.hpp>
#include <iostream>
#include <mutex>

vector<string> create_single_fp(string fold, string fp_file, string im_fold, string outfile);
void count_elems_in_single_fp(string fold, string fp_file, string im_fold, string outfile);

boost::mutex output_mutex;
mutex mut;

void create_hybrid_db(string main_fold) {
	string db_fold = main_fold + "flickr_uae/";
	string im_fold = main_fold + "test/";
	string fp_fold = main_fold + "fps/";
	string hyb_fold = main_fold + "hyb/";
	//1) size of each noise and sort them.
}

void create_prnu_worker(string img_file, string prnu_file) {
	NoiseExtractionFromImage_cls nefi;
	Mat ns = nefi.getNoise(img_file);
	NoiseNode nn; nn.writeNoise(prnu_file, ns);
}
void create_prnus(map<string, string>& dir_params, map<string, string>& img_srcs, map<string, int>& srcs) {

	size_t number_of_threads = boost::thread::hardware_concurrency() - 2;
	cout << "there are " << number_of_threads << " threads\n";

	auto cnt = 0;
	boost::asio::io_service io_service;
	boost::thread_group threads; {
		boost::scoped_ptr< boost::asio::io_service::work > work(new boost::asio::io_service::work(io_service));
		
		for (size_t t = 0; t < number_of_threads; t++) {
			threads.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
		}
		
		for (auto img : img_srcs) {
			if (srcs.find(img.second) != srcs.end()) {
				string prnu_file = dir_params["img_fold"] + img.first + ".ndat";
				if (!bfs::is_regular_file(prnu_file)) {
					io_service.post(boost::bind(create_prnu_worker, dir_params["img_fold"] + img.first, prnu_file));
					cnt++;
					//cout << cnt << "-th image is denoised\n";
				}
			}
		}
	}
	cout << cnt << " images will be denoised \n";
	io_service.stop();
	threads.join_all();
}

static int current = 0;
void dbc_main_par(string main_fold) {
	//string db_fold = main_fold;
	string db_fold = main_fold + "flickr_uae/";
	string im_fold = main_fold + "test/";
	string fp_fold = main_fold + "fps/";

	if(bfs::is_directory(db_fold)) {
		auto cnt = std::count_if(
			bfs::directory_iterator(bfs::path(db_fold)),
			bfs::directory_iterator(),
			static_cast<bool(*)(const bfs::path&)>(bfs::is_directory));

		cout << "THERE ARE " << cnt << " USERS\n";
	}
	else {
		return;
	}
	if (!bfs::is_directory(im_fold)) bfs::create_directory(im_fold);
	if (!bfs::is_directory(fp_fold)) bfs::create_directory(fp_fold);

	string outfile = main_fold + "image_counts.txt";

	size_t number_of_threads = boost::thread::hardware_concurrency();
	//size_t number_of_threads = 10;

	boost::asio::io_service io_service;
	boost::thread_group threads;{
		boost::scoped_ptr< boost::asio::io_service::work > work(new boost::asio::io_service::work(io_service));
		for (size_t t = 0; t < number_of_threads; t++){
			threads.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
		}
		cout << threads.size() << " threads\n";

		bfs::path p(db_fold);
		if (is_directory(p)) {
			//std::cout << db_fold << " contains:\n";

			string fp_file = "";
			for (auto& user : boost::make_iterator_range(bfs::directory_iterator(p), {})) {
				//std::cout << "$$ user: " << user << "\n";
				for (auto& cam : boost::make_iterator_range(bfs::directory_iterator(user), {})) {
					if (is_directory(cam)) {
						//std::cout << "**** camera: " << cam.path().filename() << "\n";

						fp_file = fp_fold + user.path().filename().string() + ";" + cam.path().filename().string() + ".ndat";
						if (!bfs::exists(fp_file)) {
							io_service.post(boost::bind(count_elems_in_single_fp, cam.path().string(), fp_file, im_fold, outfile));
						}
						else {
							current++;
							//cout << current++ << " file exist: " << fp_file << endl;
						}
					}
				}
			}
		}

	}
	io_service.stop();
	threads.join_all();
}

vector<string> create_single_fp(string fold, string fp_file, string im_fold, string outfile) {
	NoiseExtractionFromImage_cls nef;
	Functions fn;
	mut.lock();
	cout << current++ <<" create fp: " << fp_file << "\n***************\n";
	mut.unlock();

	vector <string> img_files = Func::boost_files_ina_directory(fold, "jpg");
	double THRESH = 100;
	auto m = 5, M = 100, N = 10;
	//auto M = 50, N = 10;

	vector<int> match;
	vector<string> mfiles;
	auto sz = img_files.size();

	//cout << " images in the create single fp: "  << img_files.size() << endl;
	vector<string> test;
	vector<string> train;
	cv::Mat FP;
	if (sz < m) {
		return mfiles;
	}
	else if (sz <= M ) {
		copy(img_files.begin(), img_files.begin() + sz, std::back_inserter(train));
	}
	else {
		copy(img_files.begin(), img_files.begin() + M, std::back_inserter(train));
		auto x = mymin(M + N, (int)img_files.size());
		copy(img_files.begin() + M, img_files.begin() + x, std::back_inserter(test));
	}

	//cout << train.size() << " " << test.size() << endl;

	//for (auto p : train)
	//	cout << p << endl;
	//for (auto p : test)
	//	cout << p << endl;

	FP = getFingerprint(train);
	NoiseNode nnode;
	nnode.writeNoise(fp_file, FP);

	for (auto ii = 0; ii < test.size(); ii++) {
		NoiseNode nnode;
		Mat noise = nef.getNoise(test[ii]);

		if (!fn.isprnuOK(noise)) {
			cout << "ns patladi\n";
		}
		else {
			double pce = fn.pce_op_wo_image(FP, noise);

			if (pce > THRESH) {
				match.push_back(ii + 1);
			}
			else {
				fn.rotate_CCW(noise);
				fn.rotate_CCW(noise);
				pce = fn.pce_op_wo_image(FP, noise);
				if (pce > THRESH) {
					//cout << test[ii] << " is matching: " << pce << endl;
					match.push_back(-ii - 1);
				}
				//else				cout << ii << " is NOT MATCHING " << pce << endl;
			}
		}
	}

	for (auto a : match) {
		//cout << "matchings" << a << endl;

		if (a > 0) {
			string file = test[a - 1];
			string out = im_fold + bfs::path(file).filename().string();
			//cout << file << "\n" << out << endl;
			if (!bfs::is_regular_file(bfs::path(out)))
				bfs::copy_file(file, im_fold + bfs::path(file).filename().string());

			mfiles.push_back(out);
		}
		else {
			string file = test[-a - 1];
			string out = im_fold + bfs::path(file).filename().string();
			//cout << file << "\n" << out << endl;

			cv::Mat im = cv::imread(file);
			fn.rotate_CCW(im);
			fn.rotate_CCW(im);

			std::vector<int> compression_params;
			compression_params.push_back(CV_IMWRITE_JPEG_QUALITY);
			compression_params.push_back(100);
			if (!bfs::is_regular_file(bfs::path(out)))
				imwrite(out, im, compression_params);
			mfiles.push_back(out);
		}
	}

	mut.lock();
	ofstream file; file.open(outfile, std::ios_base::app);
	for (auto mfile : mfiles) {
		file << bfs::path(fp_file).filename().string() << " " << bfs::path(mfile).filename().string() << endl;
	}
	file.close();
	mut.unlock();

	return mfiles;
}

void count_elems_in_single_fp(string fold, string fp_file, string im_fold, string outfile) {
	NoiseExtractionFromImage_cls nef;
	Functions fn;
	mut.lock();
	cout << current++ << " count: " << fp_file << "\n***************\n";
	mut.unlock();

	vector <string> img_files = Func::boost_files_ina_directory(fold, "jpg");
	double THRESH = 100;
	auto m = 5, M = 100, N = 10;
	//auto M = 50, N = 10;

	vector<int> match;
	vector<string> mfiles;
	auto sz = img_files.size();

	//cout << " images in the create single fp: "  << img_files.size() << endl;
	vector<string> test;
	vector<string> train;
	cv::Mat FP;
	if (sz < m) {
		return;
	}
	else if (sz <= M) {
		copy(img_files.begin(), img_files.begin() + sz, std::back_inserter(train));
	}
	else {
		copy(img_files.begin(), img_files.begin() + M, std::back_inserter(train));
		auto x = mymin(M + N, (int)img_files.size());
		copy(img_files.begin() + M, img_files.begin() + x, std::back_inserter(test));
	}

	mut.lock();
	ofstream file; file.open(outfile, std::ios_base::app);
	file << bfs::path(fp_file).filename().string() << " " << train.size() << endl;
	file.close();
	mut.unlock();
}


void dbc_main(string main_fold) {
	//string db_fold = main_fold + "cameras/";
	//string db_fold = main_fold;
	string db_fold = main_fold + "flickr_uae/";
	string im_fold = main_fold + "test/";
	string outfile = main_fold + "images_count.txt";

	//ofstream outfile(main_fold + "images.txt");
	bfs::path p(db_fold);
	if (is_directory(p)) {
		//std::cout << db_fold << " contains:\n";

		string fp_file = "";
		for (auto& user : boost::make_iterator_range(bfs::directory_iterator(p), {})) {
			//std::cout << "$$ user: " << user << "\n";
			for (auto& cam : boost::make_iterator_range(bfs::directory_iterator(user), {})) {
				if (is_directory(cam)) {
					std::cout << "**** camera: " << cam.path().filename() << "\n";

					fp_file = main_fold + "fps/" + user.path().filename().string() + ";" + cam.path().filename().string() + ".ndat";

					//cout << fp_file << endl;
					count_elems_in_single_fp(cam.path().string(), fp_file, im_fold, outfile);
				}
			}
		}
	}
	//outfile.close();
}

#endif
