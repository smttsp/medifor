#ifndef MAIN_CPP
#define MAIN_CPP


#include "NoiseExtractFromImage.h"
#include "MakeONFilter.h"
//#include <windows.h>

//#include "GetFingerprint.h"
//#include "PCE.h"
//#include "Functions.h"
//#include <conio.h>
//#include <time.h>
//#include "IO_operations.h"
//#include "Noise.h"
//#include <iostream>

#include <map>
#include <set>
#include "db_create_main.h"
#include "Digest.h"
#include "Hybrid.h"
#include "Functions2.h"
void readFile(string filename);
void createFP(String in_dir, String out_dir, string run_dir, string username);
void create_prnus(map<string, string>& dir_params, map<string, string>& img_srcs, map<string, int>& srcs);

#include <opencv2/core/core.hpp>
#include "Test_these.h"
#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <string>
#include "perf_analysis_main.h"

void load_txts(map<string, string> dir_params, map<string, string>& img_srcs, map<string, int>& srcs, map<string, int>& nfactor);
void set_params(map<string, string>& dir_params);
void createHybrids(map<string, string>& dir_params, map<string, int>& nfactor);

int main(int argc, char** argv){
	//if(argc == 2)
	//	dbc_main(argv[1]);
	//else {
	//	string cam_dir = "d:/datasets/_pureset/";
	//	dbc_main(cam_dir);
	//}

	map<string, string> dir_params;
	set_params(dir_params);
	map<string, string> img_srcs;
	map<string, int> srcs;
	map<string, int> nfactor;
	cout << "param: " << dir_params["img_fold"] << endl;

	load_txts(dir_params, img_srcs, srcs, nfactor);
	//createHybrids(dir_params, nfactor);

	//create_prnus(dir_params, img_srcs, srcs);


	Perf_main pm;
	pm.test_composite(dir_params,img_srcs, srcs);

	getchar();
	return 1;
}

void createHybrids(map<string, string>& dir_params, map<string, int> &nfactor) {
	string fp_fold = dir_params["fp_fold"];
	string hyb_fold = dir_params["hyb_fold"];

	if (!bfs::is_directory(hyb_fold))
		bfs::create_directory(hyb_fold);

	std::vector<std::string> vec;
	vector<string> fp_files = Func::boost_files_ina_directory(fp_fold, "ndat");
	//fp_files.erase(fp_files.begin() + 10, fp_files.end());
	cout << "num fingerprints: " << fp_files.size() << endl;

	HybridNode::compute_HybridNodes(hyb_fold, fp_files, nfactor, 50);
}

void set_params(map<string, string> &dir_params) {
	dir_params["main_fold"]= "D:/datasets/medifor/DARPA/";
	dir_params["fp_fold"] = dir_params["main_fold"] + "fps_/"; 
	dir_params["hyb_fold"] = dir_params["main_fold"] + "hyb50_norm2/";
	dir_params["img_fold"] = dir_params["main_fold"] + "test/";
	dir_params["res_fold"] = dir_params["main_fold"] + "results/";
	dir_params["img_sources"] = dir_params["img_fold"] + "images.txt";
	dir_params["nfactors"] = dir_params["fp_fold"] + "nfactor.txt";

	if (!bfs::is_directory(dir_params["fp_fold"])) bfs::create_directory(dir_params["fp_fold"]);
	if (!bfs::is_directory(dir_params["hyb_fold"])) bfs::create_directory(dir_params["hyb_fold"]);
	if (!bfs::is_directory(dir_params["img_fold"])) bfs::create_directory(dir_params["img_fold"]);
	if (!bfs::is_directory(dir_params["res_fold"])) bfs::create_directory(dir_params["res_fold"]);
}

#include "AlphaNumSort.h"

void load_txts(map<string, string> dir_params, map<string, string>& img_srcs, map<string, int>& srcs, map<string,int> &nfactor) {

	vector<string> txts4map = Func::boost_files_ina_directory(dir_params["hyb_fold"], "txt");
	std::sort(txts4map.begin(), txts4map.end(), doj::alphanum_less<std::string>());

	string file4set = dir_params["img_sources"];
	string nfact_file = dir_params["nfactors"];

	ifstream file1(file4set);
	while (file1.is_open() && !file1.eof()) {
		string line;
		std::getline(file1, line);
		boost::trim(line);
		if (line == "") continue;
		auto p = line.find_last_of(" ");
		if (img_srcs.find(line.substr(p + 1)) != img_srcs.end())
			cout << line.substr(p + 1) << endl;
		img_srcs[line.substr(p + 1)] = line.substr(0, p);
		//cout << p << " --> img:"<< line.substr(p + 1) << "- src:" << line.substr(0,p) <<"-" << endl;
	}
	cout << img_srcs.size() << "\n";
	file1.close();

	int cnt = 0;
	for (string txt : txts4map) {

		ifstream file2(txt);
		while (!file2.eof()) {
			string line;
			std::getline(file2, line);
			boost::trim(line);
			if (line == "") continue;
			boost::replace_all(line, "\"", "");
			srcs[line]= cnt;
		}
		file2.close();
		cnt++;
	}
	cout << srcs.size() << "\n";

	//load nfactors into a map
	ifstream file3(nfact_file);
	while (file3.is_open() && !file3.eof()) {
		string line;
		std::getline(file3, line);
		boost::trim(line);
		if (line == "") continue;
		auto p = line.find_last_of(" ");
		string key = line.substr(0, p);
		int value = atoi(line.substr(p + 1).c_str());

		if (nfactor.find(key) != nfactor.end())
			cout << key << endl;
		nfactor[key] = value;
		//cout << p << " --> img:"<< line.substr(p + 1) << "- src:" << line.substr(0,p) <<"-" << endl;
	}
	file3.close();

	cout << "nfactors size: " << nfactor.size() << endl;
}

void createFP(String in_dir, String out_dir, string run_dir = "", string username = "") {
	Functions fn;
	NoiseNode nnode;
	in_dir = fn.sanitize_dir(in_dir);
	out_dir = fn.sanitize_dir(out_dir);

	String folder = out_dir.substr(0, out_dir.length() - 1);
	folder = folder.substr(folder.find_last_of("\\") + 1);

	vector<string> img_files = Func::boost_files_ina_directory(in_dir, "jpg");
	if (img_files.size() == 0)
		return;

	Mat FP;
	FP = getFingerprint(img_files);
	fn.printXbyY(FP, 3, 3);
	String fp_file = out_dir + folder + ".ndat";
	//writeNoiseFloat(fp_file, FP, FP.cols, FP.rows); 
	nnode.writeNoise(fp_file, FP);
	if (run_dir != "") {
		run_dir = fn.sanitize_dir(run_dir);
		ofstream run_res(run_dir + "run_" + username + ".txt");
		run_res << fp_file << endl;
		run_res.close();
	}
	img_files.clear();
}

void readFile(string filename) {
	string line;
	ifstream myfile(filename.c_str());
	if (myfile.is_open()) {
		while (getline(myfile, line)) {
			cout << line << '\n';
		}
		myfile.close();
	}
	else cout << "Unable to open file";
}

void used_codes(){

	//NoiseNode nn;
	//string fp_file1 = "D:/datasets/medifor/DARPA/fps_/100348233@N07;sony;nex-7;4000;6000.ndat";
	//string fp_file2 = "D:/datasets/medifor/DARPA/fps_/100780671@N04;sony;ilce-7;4000;6000.ndat";
	//string ns_file = "D:/datasets/medifor/DARPA/test/14126050916.jpg.ndat";
	//Mat fp1, fp2,ns1 ;
	//nn.readNoise(fp_file1, fp1);
	//nn.readNoise(fp_file2, fp2);
	//nn.readNoise(ns_file, ns1);
	//Functions fn;

	//cout << "r1 " << fn.pce_op_wo_image(fp1, ns1) << endl;
	//cout << "r2 " << fn.pce_op_wo_image(fp2, ns1) << endl;
	//cout << "r3 " << fn.pce_op_wo_image(fp1+fp2, ns1) << endl;


	//string hyb_file = "D:/datasets/medifor/DARPA/hyb4/hybrid0.nhyb";
	//HybridNode hn;
	//HybridNode::readHybrid(hyb_file, hn);

	//cout << "r4 " << fn.pce_op_wo_image(hn.get_parent().get_prnu(), ns1) << endl;


	// EXTRACTING FINGERPRINTS IN THE DATABASE
	//if(argc == 2)
	//	dbc_main_par(argv[1]);
	//else {
	//	string cam_dir = "d:/datasets/_pureset/";
	//	dbc_main_par(cam_dir);
	//}


	//create_single_fp();
	//if(argc == 2)
	//	createFP(argv[1], argv[1]);
	//else if(argc == 3)
	//	createFP(argv[1], argv[2]);
	//else if(argc == 4)
	//	createFP(argv[1], argv[2], argv[3]);
	//else if(argc == 5)
	//	createFP(argv[1], argv[2], argv[3], argv[4]);
	//else{
	//	if(argc < 3)
	//		cerr << "not enough argument\n";
	//	else 
	//		cerr << "too many arguments\n";
	//	cout << "usage:\n    extractFP.exe \"image directory\" \"output FP directory\"\n";
	//}
}
#endif
