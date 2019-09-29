#pragma once
//HYBRID STRUCTURE: It is like 2-level composite-digest search tree in this paper:
// FAST CAMERA FINGERPRINT MATCHING IN VERY LARGE DATABASES
#ifndef HYBRID_H
#define HYBRID_H

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "opencv\cv.h"
#include "opencv\highgui.h"

#include "Noise.h"
#include <opencv2\core\mat.hpp> 
#include <opencv2\core\core.hpp>
using namespace cv;
#include "Functions2.h"
#include <algorithm> 
#include "RC_Pair.h"
#include "Digest.h"

class HybridNode {
private:
	NoiseNode parent;
	vector<DigestNode> children;

public:

	HybridNode() { parent = NoiseNode();  }


	//-------------------------setters/getters----------------------
	NoiseNode get_parent() { return parent; }
	vector<DigestNode> get_children() { return children; }

	void set_parent(NoiseNode par) { parent = par; }
	void set_children(vector<DigestNode> kids) { children = kids; }

	void toString() {
		if (parent.get_rows() <= 0) {
			cout << "no parent";
		}
		else {
			cout << "parent size is:" << parent.get_rows() << " " << parent.get_cols() << endl;
			cout << "num children is: " << children.size() << " with first size: "
				<< children[0].get_rows() << " " << children[0].get_cols() << endl;
		}
	}

	//-----------------methods--------------------------------------
	static void writeHybrid(const string hyb_file, HybridNode hybrid);
	static void readHybrid(const string hyb_file, HybridNode& hybrid);
	static HybridNode extractHybrid(const vector<string> files, map<string, int> nfactor, int beg, int num);
	static HybridNode extractHybrid(const vector<NoiseNode> noises, int beg, int num);
	static void compute_HybridNodes(string hyb_fold, vector<string>& noise_files, map<string, int> nfactor, int num);
};




#endif