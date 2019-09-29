#ifndef DIGEST_CPP
#define DIGEST_CPP

#include "Digest.h"
//#include "Noise.h"
//#include "Functions.h"

//-------------CONSTRUCTOR -----------
// I will fix this
//DigestNode::DigestNode() {
//	//Mat ns_mat = cv::Mat::zeros(rows, cols, type);
//	this->cols = 0;
//	this->rows = 0;
//	this->type = 0;
//}

//void fps2digVec() {
//	const std::string folder = "C:/Users/Samet/Desktop/main/fps/";
//	const std::string ext = "ndat";
//	std::vector<std::string> vec = boost_files_ina_directory(folder, ext);
//	cout << vec.size() << " noises" << endl;
//
//	vector<NoiseNode> noises;
//
//	for (auto file : vec) {
//		NoiseNode nn;
//		NoiseNode::readNoise(file, nn);
//		noises.push_back(nn);
//	}
//
//	for (auto nn : noises) {
//		cout << nn.get_rows() << " " << nn.get_cols() << endl;
//	}
//
//	HybridNode::groupNodes4Hybrid(noises, vec, 3);
//
//	cout << "---------\nAfter\n";
//	for (auto nn : noises) {
//		cout << nn.get_rows() << " " << nn.get_cols() << endl;
//	}
//}


void DigestNode::noise2digest(DigestNode &dig, const Mat ns, const double PERC) {
	
	if (PERC < 0 || PERC >= 100)
		return;

	int M = ns.size().area();
	int nth = (int)(M * (100.0 - PERC) / 100.0)-1;

	cv::Mat out;
	ns.copyTo(out);
	out = abs(out);

	std::vector<float> v(out.begin<float>(), out.end<float>()); // create a vector based on pointers
	std::nth_element(v.begin(), v.begin() + nth, v.end());

	threshold(out, out, v.at(nth), 1, THRESH_BINARY);
	multiply(ns, out, out);

	Mat digest (1, M - nth-1, CV_32F);
	Mat indices(1, M - nth-1, CV_32S);

	int cnt = 0;
	for (int i = 0; i < ns.rows; i++) {
		for (int j = 0; j < ns.cols; j++) {
			if (out.at<float>(i, j) != 0) {
				digest .at<float>((int)(cnt / digest.cols),  cnt % digest.cols)  = out.at<float>(i, j);
				indices.at<int>((int)(cnt / indices.cols), cnt % indices.cols) = (int)(i * ns.cols + j);
				cnt++;
			}
		}
	}
	//for (auto i = 0; i < 10; i++) {
	//	cout << digest.at<float>(indices.cols - i-1) << " " << indices.at<int>(indices.cols - i-1) << " \n";
	//}
	dig.set_rows(digest.rows);
	dig.set_cols(digest.cols);
	dig.set_orows(ns.rows);
	dig.set_ocols(ns.cols);
	dig.set_type(ns.type());
	dig.set_prnu(digest);
	dig.set_indices(indices);
}

//----------setters/getters -----------
void DigestNode::set_digest(DigestNode dig, const Mat ns, const Mat digest, const Mat indices) {
	dig.set_rows(digest.rows);
	dig.set_cols(digest.cols);
	dig.set_orows(ns.rows);
	dig.set_ocols(ns.cols);
	dig.set_prnu(digest);
	dig.set_indices(indices);
}

//-----------methods-------------------
void DigestNode::readDigest(const string dig_file, DigestNode &dig) {
	FILE* fp = fopen(dig_file.c_str(), "rb");
	if (!fp) {
		cout << "digest file " << dig_file << " not exist" << endl;
		return;
	}

	uint64 tot_bytes;
	fread(&tot_bytes, sizeof(uint64), 1, fp);

	//check the length of the file
	uint64 x = boost::filesystem::file_size(dig_file.c_str());
	if (x != tot_bytes) {
		cout << "there is an error with the file size of " << dig_file << "\n";
	}

	vector<int> header(5);
	fread(&header[0], sizeof(int), 5, fp);

	dig.set_rows(header[0]);
	dig.set_cols(header[1]);
	dig.set_orows(header[2]);
	dig.set_ocols(header[3]);
	dig.set_type(header[4]);

	Mat prnu = dig.get_prnu();
	Mat indices = dig.get_indices();

	int rows = dig.get_rows();
	int cols = dig.get_cols();

	if (prnu.rows != rows || prnu.cols != cols) {
		prnu.create(rows, cols, dig.get_type());
		indices.create(rows, cols, CV_32S);
	}

	float* buffer = new float[cols];
	for (int i = 0; i < rows; ++i) {
		fread(buffer, sizeof(float), cols, fp);
		for (int j = 0; j < cols; ++j) {
			prnu.at<float>(i, j) = buffer[j];
		}
	}

	int* buffer2 = new int[cols];
	for (int i = 0; i < rows; ++i) {
		fread(buffer2, sizeof(int), cols, fp);
		for (int j = 0; j < cols; ++j) {
			indices.at<int>(i, j) = buffer2[j]; //back inserter?
		}
	}
	dig.set_prnu(prnu);
	dig.set_indices(indices);

	fclose(fp);
	delete[] buffer;
	delete[] buffer2;
}

void DigestNode::readDigestVec(string vec_file, vector<DigestNode> &vec){
	FILE* fp = fopen(vec_file.c_str(), "rb");
	if (!fp) {
		cout << "digest file " << vec_file << " not exist" << endl;
		return;
	}

	uint64 tot_bytes;
	fread(&tot_bytes, sizeof(uint64), 1, fp);

	//check the length of the file
	uint64 x = boost::filesystem::file_size(vec_file);
	//cout << x << " " << tot_bytes << endl;

	uint64 fi_sz = (uint64)(sizeof(int) + sizeof(float));


	if (x != tot_bytes) {
		cout << "there is an error with the file size of " << vec_file << "\n";
	}
	uint64 read = sizeof(uint64);

	while (read < tot_bytes) {
		vector<int> header(5);
		fread(&header[0], sizeof(int), 5, fp);
		read += sizeof(int) * 5;

		DigestNode dig = DigestNode();
		dig.set_rows(header[0]);
		dig.set_cols(header[1]);
		dig.set_orows(header[2]);
		dig.set_ocols(header[3]);
		dig.set_type(header[4]);

		Mat prnu = dig.get_prnu();
		Mat indices = dig.get_indices();

		int rows = dig.get_rows();
		int cols = dig.get_cols();

		if (prnu.rows != rows || prnu.cols != cols) {
			prnu.create(rows, cols, dig.get_type());
			indices.create(rows, cols, CV_32S);
			rows = prnu.rows;
			cols = prnu.cols;
		}


		float* buffer = new float[cols];
		for (int i = 0; i < rows; ++i) {
			fread(buffer, sizeof(float), cols, fp);
			for (int j = 0; j < cols; ++j) {
				prnu.at<float>(i, j) = buffer[j];
			}
		}

		int* buffer2 = new int[cols];
		for (int i = 0; i < rows; ++i) {
			fread(buffer2, sizeof(int), cols, fp);
			for (int j = 0; j < cols; ++j) {
				indices.at<int>(i, j) = buffer2[j]; //back inserter?
			}
		}
		dig.set_prnu(prnu);
		dig.set_indices(indices);

		read += ((uint64)prnu.rows * prnu.cols)*fi_sz;

		vec.push_back(dig);
		delete[] buffer;
		delete[] buffer2;
	}

	fclose(fp);

}

vector<int> DigestNode::get_header(DigestNode dig) {
	int rows = dig.get_rows();
	int cols = dig.get_cols();
	int orows = dig.get_orows();
	int ocols = dig.get_ocols();
	int type = dig.get_type();
	
	vector<int> header = { rows, cols, orows, ocols, type };
	return header;
}

void DigestNode::writeDigest(string dig_file, DigestNode dig) {
	vector<int> header = get_header(dig);
	FILE* fp = fopen(dig_file.c_str(), "wb");

	if (!fp) {
		cout << "dig_file can't be created\n";
		perror("fopen");
		return;
	}
	uint64 tot_bytes = sizeof(uint64);

	//cout << " uint64 " << sizeof(uint64) << endl;
	uint64 fi_sz = (uint64)(sizeof(int) + sizeof(float));
	tot_bytes += (uint64)(5 * sizeof(int)) + ((uint64)dig.get_cols() * dig.get_rows()) * fi_sz;

	fwrite(&tot_bytes, sizeof(uint64), 1, fp);
	fwrite(&header[0], sizeof(int), header.size(), fp);

	int cols = dig.get_cols();
	int rows = dig.get_rows();

	float* buffer = new float[cols];
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			buffer[j] = dig.get_prnu().at<float>(i, j);
		}
		fwrite(buffer, sizeof(float), cols, fp);
	}

	int* buffer2 = new int[cols];
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			buffer2[j] = dig.get_indices().at<int>(i, j);
		}
		fwrite(buffer2, sizeof(int), cols, fp);
	}

	fclose(fp);
	delete[] buffer;
	delete[] buffer2;
}

void DigestNode::writeDigestVec(string vec_file, vector<DigestNode> vec){

	FILE* fp = fopen(vec_file.c_str(), "wb");
	if (!fp) {
		cout << "vec_file " << vec_file << " can't be created\n" ;
		perror("fopen");
		return;
	}

	uint64 tot_bytes = sizeof(uint64);
	uint64 fi_sz = (uint64)(sizeof(int) + sizeof(float));
	//cout << "bytes while writing: " << tot_bytes << endl;
	for (int i = 0; i < vec.size(); i++) {
		tot_bytes += (uint64)(5 * sizeof(int)) + ((uint64)vec[i].get_cols() * (uint64)vec[i].get_rows()) * fi_sz;
		//cout << i << "-th bytes while writing: " << tot_bytes << endl;
	}

	fwrite(&tot_bytes, sizeof(uint64), 1, fp);


	for (int i = 0; i < vec.size(); i++) {
		DigestNode dig = vec[i];
		vector<int> header = get_header(dig);
		fwrite(&header[0], sizeof(int), header.size(), fp);

		int cols = dig.get_cols();
		int rows = dig.get_rows();

		float* buffer = new float[cols];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				buffer[j] = dig.get_prnu().at<float>(i, j);
			}
			fwrite(buffer, sizeof(float), cols, fp);
		}

		int* buffer2 = new int[cols];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				buffer2[j] = dig.get_indices().at<int>(i, j);
			}
			fwrite(buffer2, sizeof(int), cols, fp);
		}
		delete[] buffer;
		delete[] buffer2;
	}
	fclose(fp);
}

#endif