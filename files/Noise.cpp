#ifndef NOISE_CPP
#define NOISE_CPP
#include "Noise.h"

#include <boost/filesystem.hpp>
//using namespace boost::filesystem;
//-------------CONSTRUCTOR -----------
NoiseNode::NoiseNode(int rows, int cols, int type){
	Mat ns_mat= cv::Mat::zeros(rows, cols, type);
	this -> set_noise_node(ns_mat);
}
//Noise::Noise(int rows, int cols, int depth, int channels){
//	int type = CV_MAKE_TYPE(depth, channels);
//	cout <<"type "<< type << endl;
//	Mat ns_mat= cv::Mat::zeros(rows, cols, type);
//	this -> set_noise(ns_mat);
//}


//----------setters/getters -----------
void NoiseNode::set_noise_node(Mat new_ns){
	set_rows(new_ns.rows);
	set_cols(new_ns.cols);
	set_type(new_ns.type());
	set_prnu(new_ns);
}

//-----------methods-------------------
void NoiseNode::readNoise(string ns_file, Mat &ns){
	//if file not exist, catch exception
	FILE* fp = fopen(ns_file.c_str(),"rb");
	if (!fp) {
		perror("fopen");
		return;
	}
	//else get the values from the disk and create Noise object

	int header[3];
	fread(header,sizeof(float),3, fp);
	
	int rows = header[0];
	int cols = header[1];
	int type = header[2];

	//cout << rows << " " << cols << " " << type << endl;

	if(ns.rows != rows || ns.cols != cols)
		ns.create(rows, cols, type);

	float* buffer= new float[cols];
    for(int i=0;i<rows;++i){
        fread(buffer,sizeof(float),cols,fp);
        for(int j=0;j<cols;++j){
			ns.at<float>(i,j)=buffer[j];
        }
    }
    fclose(fp);
	delete[] buffer;
}

void NoiseNode::readNoise(string ns_file, NoiseNode& nn) {
	//if file not exist, catch exception
	FILE* fp = fopen(ns_file.c_str(), "rb");
	if (!fp) {
		perror("fopen");
		return;
	}
	//else get the values from the disk and create Noise object

	int header[3];
	fread(header, sizeof(int), 3, fp);

	int rows = header[0];
	int cols = header[1];
	int type = header[2];

	//cout << rows << " " << cols << " " << type << endl;

	Mat ns(rows, cols, type);

	float* buffer = new float[cols];
	for (int i = 0; i < rows; ++i) {
		fread(buffer, sizeof(float), cols, fp);
		for (int j = 0; j < cols; ++j) {
			ns.at<float>(i, j) = buffer[j];
		}
	}

	nn.set_noise_node(ns);

	fclose(fp);
	delete[] buffer;
}

void NoiseNode::writeNoise(string ns_file, NoiseNode nn) {
	int cols = nn.cols;
	int rows = nn.rows;

	FILE* fp = fopen(ns_file.c_str(), "wb");

	if (!fp) {
		cout << ns_file << " can't be created\n";
		return;
	}

	int header[3];
	header[0] = rows;
	header[1] = cols;
	header[2] = nn.type;
	//cout << rows << " " << cols << " " << ns.type() << endl;
	fwrite(header, sizeof(int), 3, fp);

	float* buffer = new float[cols];
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			buffer[j] = nn.get_prnu().at<float>(i, j);
		}
		fwrite(buffer, sizeof(float), cols, fp);
	}
	fclose(fp);
	delete[] buffer;
}

void NoiseNode::writeNoise(string ns_file, Mat ns){
	int cols = ns.cols;
	int rows = ns.rows;

	FILE* fp = fopen(ns_file.c_str(),"wb");
	
	if (!fp) {
		cout << ns_file << " can't be created\n";
		return;
	}

	int header[3];
	header[0] = rows;
	header[1] = cols;
	header[2] = ns.type();
	//cout << rows << " " << cols << " " << ns.type() << endl;
	fwrite(header, sizeof(int), 3, fp);

	float *buffer = new float[cols];
	for(int i=0;i<rows;++i){
		for(int j=0;j<cols;++j){
			buffer[j]=ns.at<float>(i,j);     
		}
		fwrite(buffer,sizeof(float),cols,fp);
	}
	fclose(fp);
	delete[] buffer;
}

void NoiseNode::printXbyY(Mat a, int x, int y){
	if (x > a.rows)
		x = a.rows;
	if (y > a.cols)
		y = a.cols;
	cout << "x,y "<< x  << " " << y << endl;
	cout << a(cv::Rect(0,0,y,x))<<endl;
}


RC_Pair NoiseNode::get_rc_from_ns_file(string ns_file) {
	RC_Pair rc = RC_Pair();

	FILE* fp = fopen(ns_file.c_str(), "rb");
	if (!fp) {
		perror("fopen");
		return rc;
	}
	//else get the values from the disk and create Noise object

	int header[3];
	fread(header, sizeof(int), 3, fp);
	uint64 tot_bytes = 3 * sizeof(int) + sizeof(float) * header[0] * header[1];
	uint64 x = boost::filesystem::file_size(ns_file.c_str());
	if (x != tot_bytes) {
		cout << "there is an error with the file size of " << ns_file << "\n";
	}

	fclose(fp);

	rc.set_rows(header[0]);
	rc.set_cols(header[1]);
	return rc;
}

#endif