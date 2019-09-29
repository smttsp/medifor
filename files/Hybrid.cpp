#include "Hybrid.h"


void HybridNode::writeHybrid(const string hyb_file, HybridNode hybrid){
	
	FILE* fp = fopen(hyb_file.c_str(), "wb");

	if (!fp) {
		cout << hyb_file << " can't be created\n";
		return;
	}

	//1 calculate total size (8 Bytes: tot_size, 12 Bytes: noise header, R*C*(sizeof float) noise data)
	//for each digest 5*4 Bytes + digest size*4 bytes
	uint64 tot_bytes = sizeof(uint64);

	tot_bytes += 3 * sizeof(int) + sizeof(float)*hybrid.get_parent().get_rows() * hybrid.get_parent().get_cols();
	for (auto d : hybrid.get_children()) {
		tot_bytes += 5 * sizeof(int) + (sizeof(float) + sizeof(int)) * d.get_rows() * d.get_cols();
	}
	fwrite(&tot_bytes, sizeof(uint64), 1, fp);

	//2 write noise
	int rows = hybrid.get_parent().get_rows();
	int cols = hybrid.get_parent().get_cols();
	int type = hybrid.get_parent().get_type();
	vector<int> header = { rows, cols, type };

	fwrite(&header[0], sizeof(int), 3, fp);

	Mat prnu = hybrid.get_parent().get_prnu();
	float* buffer = new float[cols];
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			buffer[j] = prnu.at<float>(i, j);
		}
		fwrite(buffer, sizeof(float), cols, fp);
	}
	delete[] buffer; 

	//3 write vec<digests>
	for (int i = 0; i < hybrid.get_children().size(); i++) {
		DigestNode dig = hybrid.get_children()[i];
		vector<int> header = DigestNode::get_header(dig);
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


void HybridNode::readHybrid(const string hyb_file, HybridNode& hybrid){
	//if file not exist, catch exception
	FILE* fp = fopen(hyb_file.c_str(), "rb");
	if (!fp) {
		perror("fopen");
		return;
	}
	
	//else get the values from the disk and create Noise object
	uint64 tot_bytes;
	fread(&tot_bytes, sizeof(uint64), 1, fp);

	//check the length of the file
	uint64 x = boost::filesystem::file_size(hyb_file);
	if (x != tot_bytes) {
		cout << "there is an error with the file size of " << hyb_file << "\n";
		return;
	}

	int header[3];
	fread(header, sizeof(int), 3, fp);

	int rows = header[0];
	int cols = header[1];
	int type = header[2];

	uint64 read = sizeof(uint64) + 3 * sizeof(int);
	//cout << rows << " " << cols << " " << type << endl;


	// reading the parent
	NoiseNode nn;
	Mat ns(rows, cols, type);

	float* buffer = new float[cols];
	for (int i = 0; i < rows; ++i) {
		fread(buffer, sizeof(float), cols, fp);
		for (int j = 0; j < cols; ++j) {
			ns.at<float>(i, j) = buffer[j];
		}
	}
	delete[] buffer;
	read += sizeof(float) * rows * cols;
	nn.set_noise_node(ns);
	hybrid.set_parent(nn);


	// reading the children
	vector<DigestNode> children;
	uint64 fi_sz = (uint64)(sizeof(int) + sizeof(float));

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

		float* buffer1 = new float[cols];
		for (int i = 0; i < rows; ++i) {
			fread(buffer1, sizeof(float), cols, fp);
			for (int j = 0; j < cols; ++j) {
				prnu.at<float>(i, j) = buffer1[j];
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

		read += ((uint64)prnu.rows * prnu.cols) * fi_sz;

		children.push_back(dig);

		delete[] buffer1;
		delete[] buffer2;
	}
	hybrid.set_children(children);

	fclose(fp);
}



HybridNode HybridNode::extractHybrid(const vector<string> ns_files, map<string, int> nfactor, int beg = 0, int num=-1) {
	if (num == -1) num = (int)ns_files.size();

	vector<NoiseNode> noises;
	int numfp = 0;
	for (int i = beg; i < beg + num && i < ns_files.size(); i++) {
		NoiseNode nn;
		NoiseNode::readNoise(ns_files[i], nn);

		string nf = bfs::path(ns_files[i]).filename().string();
		//int x = nfactor[nf];
		float x = cv::mean(cv::abs(nn.get_prnu()))[0];
		//x = (x == 0) ? 1 : x;
		//cout << x << "\t";


		//cout << cv::mean(cv::abs(nn.get_prnu())) << endl;
		nn.set_prnu(nn.get_prnu() / x);
		noises.push_back(nn);
		double minval, maxval;
		minMaxIdx(nn.get_prnu(), &minval, &maxval);

		//cout << cv::mean(cv::abs(nn.get_prnu()))[0] << " " << minval << " " << maxval << endl;
		numfp++;
	}
	cout << "There are " << noises.size() << " fingerprints for hybrid extraction\n";
	return extractHybrid(noises, 0, num);

}


HybridNode HybridNode::extractHybrid(vector<NoiseNode> noises, int beg = 0, int num=-1) {
	if (num == -1) num = noises.size();

	int r = 0, c = 0; //determine row and column of the composite
	for (auto i = beg; i < noises.size() && i < num; i++) {
		r = std::max(r, noises[i].get_rows());
		c = std::max(c, noises[i].get_cols());
	}

	NoiseNode par;
	vector<DigestNode> children;

	Mat par_ns(r, c, CV_32F);
	Mat par_mask(r, c, CV_32F);
	
	for (auto i = beg; i < noises.size() && i < num; i++) {
		int r1 = noises[i].get_rows();
		int c1 = noises[i].get_cols();

		//cout << "mean: " << cv::mean(cv::abs(noises[i].get_prnu()))[0] << endl;

		DigestNode dn = DigestNode();
		cv::add(par_ns(cv::Rect(0, 0, c1, r1)), noises[i].get_prnu(), par_ns(cv::Rect(0, 0, c1, r1)));
		cv::add(par_mask(cv::Rect(0, 0, c1, r1)), 1, par_mask(cv::Rect(0, 0, c1, r1)));

		cout << "mean: " << cv::mean(cv::abs(par_ns))[0] << endl;

		DigestNode::noise2digest(dn, noises[i].get_prnu(), 1.4);
		children.push_back(dn);
	}
	cv::add(par_mask, 0.00000001, par_mask);
	cv::sqrt(par_mask, par_mask);
	cv::divide(par_ns, par_mask, par_ns);

	HybridNode hyb = HybridNode();
	par.set_noise_node(par_ns);
	hyb.set_parent(par);
	hyb.set_children(children);

	hyb.toString();
	return hyb;
}


void HybridNode::compute_HybridNodes(string hyb_fold, vector<string>& noise_files, map<string, int> nfactor, int num = 128) {
	auto len = noise_files.size();

	vector<RC_Pair> rc_vec;

	auto it = noise_files.begin();
	while (it != noise_files.end()){
		string ns_file = *it;
		RC_Pair rc = NoiseNode::get_rc_from_ns_file(ns_file);

		if (rc.get_rows() > 0) {
			rc_vec.push_back(rc);
			it++;
		}
		else {
			cout << " Remove the ns file, there is discrepency with :" << ns_file << endl;
			noise_files.erase(it);
		}
	}

	vector<size_t> sind = Func::sort_indexes(rc_vec);

	vector<string> new_ns_files;
	for (auto a : sind) {
		new_ns_files.push_back(noise_files[a]);
	}

	int cnt = 0;
	for (auto i = 0; i < new_ns_files.size(); i += num) {
		cout << i << "-th :" << endl;

		string hyb_file = hyb_fold + "hybrid" + to_string(cnt++) + ".nhyb";
		HybridNode hyb = extractHybrid(new_ns_files, nfactor, i, num);
		HybridNode::writeHybrid(hyb_file, hyb);

		ofstream txtfile(hyb_file + ".txt");
		for (auto p = i; p < new_ns_files.size() && p < i + num; p++) {
			txtfile << bfs::path(new_ns_files[p]).filename() << "\n";
		}
		txtfile.close();
		cout << endl << endl;
		if(cnt > 2)	
			break;
	}
}
