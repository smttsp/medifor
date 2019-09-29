#ifndef TEST_FUNCTION_H
#define TEST_FUNCTION_H

#include <vector>
#include <numeric>      
#include <algorithm> 

#include "opencv\cv.h"
#include "opencv\highgui.h"

#include <opencv2\core\mat.hpp> 
#include <opencv2\core\core.hpp>
using namespace cv;
using namespace std;

#include <boost/filesystem.hpp>
#define bfs boost::filesystem 


class Func {

public:
	template <typename T>
	static vector<size_t> sort_indexes(const vector<T>& v) {
		vector<size_t> idx(v.size());
		iota(idx.begin(), idx.end(), 0);

		// sort indexes based on comparing values in v
		sort(idx.begin(), idx.end(),
			[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
		return idx;
	}

	static bool matIsEqual(const cv::Mat mat1, const cv::Mat mat2) {
		// treat two empty mat as identical as well
		if (mat1.empty() && mat2.empty()) {
			return true;
		}
		// if dimensionality of two mat is not identical, these two mat is not identical
		if (mat1.cols != mat2.cols || mat1.rows != mat2.rows || mat1.dims != mat2.dims) {
			return false;
		}
		cv::Mat diff;
		cv::compare(mat1, mat2, diff, cv::CMP_NE);
		int nz = cv::countNonZero(diff);
		return nz == 0;
	}

	static std::vector<std::string> boost_files_ina_directory(const std::string path, const std::string ext) {
		auto len = ext.length();
		vector <std::string> file_list;
		if (!path.empty())
		{

			bfs::path apk_path(path);
			bfs::recursive_directory_iterator end;

			for (bfs::recursive_directory_iterator i(apk_path); i != end; ++i)
			{
				const bfs::path cp = (*i);
				if (len == cp.extension().string().length() - 1) {
					if (cp.extension() == ("." + ext))
						file_list.push_back(cp.string());
				}
				else {
					if (cp.extension() == ext)
						file_list.push_back(cp.string());
				}
			}
		}
		return file_list;
	}

};

#endif