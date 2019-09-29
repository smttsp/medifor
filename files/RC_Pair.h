#pragma once

class RC_Pair {
private:
	int rows;
	int cols;
public:
	RC_Pair() : rows(-1), cols(-1) {}
	RC_Pair(int r, int c) { rows = r; cols = c;}
	int get_rows() { return rows; }
	int get_cols() { return cols; }
	void set_rows(int r) { rows = r; }
	void set_cols(int c) { cols = c; }

	//bool operator< (const RC_Pair& other) const {
	//	if (rows * cols == other.rows * other.cols)
	//		return rows < other.rows;
	//	else
	//		return rows * cols < other.rows * other.cols;
	//}

	bool operator< (const RC_Pair& other) const {
		if (rows == other.rows )
			return rows * cols < other.rows * other.cols;
		else
			return rows < other.rows;
	}

};