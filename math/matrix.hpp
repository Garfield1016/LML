#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include "basic_macros.h"

BEGIN_MATH
using namespace functions;

static std::default_random_engine rand_engine;
static std::uniform_int_distribution<int> unifIntDis(-10, 10);
static std::uniform_real_distribution<float> unifFloatDis(1.0, 20.0);

template<floatingtype T>
class Matrix;

template<floatingtype T>
struct Eigres
{
	std::vector<T> eigenvalue;
	Matrix<T> eigenvector;
};

template<floatingtype T>
class Matrix 
{
private:
	std::vector<T> data_{};
	int rows_{ 0 };
	int cols_{ 0 };
	int size_{ 0 };

public:
	Matrix() = default;
	Matrix(int rows, int cols) :rows_(rows), cols_(cols), size_(rows * cols) { data_.resize(size_); }
	Matrix(int rows, int cols, T val) :rows_(rows), cols_(cols), size_(rows * cols) { data_.resize(size_, val); }
	Matrix(int rows,int cols, std::vector<T> vec) :data_(vec),rows_(rows),cols_(cols) { }
	Matrix(const Matrix& m)
	{
		
		data_ = m.data_;
		rows_ = m.rows_;
		cols_ = m.cols_;
	}
	
	Matrix(Matrix&& m) :data_(std::move(m.data_)), rows_(m.rows_), cols_(m.cols_) { m.rows_ = 0; m.cols_ = 0; }
	int rows() const { return rows_; }
	int cols() const { return cols_; }
public:
	Matrix operator=(Matrix& m1)
	{
		cols_ = m1.cols_;
		rows_ = m1.rows_;
		//data_.begin() = m1.data_.begin(); data_.end() = m1.data_.end();
		data_ = m1.data_;	
		return *this;
	}

	Matrix operator=(const Matrix&& m1)
	{
		cols_ = m1.cols_;
		rows_ = m1.rows_;
		data_=std::move(m1.data_);
		return *this;
	}
	
	Matrix operator-(const Matrix& m) const
	{
		if (cols_ != m.cols_ || rows_ != m.rows_)
			LML_ERROR("dimension error");
		Matrix result(rows_, cols_);
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
			{
				T val = at(i,j)-m.at(i,j);
				if (std::abs(val) < epsilon)
					val = 0;
				result.set(i, j, val);
			}
		}
		return result;
	}
	Matrix& operator-=(const Matrix& m) 
	{
		if (cols_ != m.cols_ || rows_ != m.rows_)
			LML_ERROR("dimension error");

		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
			{
				T val = at(i, j) - m.at(i, j);
				if (std::abs(val) < epsilon)
					val = 0;
				set(i, j, val);
			}
		}
		return *this;
	}
	Matrix operator+(const Matrix& m) const
	{
		if (cols_ != m.rows_ || rows_ != m.cols_)
			LML_ERROR("dimension error");
		Matrix result(cols_, cols_);
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
			{
				T val = at(i, j) + m.at(i, j);
				result.set(i, j, val);
			}
		}
		return result;
	}
	Matrix operator+=(const Matrix& m) const
	{
		if (cols_ != m.cols_ || rows_ != m.rows_)
			LML_ERROR("dimension error");

		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
			{
				T val = at(i, j) + m.at(i, j);
				set(i, j, val);
			}
		}
		return *this;
	}
	Matrix operator*(const Matrix& m) const
	{
		if (cols_ != m.rows_)
			LML_ERROR("dimension error");
		Matrix result(rows_, m.cols_);
		for (int i = 0; i < rows_; i++)
		{
			std::vector temp_row = std::move(getRow(i));
			for (int j = 0; j < m.cols_; j++)
			{
				std::vector temp_col = std::move(m.getCol(j));
				T val = functions::dot(temp_row, temp_col);
				if (std::abs(val) < epsilon)
					val = 0;
				result.set(i, j, val);
			}
		}
		//result.floorZero();
		return result;
	}
	Matrix operator*=(const Matrix& m) 
	{
		if (cols_ != m.rows_)
			LML_ERROR("dimension error");

		for (int i = 0; i < rows_; i++)
		{
			std::vector temp_row = std::move(getRow(i));
			for (int j = 0; j < m.cols_; j++)
			{
				std::vector temp_col = std::move(m.getCol(j));
				T val = functions::dot(temp_row, temp_col);
				if (std::abs(val) < epsilon)
					val = 0;
				set(i, j, val);
			}
		}
		//result.floorZero();
		return *this;
	}
	T& operator()(int i, int j)
	{
		if (i<0 || j<0 || i>rows_ - 1 || j>cols_ - 1)
			LML_ERROR("index error");
		int index = i * cols_ + j;
		return data_[index];
	}
public:
	Matrix clone() const
	{
		Matrix result;
		result.rows_ = rows_;
		result.cols_ = cols_;
		result.data_ = data_;
		return result;
	}
	void reshape(int rows,int cols)
	{
		if (data_.empty())
		{
			rows_ = rows;
			cols_ = cols;
			data_.resize(rows_ * cols_);
			return;
		}
		if (rows_ * cols_ != rows * cols)
			LML_ERROR("dimension error");
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
				set(i, j, data_[i]);
		}
	}
	T at(int i, int j) const
	{
		if (i<0 || j<0 || i>rows_ - 1 || j>cols_ - 1)
			LML_ERROR("index error");
		int index = i * cols_ + j;
		return data_[index];
	}
	void set(int i, int j, T val) 
	{
		if (i<0 || j<0 || i>rows_ - 1 || j>cols_ - 1)
			LML_ERROR("index error");
		int index = i * cols_ + j;
		data_[index] = val;
	}
	std::vector<T> getRow(int i) const
	{
		if (i<0 || i>rows_ - 1)
			LML_ERROR("index error");
		std::vector<T> result(cols_);
		std::copy(data_.begin() + i * cols_, data_.begin() + (i + 1) * cols_, result.begin());
		return result;
	}
	std::vector<T> getCol(int j) const
	{
		if (j<0 || j>cols_ - 1)
			LML_ERROR("index error");
		std::vector<T> result(rows_);
		for (int i = 0; i < rows_; i++)
		{
			int index = i * cols_ + j;
			result.at(i) = data_[index];
		}
		return result;
	}
	std::vector<T> getCol(int j,int begin,int end) const
	{
		if (j<0 || j>cols_ - 1)
			LML_ERROR("index error");
		std::vector<T> result(end-begin);
		for (int i = begin; i < end; i++)
		{
			int index = i-begin;
			result[index] = data_[i*cols_+j];
		}
		return result;
	}
	void setRow(int i, T val)
	{
		if (i<0 || i>rows_ - 1)
			LML_ERROR("index error");
		for (int j = 0; j< cols_; j++)
		{
			int index = i * cols_ + j;
			data_[index] = val;
		}	
	}
	void setRowCoeff(int i, T coeff)
	{
		if (i<0 || i>rows_ - 1)
			LML_ERROR("index error");
		for (int j = 0; j < cols_; j++)
		{
			int index = i * cols_ + j;
			data_[index] *=coeff;
		}
	}
	void setCol(int j, T val)
	{
		if (j<0 || j>cols_ - 1)
			LML_ERROR("index error");
		for (int i = 0; i < rows_; i++)
		{
			int index = i * cols_ + j;
			data_[index] = val;
		}
	}
	void setRow(int i, const std::vector<T>& arr)
	{
		if (i<0 || i>rows_ - 1)
			LML_ERROR("index error");
		if (arr.size() != cols_)
			LML_ERROR("dimension error");
		for (int j = 0; j < cols_; j++)
			set(i, j, arr[j]);
	}
	void setCol(int j, std::vector<T> arr)
	{
		if (j<0 || j>cols_ - 1)
			LML_ERROR("index error");
		if (arr.size() != rows_)
			LML_ERROR("dimension error");
		for (int i = 0; i < rows_; i++)
			set(i, j, arr[i]);
	}
	void swapRow(int i1, int i2)
	{
		if (i1 == i2)
			return;
		std::vector<T> temp(getRow(i1));
		setRow(i1, getRow(i2));
		setRow(i2, temp);
	}
	void swapRow(int i1, int i2, int begin, int end)
	{
		if (i1 == i2||begin>=end)
			return;

		std::vector<T> temp(end - begin);
		for (int j = 0; j < end - begin; j++)
			temp.at(j) = (*this).at(i1, begin+j);
		for (int j = 0; j < end - begin; j++)
		{
			(*this).set(i1, begin + j, (*this).at(i2, begin + j));
			(*this).set(i2, begin + j, temp[j]);
		}
	}
	void swapCol(int j1, int j2)
	{
		if (j1 == j2) return;
		std::vector<T> temp(getCol(j1));
		setCol(j1, getCol(j2));
		setCol(j2, temp);
	}
	Matrix transpose() const
	{
		Matrix result(cols_, rows_);
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
				result(j, i) = this->at(i, j);
		}
		return result;
	}
	void print() const
	{
		std::cout << "----------------" << std::endl;
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
			{
				std::cout << at(i, j)<<" ";
			}
			std::cout << std::endl;
		}
		std::cout << "----------------" << std::endl;
	}
	Matrix scale(T a)
	{
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
			{
				set(i, j, at(i, j) * a);
			}
		}
		return *this;
	}
public:
	static T randVal()
	{
		rand_engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
		return T(unifIntDis(rand_engine));
	}
	static Matrix randomMN(int rows, int cols)
	{
		Matrix result(rows, cols);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
				result.set(i, j, randVal());
		}
		return result;
	}
	static Matrix I(int n)
	{
		Matrix result(n, n, 0);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)
					result.set(i, j, T(1));
			}
		}
		return result;
	}
	static Matrix ones(int rows, int cols)
	{
		Matrix result(rows, cols);
		for (int i = 0; i < rows; i++)
		{
			result.setRow(i, 1);
		}
	}
	void setZero()
	{
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
				set(i, j, 0);
		}
	}
	void floorZero()
	{
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
			{
				if (std::abs(at(i, j)) < epsilon)
					set(i, j, 0);
			}
		}
	}
public:
	Matrix echlon() const
	{
		Matrix mat(rows_,cols_,0);
		if (mat.rows_ == 1||mat.cols_==1)
			return mat;
		//按照每行第一个非零元素的位置排序，按照行向量排序
		std::vector<int> temp_arr(rows_,int(-1));					//记录此行是否已经找过
		std::multimap<int, int> col_row;
		for (int j = 0; j < cols_; j++)
		{
			for (int i = 0; i < rows_; i++)
			{
				if (at(i, j) != 0 && temp_arr[i] == int(-1))	// &&col_row.size()<(i+1) && temp_arr[i] == int(-1)
				{
					temp_arr[i] = j; 
					col_row.insert(std::make_pair(j, i));
				}
			}
			if (*std::min_element(temp_arr.begin(),temp_arr.end()) != int(-1))
				break;
		}
		//遍历multimap，对mat进行排序
		int row = 0;
		for (auto it = col_row.begin(); it != col_row.end(); it++)
		{
			mat.setRow(row, this->getRow(it->second)); 
			row++;
		}
		//计算
		for (int i = 0; i < rows_-1; i++)
		{
			for (int i_level1 = i + 1; i_level1 < rows_; i_level1++)
			{
				if (mat.at(i, i) == 0)	//&&mat.at(i_level1,i)==0
					continue;
				T coeff = mat.at(i_level1, i)/mat.at(i, i);			//第一行 a(0,0),第二行a(1,1)...
				std::vector<T> temp_row = mat.getRow(i_level1)-functions::coef(mat.getRow(i), coeff);
				mat.setRow(i_level1, temp_row);
			}
		}
		return mat;
	}
	int rank() const
	{
		int r = 0;
		Matrix m = echlon();
		for (int i = 0; i < m.rows_; i++)
		{
			for (int j = 0; j < m.cols_; j++)
			{
				if (m.at(i, j) != 0)
				{
					r++;
					break;
				}
			}
		}
		return r;
	}
	std::vector<Matrix> LUP() const
	{
		Matrix L(std::move(I(rows_)));
		Matrix U(*this);
		Matrix P(std::move(I(rows_)));
		//U.print();
		for (int j = 0; j < cols_; j++)
		{
			//选主元
			int maxval_line = j;
			T maxval = std::abs(U.at(j, j));
			for (int k = j + 1; k < rows_; k++)
			{
				if (std::abs(U.at(k, j)) > maxval)
				{
					maxval = std::abs(U.at(k, j));
					maxval_line = k;
				}
			}
			if (maxval < epsilon)
			{
				std::cerr << " no decompose" << std::endl;
				throw "no decompose";
			}
			//L.print();
			L.swapRow(j, maxval_line,0,j);	//j,i,1:j-1
			//L.print();
			U.swapRow(j, maxval_line,j,cols_);	//j,i,j:n
			P.swapRow(j, maxval_line);	//j,i,:
			//end
			for (int i = j + 1; i < rows_; i++)
			{
				T coeff = U.at(i, j) / U.at(j, j);
				L.set(i, j, coeff);
				for (int k = j; k < cols_; k++)
				{
					T new_val = U.at(i, k) - coeff * U.at(j, k);
					U.set(i, k, new_val);
				}
				
			}
		}
		L.floorZero(); U.floorZero();
		std::vector<Matrix> lup(3);
		lup[0] = L;
		lup[1] = U;
		lup[2] = P;

		return lup;
	}
	Matrix inverse() const
	{
		if (rows_ != cols_)
			LML_ERROR("row != col");
		if (rank() != rows_)
			LML_ERROR("rank error");
		std::vector<Matrix<T>> lup = std::move(LUP());
		Matrix L(lup[0]);
		Matrix U(lup[1]);
		Matrix P(lup[2]);
		Matrix L_inv(I(rows_));
		Matrix U_inv(I(rows_));
		Matrix P_inv(I(rows_));
		//L_inv
		for (int j = 0; j < cols_; j++)
		{
			for (int i = j + 1; i < rows_; i++)
			{
				T coeff = L.at(i, j)/L.at(j,j);
				for (int k = 0; k < cols_; k++)
				{
					T new_val = L.at(i, k) - coeff * L.at(j,k);
					T new_val_inv = L_inv.at(i, k) - coeff * L_inv.at(j, k);
					L.set(i, k, new_val);
					L_inv.set(i, k, new_val_inv);
				}
			}
		}
		//U_inv
		for (int j=cols_-1; j >= 0; j--)
		{
			for (int i = j - 1; i >= 0; i--)
			{
				T coeff = U.at(i, j) / U.at(j, j);
				for (int k = cols_-1; k >= 0; k--)
				{
					T new_val = U.at(i, k) - coeff * U.at(j, k);
					T new_val_inv = U_inv.at(i, k) - coeff * U_inv.at(j, k);
					U.set(i, k, new_val);
					//L.print();
					U_inv.set(i, k, new_val_inv);
				}
			}

		}
		for (int i = 0; i < rows_; i++)
		{
			T coeff = U.at(i, i);
			for (int j = 0; j < cols_; j++)
			{
				T new_val_inv = U_inv.at(i, j) / coeff;
				U_inv.set(i, j, new_val_inv);
			}
		}
		//P_inv
		for (int j = 0; j < cols_; j++)
		{
			int maxval_line = j;
			for (int k = j + 1; k < rows_; k++)
			{
				if (P.at(k, j) == 1)
					maxval_line = k;
			}
			P.swapRow(j, maxval_line);
			P_inv.swapRow(j, maxval_line);
		}

		Matrix inv = U_inv*L_inv*lup[2];
		return inv;
	}
	
	std::vector<Matrix> QRbyHouseholder() const
	{
		std::vector<Matrix> result(2);
		Matrix Q = I(rows_);
		Matrix R = *this;
		Matrix H(rows_,rows_);
		std::vector<Matrix> H_arr;
		for (int j = 0; j < cols_-1; j++)
		{
			int H_size = rows_ - j;
			Matrix H_temp = I(H_size);
			std::vector<T> a1 = R.getCol(j,j,rows_);
			double norm_a1 = functions::norm(a1);
			std::vector<T> norma_a1(H_size);
			norma_a1[0]= functions::norm(a1);
			double base = functions::norm(a1 - norma_a1);
			std::vector<double> temp1 = coef(a1-norma_a1, 1 / base);
			Matrix u1(H_size, 1, std::move(temp1));
			Matrix u1t = u1.transpose();
			H_temp = I(H_size) -(u1 * u1t).scale(2);
			if (H_temp.rows_ != R.rows_)
			{
				int diff = R.rows_ - H_temp.rows_;
				for (int i_temp = 0; i_temp < H_temp.rows_; i_temp++)
				{
					for (int j_temp = 0; j_temp < H_temp.cols_; j_temp++)
						H(i_temp + diff, j_temp + diff)= H_temp.at(i_temp, j_temp);
				}
				for (int i_temp = 0; i_temp < diff; i_temp++)
					H(i_temp, i_temp)= 1;
			}
			else
				H = H_temp;
			
			Q=std::move(Q*H);
			R = std::move(H*R);
			std::vector<T> row_bottom = std::move(R.getRow(rows_ - 1));
			//极小值置为0
			T sumval = 0;
			for (int index = 0; index < row_bottom.size(); index++)
				sumval += std::abs(row_bottom[index]);
			if (sumval < epsilon)
				break;
			H.setZero();
		}
		result[0] = std::move(Q);
		result[1] = std::move(R);
		return result;
	}
	std::vector<Matrix> QRbyGivens() const
	{
		Matrix P=I(rows_);
		Matrix A(*this);
		for (int j = 0; j < cols_; j++)
		{
			for (int i = j+1; i <rows_; i++)
			{
				T c = A(j, j) / std::sqrt(A(j, j) * A(j, j) + A(i, j) * A(i, j));
				T s = A(i, j) / std::sqrt(A(j, j) * A(j, j) + A(i, j) * A(i, j));
				Matrix P_ = I(rows_);
				P_(j, j) = c;
				P_(j, i) = s;
				P_(i, j) = -s;
				P_(i, i) = c;
				P =std::move(P_*P);
				A = std::move(P_ * A);
			}
		}
		std::vector<Matrix> res(2);
		res[0] = std::move(P.transpose());
		res[1] = std::move(A);
		return res;
	}

	Eigres<T> eigen(int num) const
	{
		Matrix A(*this);
		Matrix Q(std::move(I(rows_)));
		Eigres<T> eigres;
		std::vector<T> e;
		std::vector<Matrix<double>> qr;
		for (int i = 0; i < num; i++)
		{
			qr = A.QRbyHouseholder();
			Q = Q * qr[0];
			A = qr[1] * qr[0];
		}
		Matrix Ak = qr[0] * qr[1];
		
		for (int i = 0; i < rows_; i++)
			e.emplace_back(Ak.at(i, i));
		for (int i = 0; i < rows_ - 1; i++)
		{
			T minval = e[i];
			for (int j = i + 1; j < rows_; j++)
			{
				if (e[i] < minval)
				{
					T temp = e[i];
					e[i] = e[j];
					e[j] = temp;
					Q.swapCol(i, j);
				}
			}
		}
		eigres.eigenvalue = std::move(e);
		eigres.eigenvector = std::move(Q);
		return eigres;
	}
	T determinant() const
	{
		T res=1;
		Matrix A(*this);
		for (int j = 0; j < cols_; j++)
		{
			//
			int maxval_line = j;
			T maxval = std::abs(A.at(j, j));
			for (int k = j + 1; k < rows_; k++)
			{
				if (std::abs(A.at(k, j)) > maxval)
				{
					maxval = std::abs(A.at(k, j));
					maxval_line = k;
				}
			}
			if (maxval < epsilon)
			{
				res = 0;
				return res;
			}
			if (maxval_line != j)
			{
				res *= -1;
				A.swapRow(j, maxval_line, j, cols_);	//j,i,j:n
			}
			for (int i = j + 1; i < rows_; i++)
			{
				T coeff = A.at(i, j) / A.at(j, j);
				for (int k = j; k < cols_; k++)
				{
					T new_val = A.at(i, k) - coeff * A.at(j, k);
					A.set(i, k, new_val);
					A(i, k)= new_val;
				}
			}
		}
		for (int i = 0; i < rows_; i++)
			res *= A.at(i, i);
		return res;
	}
	std::vector<Matrix<T>> SVD() const
	{
		Matrix AAT=(*this)*(this->transpose()); //U
		Matrix ATA = (this->transpose()) * (*this); //E
		Eigres<T> AAT_eigres = AAT.eigen(100);
		Eigres<T> ATA_eigres = ATA.eigen(100);

		Matrix U = std::move(AAT_eigres.eigenvector);
		Matrix V = std::move(ATA_eigres.eigenvector);
		Matrix E(rows_, cols_);
		std::vector<Matrix<T>> UEVt(3);
		for (int j = 0;j<cols_; j++)
		{
			Matrix Vi(V.rows_, 1, V.getCol(j));
			Matrix val1 = (*this) * Vi;	//A*Vi/Ui
			T avi = norm(((*this) * Vi).getCol(0));
			T ui = norm(U.getCol(j));
			E(j, j)= avi / ui;
		}
		Matrix VT = std::move(V.transpose());
		UEVt[0] = std::move(U);
		UEVt[1] = std::move(E);
		UEVt[2] = std::move(VT);
		return UEVt;
	} 
	std::vector<Matrix<T>> cholesky() const
	{
		Matrix A(*this);
		Matrix L(rows_, rows_);
		Matrix Lt(rows_, rows_);
		int n = rows_;

		for (int j = 0; j < cols_; j++)
		{
			L.set(j, j,T(std::sqrt(double(A.at(j, j)))));
			for (int i = j + 1; i < rows_; i++)
			{
				T val_l = T(1 / L.at(j, j)) * A.at(i, j);
				L(i, j)= val_l;
				L(j, i)= val_l;
			}	
			for (int i_sub = j + 1; i_sub < rows_; i_sub++)
			{
				for (int j_sub = j + 1; j_sub < cols_; j_sub++)
				{
					T A_new_i_j = A.at(i_sub, j_sub) - L.at(i_sub, j) * L.at(j, j_sub);
					A(i_sub, j_sub)= A_new_i_j;
				}
			}
		}
		//拆分开
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < cols_; j++)
			{
				if (i < j)
				{
					Lt(i, j)= L.at(i, j);
					L(i, j)= 0;
				}
				if(i==j)
					Lt(i, j)= L.at(i, j);
			}
		}
		std::vector<Matrix<T>> LLt(2);
		LLt[0] = std::move(L);
		LLt[1] = std::move(Lt);
		return LLt;
	}
};

template<typename T>
vec::Vector<T> sovleLU(const Matrix<T>& A, vec::Vector<T>& y)
{
	if (A.rows() != y.size())
		LML_ERROR("dimension error");
	vec::Vector<T> x(y.size());

	std::vector<Matrix<T>> lup = A.LUP();
	Matrix<T> PA = std::move(lup[2] * A);
	
	Matrix L=lup[0];
	Matrix U = lup[1];
	x(0) = y(0) / L(0, 0);
	for (int i = 1; i < L.rows(); i++)
	{
			T sumval = 0;
			for (int j_temp = 0; j_temp < i; j_temp++)
				sumval += x(j_temp) * L(i, j_temp);
			T val = (y(i) - sumval);
			x(i) = (y(i) - sumval) ;
	}
	vec::Vector<T> y2 = x;
	x(U.rows() - 1) = y2(U.rows() - 1) / U(U.rows() - 1, U.rows() - 1);
	T v= y2(U.rows() - 1) / U(U.rows() - 1, U.rows() - 1);
	for (int i = U.rows() - 2; i >= 0; i--)
	{
		T sumval = 0;
		for (int j_temp = U.cols()-1; j_temp > i; j_temp--)
			sumval += x(j_temp) * U(i, j_temp);
		T val = (y2(i) - sumval)/U(i,i);
		x(i) = (y2(i) - sumval)/U(i,i);
	}
	return x;
}



void t_general()
{
	using namespace functions;
	std::vector<double> a1{ 0,1,0,1 };
	std::vector<double> e1{ 1,0,0,0 };
	double norm_a1 = functions::norm(a1);
	double base = functions::norm(a1 - coef(e1, norm_a1));
	std::vector<double> u1 = coef((a1 - coef(e1, norm_a1)), 1 / base);
	int a = 1;
}

template<typename T>
class mat
{
private:
	std::shared_ptr<T[]> data_;
	int size_;
public:
	mat() = default;
	mat(int n):size_(n),data_(std::make_shared<T[]>(n)){}
	void set(int i, T val)
	{
		data_[i] = val;
	}
	T at(int i) { return data_[i]; }
};
 


END_MATH
