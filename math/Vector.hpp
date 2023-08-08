#include <concepts>
#include <vector>
#include <concepts>

#include "basic_macros.h"
#include "functions.hpp"

BEGIN_MATH
BEGIN_VEC
using namespace functions;
template<floatingtype T>
class Vector
{
private:
	int size_{ 0 };
	std::vector<T> vec_{};
public:
	Vector() = default;
	Vector(int size) :size_(size) { vec_.resize(size);}
	Vector(int size, T val):size_(size) { vec_.resize(size,val);}
	Vector(std::vector<T> vec):vec_(vec),size_(vec.size()){}
	Vector(T* data, int n) :size_(n) { vec_.insert(vec_.begin(), data, data + n); }
	int size() const { return size_; }
	//拷贝重载,深拷贝
	Vector& operator=(const Vector& a)
	{
		vec_.resize(a.size());
		vec_.insert(vec_.begin(), vec_.end(), a.getVec());
		return *this;
	}
	T& operator()(int i)
	{
		return vec_[i];
	}
	void print()
	{
		std::cout << "----------------" << std::endl;
		for (int i = 0; i < size_; i++)
			std::cout << vec_[i] << std::endl;
		std::cout << "----------------" << std::endl;
	}
	std::vector<T> getVec() const { return vec_; }
	//模的平方
	T square() const { return functions::square(vec_); }
	//模
	T norm() const { return static_cast<T>(std::sqrt(square())); }
	Vector clone() const
	{
		return Vector(vec_);
	}
	
	//标准化
	Vector normalize() const
	{
		
		if (norm() < epsilon)
		{
			std::vector<T> result(size_);
			result.insert(result.begin(), vec_.begin(), vec_.end());
		}
		else
			return Vector(coef(vec_, T(1 / norm())));
	}
};
END_VEC
END_MATH



