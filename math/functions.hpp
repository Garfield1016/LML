#pragma once
#include <vector>

#include <omp.h>
#include <concepts>
#include "basic_macros.h"

BEGIN_MATH
BEGIN_FUNCTIONS
void LML_ERROR(const std::string& message)
{
	std::cerr << message << std::endl;
	throw message;
}

//+
template<typename T>//, std::constructible_from<T> TF>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
	if (a.empty() || b.empty())
		throw "vector empty";
	if (a.size() != b.size())
		throw "vector size not equal";
	std::vector<T> result(a.size());
//#pragma omp parallel num_threads(4)
//		hello();
//	omp_set_num_threads(2);
//#pragma omp parallel for
	for (int i = 0; i < a.size(); i++)
		result[i] = a[i] + b[i];
	return result;
}
//-
template<typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
	static_assert(std::is_arithmetic<T>::value==true,"type error");
	if (a.empty() || b.empty())
		throw "vector empty";
	if (a.size() != b.size())
		throw "vector size not equal";
	std::vector<T> result(a.size());

	for (int i = 0; i < a.size(); i++)
		result[i] = a[i] - b[i];
	return result;
}
template<typename T>
void operator-=(std::vector<T>& a, const std::vector<T>& b)
{
	static_assert(std::is_arithmetic<T>::value == true, "type error");
	if (a.empty() || b.empty())
		throw "vector empty";
	if (a.size() != b.size())
		throw "vector size not equal";
	for (int i = 0; i < a.size(); i++)
		a[i] -= b[i];
}
//叉乘
template<typename T>
std::vector<T> cross(const std::vector<T>& a, const std::vector<T>& b)
{
	static_assert(std::is_arithmetic<T>::value == true, "type error");
	if (a.empty() || b.empty())
		throw "vector empty";
	if (a.size() != b.size())
		throw "vector size not equal";
	std::vector<T> result(a.size());

	for (int i = 0; i < a.size(); i++)
		result[i] = a[i] * b[i];
	return result;
}
//除
template<typename T>
std::vector<T> operator/(const std::vector<T>& a, const std::vector<T>& b)
{
	static_assert(std::is_arithmetic<T>::value == true, "type error");
	if (a.empty() || b.empty())
		throw "vector empty";
	if (a.size() != b.size())
		throw "vector size not equal";
	std::vector<T> result(a.size());
	for (int i = 0; i < a.size(); i++)
		result[i] = a[i] / b[i];
	return result;
}
//点乘
template<typename T>
T dot(const std::vector<T>& a, const std::vector<T>& b)
{
	static_assert(std::is_arithmetic<T>::value == true, "type error");
	if (a.empty() || b.empty())
		throw "vector empty";
	if (a.size() != b.size())
		throw "vector size not equal";
	T result=0;
	for (int i = 0; i < a.size(); i++)
		result += a[i] * b[i];
	return result;
}
//倍乘
template<typename T>
std::vector<T> coef(const std::vector<T>& a, const T k)
{
	static_assert(std::is_arithmetic<T>::value == true, "type error");
	if (a.empty())
		throw "vector empty";
	std::vector<T> result(a.size());
	for (int i = 0; i < a.size(); i++)
		result[i] = a[i] * k;
	return result;
}
//计算2vector元素是否成比例，如果成比例，则返回比例，否则返回0；
template<typename T>
T proport(const std::vector<T>& a, const std::vector<T>& b)
{
	static_assert(std::is_arithmetic<T>::value == true, "type error");
	static_assert(a.size() == b.size()&&!a.empty()&&!b.empty(), "vector error");
	if (b[0] == T(0))
		return T(0);
	T coef = a[0] / b[0];
	for (int i = 0; i < a.size(); i++)
	{
		if (a[i] != b[i] * coef)
			return T(0);
	}
	return coef;
}

template<typename T>
T square(const std::vector<T>& a)
{
	static_assert(std::is_arithmetic<T>::value == true, "type error");
	T result{};
	for (auto& val : a)
		result += val * val;
	return result;
}

template<typename T>
T norm(const std::vector<T>& a)
{
	return std::sqrt(square(a));
}


END_FUNCTIONS
END_MATH


										


												