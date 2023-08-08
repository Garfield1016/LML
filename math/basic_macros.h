#pragma once
#include <concepts>
#include <iostream>
#include <string>
#define BEGIN_MATH		\
namespace math{			\

#define END_MATH }

#define BEGIN_FUNCTIONS	\
namespace functions{		\

#define END_FUNCTIONS }

#define BEGIN_VEC	\
namespace vec{		\

#define END_VEC }

#define BEGIN_NUMERICAL namespace numerical{ 
#define END_NUMERICAL }

BEGIN_MATH
constexpr double epsilon = 1e-15;
template<typename T>
concept arithmetic = std::is_arithmetic_v<T>;
template<typename T>
concept floatingtype = std::is_floating_point_v<T>;
END_MATH


