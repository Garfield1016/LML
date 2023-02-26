#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <memory>

#include "math/functions.hpp"
#include "math/Vector.hpp"
#include "cuutils/cputimer.h"
#include "math/matrix.hpp"
using namespace math;
using namespace functions;
using namespace vec;
void t_eigen()
{

	Eigen::Matrix<double, 3, 3> m1;
	m1 << 1, 2, 3,
		4, 5, 6,
		7, 8, 10;
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd(m1);
	//Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd> > lu(m1);
	std::cout << m1 << std::endl;
	
	//int rrrr = svd.rank();
	int aaa = 1;
}
class base
{
private:
	std::vector<int> data;
public:
	base() = default;
	base(int size) { data.resize(size); std::cout << "construct" << std::endl; }
	base(const base& other) { data = other.data; }
	base(base&& other) { data=(std::move(other.data)); }
	base operator=(base&& b1)
	{
		data = (std::move(b1.data));
		return *this;
	}
};
void t_base(const base& b)
{
	int a = 1;
}


