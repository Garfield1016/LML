#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <memory>

#include "math/functions.hpp"
#include "math/Vector.hpp"
#include "cuutils/cputimer.h"
#include "math/matrix.hpp"
#include "math/numerical.hpp"
using namespace math;
using namespace functions;
using namespace vec;
using namespace numerical;
void t_eigen()
{

	Eigen::Matrix<double, 3, 3> m1;
	m1 << 1, 2, 3,
		4, 5, 6,
		7, 8, 10;
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd(m1);
	//Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd> > lu(m1);
	//std::cout << m1 << std::endl;
	
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
void t_math()
{
	TimerClock timer;
	timer.reset();
	std::vector<int> a(10, 1);
	std::vector<int> b(10, 1);
	std::vector<int> res = a / b;
	timer.update();
	//timer.print();
	std::vector<int> ssss(5, 2);

	//Vector<int> vec(5);
	//Vector<int> vec2(5, 2);
	double arr[] = { 1,1,1,1,1 };
	Vector<double> vec3(arr, 5);
	Vector<double> vec4 = vec3;


	Vector<double> vec5 = vec3.normalize();
	Vector<double> vec6({ 1,1,1,1,1 });
	double aaaaaaaa = std::numeric_limits<double>::epsilon();
	std::vector<int> vec7(2, 12);
	std::vector<int> vec8;
	vec8 = vec7;
	vec8[0] = 9;

	//Data data;
	std::shared_ptr<int[]> data;
	data = std::make_shared<int[]>(2);
	data[0] = 1;
	data[1] = 2;
	int aw = data[1];
	//data = std::make_shared<int>(2);

	std::vector<int> aa = { 1,2,3 };
	std::vector<int> asss=(std::move(aa));
	Matrix<double> m4(3, 3,{1, 2, 3, 4, 5, 6, 7, 8, 10});
	m4(0, 1) = 77;
	double vt = m4(0, 1);
	
	Matrix<double> m4qr(4, 4, {1,2,3,4,5,6,7,8,9,10,11,12,8,9,11,22});
	static Matrix<double> ran = Matrix<double>::randomMN(512, 512); 
	//ran.print();
	TimerClock t_clock;
	t_clock.reset();
	std::vector<Matrix<double>> lup = ran.LUP();
	t_clock.update();
	float tcost = t_clock.getTimerMilliSec();

	Matrix<double> r = lup[0] * lup[1];
	Matrix<double> r2 = lup[2] * ran;
	r2.print();
	r.print();
	/*Matrix<double> m44(3, 3, { 1, 2, 3, 4, 5, 6, 7, 8, 99.7 });
	Matrix<double> km(3, 3, { 2,-1,0,-5,2,-1,0,-1,2 });
	Matrix<double> me(4, 4, { 1,2,3,4,2,1,2,3,3,2,1,2,4,3,2,1 });
	Matrix<double>me2(3, 3, { 1,2,3,3,2,5,1,8,10 });
	Matrix<double>me21(3, 3, { 1,2,3,3,2,5,1,8,10 });
	Matrix<double>me3(2, 2, {4,2,1,5 });
	std::vector<Matrix<double>> qr = me.QRbyHouseholder();
	Matrix tqr = qr[0] * qr[1];


	Matrix<double> inv = me.inverse();
	double det=m4.determinant();
	Matrix<double> ms(3, 2, { 0,1,1,1,1,0 });
	std::vector<Matrix<double>> uevt=ms.SVD();
	Matrix<double>t = uevt[0] * uevt[1] * uevt[2];
	
	Matrix<double> tc(3, 3, { 4,12,-16,12,37,-43,-16,-43,98 });
	std::vector<Matrix<double>>t_c=tc.cholesky();
	Matrix<double> tchol = t_c[0] * t_c[1];
	Matrix<double> met = me2 - me21;
	me2.QRbyGivens();
	Vector<double> y({ 1,8,2 });
	Vector<double> x = sovleLU(m4,y);*/
	int aaad = 1;
}

void t_numerical()
{
	int n = 4;
	std::vector<float> a = { -1, 5, -3, 3, 2 };
	float x = 0.5;
	std::vector<float> b = { 0,0,0,0,0 };
	float val = polyVal<float>(n,a,x,b);
	float x2 = solve(0.0005, 0., 1.);
}
