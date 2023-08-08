#include "basic_macros.h"

BEGIN_MATH
BEGIN_NUMERICAL

/*
* @param n 阶数
* @param a 系数 低->高
* @param b 插值项 默认为0
*/
template <floatingtype T>
T polyVal(int n,std::vector<T>& a,T x,std::vector<T>& b)
{
	T y = a[n];
	for (int i = n - 1; i >= 0; i--)
	{
		y = y * (x - b[i]) + a[i];
		int aa = 1;
	}
		
	return y;
}

template <floatingtype T>
T solve(T tol,T a,T b)
{
	std::vector<T> a_arr = { -1,1,0,1 };
	std::vector<T> b_arr = { 0,0,0,0 };
	int n = 3;
	
	T a_temp = a;
	T b_temp=b;
	T x;
	while (true)
	{
		x = (a_temp + b_temp) * 0.5;
		T y = polyVal(n, a_arr, x, b_arr);
		if ((b_temp - a_temp) < tol)
			break;
		if (y == 0)
			break;
		if (polyVal(n, a_arr, a_temp, b_arr)*y < 0)
			b_temp = x;
		else
			a_temp = x;
	}
	x = (a_temp + b_temp) / 2;
	return x;
}


END_NUMERICAL
END_MATH