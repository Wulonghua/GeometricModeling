#include "Curve.h"

Curve::Curve():n_ctls(0),n_points(0),m_curveType(Bezier),m_contrlType(ADD)
{
	m_ctls = Eigen::Matrix3Xd::Zero(3, 1000);      // assume that there will be no more than 1000 control points.
	m_points = Eigen::Matrix3Xd::Zero(3, 10000);   // assume that there will be no more than 10000 rendering points.
}


Curve::~Curve()
{

}

void Curve::reset()
{
	n_ctls = 0;
	n_points = 0;
}

void Curve::addControlPoints(double p0[3], double p1[3])
{
	/************line equation*****************
	x = x0 + (x1-x0)*t
	y = y0 + (y1-y0)*t
	z = z0 + (z1-z0)*t
	******************************************/
	double t = p0[2] / (p0[2] - p1[2]);
	double x = p0[0] + (p1[0] - p0[0])*t;
	double y = p0[1] + (p1[1] - p0[1])*t;
	m_ctls(0, n_ctls) = x;
	m_ctls(1, n_ctls) = y;
	n_ctls++;

	if (m_curveType == Bezier) 
	{
		n_points = 100;
		double itv = 1.0 / (n_points - 1);
		int n = n_ctls - 1;
		computeBersteins(n);
		for (int k = 0; k < n_points + 1; ++k)
		{
			double u = itv*double(k);
			Eigen::Vector3d p = Eigen::Vector3d::Zero();
			for (int i = 0; i < n_ctls; ++i)
			{
				p += m_ctls.col(i)*m_bernPoly[i] * std::pow(u, i)*std::pow(1 - u, n - i);
			}
			m_points.col(k) = p;
		}
	}
}

int Curve::getBinomialCoeff(int n, int i)
{
	if(i==0 || i==n)
		return 1;
	int s = 1;
	for (int k = 1; k <= i; ++k)
		s *= (n + 1 - k) / k;
	return s;
}

void Curve::computeBersteins(int n)
{
	m_bernPoly.resize(n + 1, 0.0);
	for (int i = 0; i < n + 1; ++i)
		m_bernPoly[i] = getBinomialCoeff(n, i);
}

