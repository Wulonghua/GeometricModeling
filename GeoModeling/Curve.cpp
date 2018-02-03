#include "Curve.h"

Curve::Curve():n_ctls(0),m_curveType(Bezier),m_contrlType(ADD)
{
	m_ctls = Eigen::Matrix3d::Zero(3, 1000);      // assume that there will be no more than 1000 control points.
}


Curve::~Curve()
{

}

void Curve::reset()
{
	n_ctls = 0;
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
}

