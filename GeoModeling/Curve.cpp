#include "Curve.h"

Curve::Curve():
	n_ctls(0),n_points(0),n_precision(5),
	m_curveType(Bezier),
	m_contrlType(ADD),
	m_genType(SAMPLE)
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
		generateBezierPoints();
	}
	else if (m_curveType == Quadric_B_spline)
	{
		generateQuadBspline();
	}
	else if (m_curveType == Cubic_B_spline)
	{
		generateCubicBspline();
	}
}

double Curve::getBinomialCoeff(int n, int i)
{
	if(i==0 || i==n)
		return 1;
	double s = 1;
	for (int k = 1; k <= i; ++k)
		s *= double(n + 1 - k) / double(k);
	return s;
}

void Curve::computeBersteins(int n)
{
	m_bernPoly.resize(n + 1, 0.0);
	for (int i = 0; i < n + 1; ++i)
		m_bernPoly[i] = getBinomialCoeff(n, i);
}

std::vector<Eigen::Vector3d> Curve::cancatenatePoints(std::vector<Eigen::Vector3d> poly1, std::vector<Eigen::Vector3d> poly2)
{	
	if (poly1.size() == 0 && poly2.size() != 0)
		return poly2;
	if (poly2.size() == 0 && poly1.size() != 0)
		return poly1;
	std::vector<Eigen::Vector3d> poly;
	if (poly1.back() == poly2.front())
	{
		poly.insert(poly.end(),poly1.begin(),poly1.end()-1);
		poly.insert(poly.end(),poly2.begin(),poly2.end());
	}
	else
	{
		poly.insert(poly.end(), poly1.begin(), poly1.end());
		poly.insert(poly.end(), poly2.begin(), poly2.end());
	}
	return poly;
}

std::vector<Eigen::Vector3d> Curve::OneSubdivide(std::vector<Eigen::Vector3d> points, std::vector<Eigen::Vector3d>& poly1, std::vector<Eigen::Vector3d>& poly2, double u)
{
	std::vector<Eigen::Vector3d> poly;
	int n = points.size() - 1;
	if (n == 0) 
	{
		poly = poly1;
		poly.push_back(points[0]);
		return cancatenatePoints(poly, poly2);
	}
	else
	{
		Eigen::Vector3d p0 = points.front();
		Eigen::Vector3d pn = points.back();
		poly1.push_back(p0);
		poly2.insert(poly2.begin(), pn);
		std::vector<Eigen::Vector3d> pts;
		pts.resize(n);
		for (int i = 0; i < n; ++i) 
		{
			pts[i] = points[i] + u * (points[i + 1] - points[i]);
		}
		return OneSubdivide(pts, poly1, poly2, u);
	}
}

std::vector<Eigen::Vector3d> Curve::getControlPointsVector()
{
	std::vector<Eigen::Vector3d> v;
	v.resize(n_ctls);
	for (int i = 0; i < n_ctls; ++i)
	{
		v[i] = m_ctls.col(i);
	}
	return v;
}

void Curve::setRenderPoints(const std::vector<Eigen::Vector3d>& v)
{
	n_points = v.size();
	for (int i = 0; i < n_points; ++i)
	{
		m_points.col(i) = v[i];
	}
}

void Curve::generateBezierPoints()
{
	if (n_ctls < 2) return;
	if (m_genType == SAMPLE)
	{
		n_points = (n_ctls - 1)*n_precision + 1;
		double itv = 1.0 / (n_points - 1);
		int n = n_ctls - 1;
		computeBersteins(n);
		for (int k = 0; k < n_points; ++k)
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
	else if (m_genType == SUBDIVISION)
	{
		std::vector<Eigen::Vector3d> pts = getControlPointsVector();
		std::vector<Eigen::Vector3d> points = Subdivide(pts, n_precision, 0.5);
		setRenderPoints(points);
	}
}

void Curve::generateQuadBspline()
{
	if (n_ctls < 3) return;
	int n = n_ctls - 2;
	n_points = n*n_precision;
	double itv = 1.0 / n_precision;
	double u;
	Eigen::Matrix3d M, P;
	Eigen::Vector3d U;
	M << 0.5, -1, 0.5,
		-1, 1, 0.5,
		0.5, 0, 0;
	// use column order matrix/vector, so the matrix form is transposed to 
	// that of the form in the course note.
	for (int i = 0; i < n; ++i)
	{
		P = m_ctls.block<3, 3>(0,i);
		for (int j = 0; j < n_precision; ++j)
		{
			int k = i*n_precision + j;
			u = /*double(i) + */itv*j;
			U[0] = u*u; U[1] = u; U[2] = 1;
			m_points.col(k) = P*M*U;
		}
	}
}

void Curve::generateCubicBspline()
{
	if (n_ctls < 4) return;
	int n = n_ctls - 3;
	n_points = n*n_precision;
	double itv = 1.0 / n_precision;
	double u;
	Eigen::Matrix4d M;
	Eigen::Matrix3Xd P;
	Eigen::Vector4d U;
	M << -1, 3, -3, 1,
		3, -6, 0, 4,
		-3, 3, 3, 1,
		1, 0, 0, 0;
	M = M / 6.0;
	// use column order matrix/vector, so the matrix form is transposed to 
	// that of the form in the course note.
	for (int i = 0; i < n; ++i)
	{
		P = m_ctls.block<3, 4>(0, i);
		for (int j = 0; j < n_precision; ++j)
		{
			int k = i*n_precision + j;
			u = /*double(i) + */itv*j;
			U[0] = u*u*u; U[1] = u*u; U[2] = u; U[3] = 1;
			m_points.col(k) = P*M*U;
		}
	}
}

void Curve::generateCurves()
{
	if (m_curveType = Bezier)
		generateBezierPoints();
	else if (m_curveType = Quadric_B_spline)
		generateQuadBspline();
	else if (m_curveType = Cubic_B_spline)
		generateCubicBspline();
}

std::vector<Eigen::Vector3d> Curve::Subdivide(std::vector<Eigen::Vector3d> points, int m, double u)
{
	std::vector<Eigen::Vector3d> poly1;
	std::vector<Eigen::Vector3d> poly2;
	if (m == 1)
	{
		return	OneSubdivide(points, poly1, poly2, u);
	}
	else
	{
		std::vector<Eigen::Vector3d> poly_half1;
		std::vector<Eigen::Vector3d> poly_half2;
		std::vector<Eigen::Vector3d> pts;
		pts = OneSubdivide(points, poly1, poly2, u);
		int n = (pts.size() - 1) / 2;
		std::vector<Eigen::Vector3d> pts1(pts.begin(),pts.begin()+n+1);
		std::vector<Eigen::Vector3d> pts2(pts.begin()+n,pts.end());
		poly_half1 = Subdivide(pts1, m - 1, u);
		poly_half2 = Subdivide(pts2, m - 1, u);
		return cancatenatePoints(poly_half1, poly_half2);
	}
}

