#pragma once
#include "Eigen\Core"
#include <vector>
#include <cmath>
#include <iostream>

class Curve
{
public:
	Curve();
	~Curve();

	void reset();
	void addControlPoints(double p0[3], double p1[3] , int plane=0); // add point from the selecting line intersecting at the xy plane (plane=0), yz plane(plane =1)
	void setControlPoint(double p0[3], double p1[3], int id, int plane = 0);
	int  pickControlPoint(double p0[3], double p1[3]);
	void generateBezierPoints();
	void generateQuadBspline();
	void generateCubicBspline();
	void generateCurves();
	std::vector<Eigen::Vector3d> Subdivide(std::vector<Eigen::Vector3d> points, int m, double u);

	void changeCloseStatus(bool closed);
	bool isValidYrevolve();

private:
	double getBinomialCoeff(int n, int i);    // compute C_n^i
	void computeBersteins(int n);
	std::vector<Eigen::Vector3d> cancatenatePoints(std::vector<Eigen::Vector3d> poly1, std::vector<Eigen::Vector3d> poly2);
	std::vector<Eigen::Vector3d> OneSubdivide(std::vector<Eigen::Vector3d> points,
		std::vector<Eigen::Vector3d>& poly1, std::vector<Eigen::Vector3d>& poly2, double u);
	std::vector<Eigen::Vector3d> getControlPointsVector();
	void setRenderPoints(const std::vector<Eigen::Vector3d>& v);

public:
	enum Curve_Type {Bezier, Cubic_B_spline, Quadric_B_spline};
	enum Control_Type {ADD, MOVE, VIEW};
	enum Gen_Type { SAMPLE, SUBDIVISION };

	Curve_Type		m_curveType;
	Control_Type	m_contrlType;
	Gen_Type		m_genType;

	int n_ctls;                // number of control points
	int n_points;              // number of rendering points
	int n_precision;           // number of segments in one interval

	Eigen::Matrix3Xd m_ctls;   // control points
	Eigen::Matrix3Xd m_points; // rendering points
	std::vector<double> m_bernPoly;

	bool m_closed;

};

