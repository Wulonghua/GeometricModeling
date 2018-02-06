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
	void addControlPoints(double p0[3], double p1[3]); // add point from the selecting line intersecting at the xy plane
	void generateBezierPoints();
	void generateQuadBspline();
	void generateCubicBspline();

private:
	double getBinomialCoeff(int n, int i);    // compute C_n^i
	void computeBersteins(int n);


public:
	enum Curve_Type {Bezier, Cubic_B_spline, Quadric_B_spline};
	enum Control_Type {ADD, MOVE, VIEW};

	Curve_Type		m_curveType;
	Control_Type	m_contrlType;

	int n_ctls;                // number of control points
	int n_points;              // number of rendering points
	int n_precision;           // number of segments in one interval

	Eigen::Matrix3Xd m_ctls;   // control points
	Eigen::Matrix3Xd m_points; // rendering points
	std::vector<double> m_bernPoly;

};

