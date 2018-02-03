#pragma once
#include "Eigen\Core"

class Curve
{
public:
	Curve();
	~Curve();

	void reset();
	void addControlPoints(double p0[3], double p1[3]); // add point from the selecting line intersecting at the xy plane
public:
	enum Curve_Type {Bezier, Cubic_B_spline, Quadric_B_spline};
	enum Control_Type {ADD, MOVE, VIEW};

	Curve_Type		m_curveType;
	Control_Type	m_contrlType;

	int n_ctls;                // number of control points
	Eigen::Matrix3Xd m_ctls;   // control points
};

