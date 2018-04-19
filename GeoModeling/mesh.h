
#ifndef MESHCLASS
#define	MESHCLASS

#include "Eigen\Core"
#include "Eigen\geometry"
#include <stdio.h>
#include <vector>
#include <set>
#include <cmath>
#include <memory>
#include <algorithm>
#include <iostream>
#include <QFile>
#include <QTextStream>
#include <QString>

#define PI2	6.2831853071796
// Change to double for more precision
typedef	double	datatype;
using namespace std;

// ------------------------------------------------------------
// GeomVert: this class holds the geometric coordinates of a vertex
// ------------------------------------------------------------
class GeomVert {
public:
	GeomVert() { mCo[0] = 0.0; mCo[1] = 0.0; mCo[2] = 0.0; }
	GeomVert(datatype x, datatype y, datatype z) { mCo[0] = x; mCo[1] = y; mCo[2] = z; }
	GeomVert(datatype x, datatype y) { mCo[0] = x; mCo[1] = y; mCo[2] = 0.0; }

	datatype      GetCo(int axis)                { return mCo[axis]; }
	datatype	  Norm() { return std::sqrt(mCo[0]*mCo[0]+mCo[1]*mCo[1]+mCo[2]*mCo[2]); }

	bool operator == (GeomVert &A)               { 
		return ( (mCo[0] == A.GetCo(0)) && (mCo[1] == A.GetCo(1)) && (mCo[2] == A.GetCo(2)) );		
	}
 
	GeomVert operator + (GeomVert &A) {
		return GeomVert(A.GetCo(0) + mCo[0], A.GetCo(1) + mCo[1], A.GetCo(2) + mCo[2]);
	}

	GeomVert operator - (GeomVert &A) {
		return GeomVert(-A.GetCo(0) + mCo[0], -A.GetCo(1) + mCo[1], -A.GetCo(2) + mCo[2]);
	}

	GeomVert operator / (int n) {
		datatype nf = datatype(n);
		return GeomVert(mCo[0] / nf, mCo[1] / nf, mCo[2] / nf);
	}

	GeomVert operator * (datatype s) {
		return GeomVert(mCo[0] * s, mCo[1] * s, mCo[2] * s);
	}

	GeomVert operator * (int s) {
		return GeomVert(mCo[0] * s, mCo[1] * s, mCo[2] * s);
	}
	
	GeomVert& operator = (GeomVert& A) {
		for (int i = 0; i < 3; ++i)
			mCo[i] = A.GetCo(i);
		return *this;
	}

	GeomVert& operator += (GeomVert& A) {
		for (int i = 0; i < 3; ++i)
			mCo[i] += A.GetCo(i);
		return *this;
	}

private:
	datatype	mCo[3];
};
// ------------------------------------------------------------

// ------------------------------------------------------------
// TopoVert: this class holds all of a vertex's topological (connectivity) information
// ------------------------------------------------------------
class TopoVert {
public:
	TopoVert()                      { };
	~TopoVert()                     { mIncVerts.clear(); mIncEdges.clear(); mIncFacets.clear(); }
	void AddIncVert (int vert_ind ) { mIncVerts.insert( vert_ind ); }
	void AddIncEdge (int edge_ind ) { mIncEdges.push_back( edge_ind ); }
	void AddIncFacet(int facet_ind) { mIncFacets.push_back( facet_ind ); }
	
	int GetNumberIncVertices()		{ return mIncVerts.size(); }
	int GetIncVertex(int vert_ind)  { 
		set<int>::iterator sit = mIncVerts.begin(); for (int i = 0; i < vert_ind; i++) sit++;		
		return *sit;
	}
	int GetNumberIncEdges()			{ return mIncEdges.size(); }
	int GetIncEdge(int edge_ind)    { return mIncEdges[edge_ind]; }
	int GetNumberIncFacets()	    { return mIncFacets.size(); }
	int GetIncFacet(int facet_ind)  { return mIncFacets[facet_ind]; }

private:
	set<int>    mIncVerts;
	vector<int> mIncEdges;
	vector<int> mIncFacets;
	
};
// ------------------------------------------------------------





// ------------------------------------------------------------
// TopoEdge
// ------------------------------------------------------------
class TopoEdge {
public:
	TopoEdge()                      { v1 = v2 = -1; }
	~TopoEdge()                     { mIncFacets.clear(); }

	bool operator == (TopoEdge &A)  {
		return (  ((v1 == A.GetVertex(0)) && (v2 == A.GetVertex(1))) || ((v2 == A.GetVertex(0)) && (v1 == A.GetVertex(1))) );
	}

	int  GetVertex(int ind)         { if (ind == 0) return v1;  return v2; }
	void SetVertex(int ind, int v)  { if (ind == 0) { v1 = v; } else { v2 = v; } }
	
	void AddIncFacet(int facet_ind) { mIncFacets.push_back(facet_ind); }
	int  GetNumberIncFacets()       { return mIncFacets.size(); }
	int  GetIncFacet(int facet_ind) { return mIncFacets[facet_ind]; }
	int  GetOtherIncFacet(int fi) {
		if (mIncFacets.size() < 2) return -1;
		if (mIncFacets[0] == fi) {
			return mIncFacets[1];
		}
		else if (mIncFacets[1] == fi) {
			return mIncFacets[0];
		}
		else {
			return -1;
		}
	}

private:
	int v1, v2;
	vector<int> mIncFacets;
};
// ------------------------------------------------------------





// ------------------------------------------------------------
// TopoFacet:  this class holds a facet's topological connectivity) information
//             Facets are represented as a list of vertex indices
// ------------------------------------------------------------
class TopoFacet {
public:
	TopoFacet()                     { };
	~TopoFacet()                    { mIncVerts.clear(); mIncEdges.clear();  mIncFacets.clear(); }
	void AddIncVertex(int v_ind)    { mIncVerts.push_back( v_ind ); }
	void AddIncEdge(int e_ind)      { mIncEdges.push_back( e_ind ); }
	void AddIncFacet(int f_ind)     { mIncFacets.insert( f_ind ); }
	int  GetNumberVertices()        { return mIncVerts.size(); }
	int  GetVertexInd(int vert_ind) { return mIncVerts[vert_ind]; }
	int  GetVertexOrder(int id) {
		int order = std::find(mIncVerts.begin(), mIncVerts.end(), id)-mIncVerts.begin();
		return order >= mIncVerts.size() ? -1 : order;
	}

	int  GetNumberEdges()		    { return mIncEdges.size(); }
	int  GetIncEdge(int edge_ind)   { return mIncEdges[edge_ind]; }
	int  GetNumberFacets()		    { return mIncFacets.size(); }
	int  GetIncFacet(int facet_ind) { 
		set<int>::iterator sit = mIncFacets.begin(); for (int i = 0; i < facet_ind; i++) sit++;		
		return *sit;
	}

private:
	vector<int> mIncVerts;	
	vector<int> mIncEdges;
	set<int>    mIncFacets;
};
// ------------------------------------------------------------

// ------------------------------------------------------------
// Mesh:  This class uses all the preceding classes to represent a mesh with
//        adjacency.connectivity information
// ------------------------------------------------------------
class Mesh {
public:
	Mesh():n_slice(3), m_depth(0.5),m_buildType(DONOTHING) { };
	~Mesh() { Erase(); };

	void	AddFacet(datatype x1, datatype y1, datatype z1, datatype x2, datatype y2, datatype z2, datatype x3, datatype y3, datatype z3);
	void	AddFacet(datatype x1, datatype y1, datatype z1, datatype x2, datatype y2, datatype z2,
		datatype x3, datatype y3, datatype z3, datatype x4, datatype y4, datatype z4);
	void	AddFacet(vector<GeomVert> geomfacet);
	void	AddFacet(TopoFacet topofacet); // use for adding facet from *.off file, alreay got vertices posision
	void	AddNewVertex(GeomVert v) { mGeomVerts.push_back(v); mTopoVerts.push_back(TopoVert()); }
	int		GetNumberVertices()		  { return mGeomVerts.size(); }
	int		GetNumberEdges()			  { return mTopoEdges.size(); }
	int		GetNumberFacets()           { return mTopoFacets.size(); }

	TopoVert  GetVertex(int vert_ind)     { return mTopoVerts[vert_ind]; }
	TopoEdge  GetEdge(int edge_ind)       { return mTopoEdges[edge_ind]; }
	TopoFacet GetFacet(int facet_ind)     { return mTopoFacets[facet_ind]; }
	
	GeomVert  GetGeomVertex(int vert_ind) { return mGeomVerts[vert_ind]; }

	void RevolveYaxis(const Eigen::MatrixXd& curve_pts,int n_curve_pts);
	void ExtrusionZaxis(const Eigen::MatrixXd& curve_pts, int n_curve_pts);
	void Sweep(const Eigen::MatrixXd& pts, int n_pts, const Eigen::MatrixXd& traj, int n_traj, bool closed);



private:
	int		FindGeomVertex(GeomVert v);
	int		FindTopoEdge(TopoEdge e);
	void	Erase();
	void	computeFaceEdgeCenters();
	bool    isAdjFacesRightOrder(const std::vector<int>& fs, int vi);

	vector<GeomVert>  mGeomVerts;
	vector<TopoVert>  mTopoVerts;
	vector<TopoEdge>  mTopoEdges;
	vector<TopoFacet> mTopoFacets;

	vector<vector<GeomVert>> mFEcenters;  // Edge centers in each face
	vector<GeomVert>		 mEcenters;
	vector<GeomVert>		 mFcenters;

	vector<Eigen::Vector3f>	 mFaceNormals;
	vector<Eigen::Vector3f>	 mVertNormals;


public:
	enum Build_Type { DONOTHING, REVOLUTION, EXTRUSION, SWEEP };
	void reset() { Erase(); m_buildType = DONOTHING; }
	void prepareRender();
	void saveMesh(QString filename);
	void LoadModel(QString filepath);
	void SubDooSabin(std::shared_ptr<Mesh> mesh);
	void SubCatmullClark(std::shared_ptr<Mesh> mesh);
	void SubLoop(std::shared_ptr<Mesh> mesh);

	bool NNCrust(std::vector<int>& graph);

	// for rendering purpose
	vector<float> renderVerts;
	vector<float> renderNormals;
	int		n_slice;
	double	m_depth;
	Build_Type m_buildType;



};
// ------------------------------------------------------------

#endif