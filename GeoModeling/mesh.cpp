// ------------------------------------------------------------
// Mesh.cpp:  Implementation of mesh class
// ------------------------------------------------------------

#include "Mesh.h"

// ------------------------------------------------------------
// AddFacet:  Adds a triangle to the mesh.
//            This is one of 2 functions that can be used to build a mesh
// ------------------------------------------------------------
void Mesh::AddFacet(datatype x1, datatype y1, datatype z1, datatype x2, datatype y2, datatype z2, datatype x3, datatype y3, datatype z3) {
	
	vector<GeomVert> geomfacet;
	geomfacet.push_back( GeomVert(x1, y1, z1) );
	geomfacet.push_back( GeomVert(x2, y2, z2) );
	geomfacet.push_back( GeomVert(x3, y3, z3) );

	AddFacet( geomfacet );
}
void Mesh::AddFacet(datatype x1, datatype y1, datatype z1, datatype x2, datatype y2, datatype z2, datatype x3, datatype y3, datatype z3, datatype x4, datatype y4, datatype z4)
{
	vector<GeomVert> geomfacet;
	geomfacet.push_back(GeomVert(x1, y1, z1));
	geomfacet.push_back(GeomVert(x2, y2, z2));
	geomfacet.push_back(GeomVert(x3, y3, z3));
	geomfacet.push_back(GeomVert(x4, y4, z4));

	AddFacet(geomfacet);

	//for rendering gldrawarrays
	renderVerts.push_back(float(x1)); renderVerts.push_back(float(y1)); renderVerts.push_back(float(z1));
	renderVerts.push_back(float(x2)); renderVerts.push_back(float(y2)); renderVerts.push_back(float(z2));
	renderVerts.push_back(float(x3)); renderVerts.push_back(float(y3)); renderVerts.push_back(float(z3));
	renderVerts.push_back(float(x1)); renderVerts.push_back(float(y1)); renderVerts.push_back(float(z1));
	renderVerts.push_back(float(x3)); renderVerts.push_back(float(y3)); renderVerts.push_back(float(z3));
	renderVerts.push_back(float(x4)); renderVerts.push_back(float(y4)); renderVerts.push_back(float(z4));

	double a[3], b[3];
	float c[3];
	a[0] = x2 - x1; a[1] = y2 - y1; a[2] = z2 - z1;
	b[0] = x3 - x1; b[1] = y3 - y1; b[2] = z3 - z1;
	c[0] = float(a[1] * b[2] - a[2] * b[1]);
	c[1] = float(a[2] * b[0] - a[0] * b[2]);
	c[2] = float(a[0] * b[1] - a[1] * b[0]);
	Eigen::Vector3f cv(c[0], c[1], c[2]);
	cv.normalize();
	for (int i = 0; i < 6; ++i)
	{
		renderNormals.push_back(cv[0]);
		renderNormals.push_back(cv[1]);
		renderNormals.push_back(cv[2]);
	}
}
// ------------------------------------------------------------





// ------------------------------------------------------------
// AddFacet:  Adds a facet with arbitrary number of vertices to mesh
//            This is one of 2 functions that can be used to build a mesh
// ------------------------------------------------------------
void Mesh::AddFacet(vector<GeomVert> geomfacet) {
	int i;

	// --------------
	// Create topo facet (list of geom vertex indices)
	TopoFacet topofacet;
	// Look for facet vertices in mesh - if they don't already exist in mesh then add them
	for (i = 0; i < geomfacet.size(); i++) {
		int v_ind = FindGeomVertex( geomfacet[i] );
		if (v_ind == -1) {
			// New vertex:  add geomtric vertex
			v_ind = mGeomVerts.size();
			mGeomVerts.push_back( geomfacet[i] );

			// Add topo vertex
			TopoVert topovert;
			mTopoVerts.push_back( topovert );
		}

		// Add vertex indice to topo facet
		topofacet.AddIncVertex( v_ind );
	}

	// Add this new topo facet to mesh	
	int facet_ind = mTopoFacets.size();
	mTopoFacets.push_back( topofacet );


	// Add edges of facet to mesh, again checking if they already exist
	for (i = 0; i < topofacet.GetNumberVertices(); i++) {
		int prev = (i == 0) ? topofacet.GetNumberVertices() - 1  :  i - 1;
		
		// Create edge
		TopoEdge e;
		e.SetVertex(0, topofacet.GetVertexInd(prev) );
		e.SetVertex(1, topofacet.GetVertexInd(i) );

		// Check if exists
		int e_ind = FindTopoEdge( e );
		
		if (e_ind == -1) {
			// Didn't exist, add to mesh
			e_ind = mTopoEdges.size();
			mTopoVerts[ e.GetVertex(0) ].AddIncEdge( e_ind );
			mTopoVerts[ e.GetVertex(1) ].AddIncEdge( e_ind );
			mTopoEdges.push_back( e );			
		}

		// Point edge to this facet
		mTopoEdges[e_ind].AddIncFacet( facet_ind );

		// Point facet to this edge
		mTopoFacets[ facet_ind ].AddIncEdge( e_ind );
	}
	// --------------
		

	
	// Compute other connectivity
	for (i = 0; i < topofacet.GetNumberVertices(); i++) {
		// Add vertex-facet topology
		mTopoVerts[  topofacet.GetVertexInd(i) ].AddIncFacet( facet_ind );

		// Add vertex-vertex (edge) topology
		int prev = (i == 0) ? topofacet.GetNumberVertices() - 1  :  i - 1;
		int next = (i == topofacet.GetNumberVertices() - 1) ? 0 : i + 1;

		mTopoVerts[ topofacet.GetVertexInd(i) ].AddIncVert( topofacet.GetVertexInd( prev ) );
		mTopoVerts[ topofacet.GetVertexInd(i) ].AddIncVert( topofacet.GetVertexInd( next ) );
	}
	
	// Facet-facet adjacency...
	for (i = 0; i < mTopoFacets[ facet_ind ].GetNumberEdges(); i++) {		
		TopoEdge edge = mTopoEdges[ mTopoFacets[ facet_ind ].GetIncEdge(i) ];
		for (int j = 0; j < edge.GetNumberIncFacets(); j++) {
			if (edge.GetIncFacet(j) != facet_ind) {
				mTopoFacets[ facet_ind ].AddIncFacet( edge.GetIncFacet(j) );
				mTopoFacets[ edge.GetIncFacet(j) ].AddIncFacet( facet_ind );
			}
		}
	}
}
// ------------------------------------------------------------





// ------------------------------------------------------------
// Erase:  Releases all memory used by object
// ------------------------------------------------------------
void Mesh::Erase() {
	mGeomVerts.clear();
	mTopoVerts.clear();
	mTopoEdges.clear();
	mTopoFacets.clear();
	renderVerts.clear();
	renderNormals.clear();
}
// ------------------------------------------------------------

void Mesh::RevolveYaxis(const Eigen::MatrixXd & curve_pts, int n_curve_pts)
{
	Eigen::MatrixXd s0 = curve_pts.leftCols(n_curve_pts);
	Eigen::MatrixXd s1 = s0;
	Eigen::MatrixXd s2 = s0;
	double ang0 = PI2 / n_slice;
	for (int i = 1; i < n_slice; ++i) 
	{
		double ang = ang0*i;
		double cos_ang = cos(ang);
		double sin_ang = sin(ang);
		for (int j = 0; j < n_curve_pts; ++j) 
		{
			double r = s0(0, j);
			s2(0, j) = r*cos_ang;
			s2(2, j) = r*sin_ang;
		}

		//add facet
		for (int k = 0; k < n_curve_pts-1; ++k) 
		{
			AddFacet(s1(0, k), s1(1, k), s1(2, k),
					 s2(0, k), s2(1, k), s2(2, k),
					 s2(0, k + 1), s2(1, k + 1), s2(2, k + 1),
					 s1(0, k + 1), s1(1, k + 1), s1(2, k + 1));
		}
		s1 = s2;
	}
	s2 = s0;
	for (int k = 0; k < n_curve_pts - 1; ++k)
	{
		AddFacet(s1(0, k), s1(1, k), s1(2, k),
			s2(0, k), s2(1, k), s2(2, k),
			s2(0, k + 1), s2(1, k + 1), s2(2, k + 1),
			s1(0, k + 1), s1(1, k + 1), s1(2, k + 1));
	}
}

// ------------------------------------------------------------
// FindGeomVertex:  Searches for a geometric vertex in the mesh,
//                  returning its indice if found, -1 otherwise
// ------------------------------------------------------------
int Mesh::FindGeomVertex(GeomVert v) {
	for (int i = 0; i < mGeomVerts.size(); i++) {
		if (mGeomVerts[i] == v) return i;
	}
	return -1;
}
// ------------------------------------------------------------

// ------------------------------------------------------------
// FindTopoEdge:  Searches for an edge in the mesh, returing 
//                its indice if found, -1 otherwise
// ------------------------------------------------------------
int	Mesh::FindTopoEdge(TopoEdge e) {
	for (int i = 0; i < mTopoEdges.size(); i++) {
		if (mTopoEdges[i] == e) return i;
	}
	return -1;
}
// ------------------------------------------------------------

