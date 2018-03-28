// ------------------------------------------------------------
// Mesh.cpp:  Implementation of mesh class
// ------------------------------------------------------------

#include "Mesh.h"

void Mesh::LoadModel(QString filepath)
{
	QFile file(filepath);
	std::cout << filepath.toStdString() << std::endl;
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cerr << "Cannot load model file." << std::endl;
		exit(1);
	}
	QTextStream fin(&file);
	QString line = fin.readLine().simplified();
	if (line != QString("OFF"))
	{
		std::cerr << "File does not start with OFF" << std::endl;
		exit(1);
	}
	int n_verts, n_faces, n_tmp;
	fin >> n_verts >> n_faces >> n_tmp;
	std::cout << "#verts:" << n_verts << " #faces:" << n_faces << std::endl;
	// add vertices.
	for (int i = 0; i < n_verts; ++i)
	{
		datatype x, y, z;
		fin >> x >> y >> z;
		mGeomVerts.push_back(GeomVert(x, y, z));
		TopoVert topovert;
		mTopoVerts.push_back(topovert);
	}
	// add facets
	for (int i = 0; i < n_faces; ++i)
	{
		TopoFacet topofacet;
		int nv,v_id;
		fin >> nv;
		for (int j = 0; j < nv; ++j)
		{
			fin >> v_id;
			topofacet.AddIncVertex(v_id);
		}
		AddFacet(topofacet);
	}
	prepareRender();
}

void Mesh::SubDooSabin(std::shared_ptr<Mesh> mesh)
{
	computeFaceEdgeCenters();
	vector<vector<GeomVert>> v_new;
	vector<vector<int>>		 fvi_new;
	int id = 0;
	int n_f = mTopoFacets.size();
	v_new.resize(n_f);
	fvi_new.resize(n_f);
	// add F-face
	for (int i = 0; i < n_f; ++i)
	{
		int n_fv = mTopoFacets[i].GetNumberVertices();
		TopoFacet topofacet;
		for (int j = 0; j < n_fv; ++j)
		{
			int k = (j == (n_fv - 1) ? 0 : j + 1);
			GeomVert v = (mGeomVerts[mTopoFacets[i].GetVertexInd(j)] +
				mFEcenters[i][j] +
				mFEcenters[i][k] +
				mFcenters[i])*0.25;
			v_new[i].push_back(v);
			topofacet.AddIncVertex(id);
			fvi_new[i].push_back(id++);
			mesh->AddNewVertex(v);
		}
		mesh->AddFacet(topofacet);
	}

	// add E-face
	int n_e = mTopoEdges.size();
	for (int i = 0; i < n_e; ++i)
	{
		TopoFacet topofacet;
		if (mTopoEdges[i].GetNumberIncFacets() < 2) continue;
		int v0 = mTopoEdges[i].GetVertex(0);
		int v1 = mTopoEdges[i].GetVertex(1);
		int f0 = mTopoEdges[i].GetIncFacet(0);
		int f1 = mTopoEdges[i].GetIncFacet(1);
		int a, b, c, d;
		a = mTopoFacets[f0].GetVertexOrder(v0);
		b = mTopoFacets[f0].GetVertexOrder(v1);
		c = mTopoFacets[f1].GetVertexOrder(v0);
		d = mTopoFacets[f1].GetVertexOrder(v1);
		if (b - a == 1 || a - b == mTopoFacets[f0].GetNumberVertices() - 1)
		{
			topofacet.AddIncVertex(fvi_new[f0][b]);
			topofacet.AddIncVertex(fvi_new[f0][a]);
			topofacet.AddIncVertex(fvi_new[f1][c]);
			topofacet.AddIncVertex(fvi_new[f1][d]);
		}
		else
		{
			topofacet.AddIncVertex(fvi_new[f0][a]);
			topofacet.AddIncVertex(fvi_new[f0][b]);
			topofacet.AddIncVertex(fvi_new[f1][d]);
			topofacet.AddIncVertex(fvi_new[f1][c]);
		}
	}

	//Todo: add V-face

}

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

void Mesh::AddFacet(TopoFacet topofacet)
{
	int i;

	// Add this new topo facet to mesh	
	int facet_ind = mTopoFacets.size();
	mTopoFacets.push_back(topofacet);

	// Add edges of facet to mesh, again checking if they already exist
	for (i = 0; i < topofacet.GetNumberVertices(); i++) {
		int prev = (i == 0) ? topofacet.GetNumberVertices() - 1 : i - 1;

		// Create edge
		TopoEdge e;
		e.SetVertex(0, topofacet.GetVertexInd(prev));
		e.SetVertex(1, topofacet.GetVertexInd(i));

		// Check if exists
		int e_ind = FindTopoEdge(e);

		if (e_ind == -1) {
			// Didn't exist, add to mesh
			e_ind = mTopoEdges.size();
			mTopoVerts[e.GetVertex(0)].AddIncEdge(e_ind);
			mTopoVerts[e.GetVertex(1)].AddIncEdge(e_ind);
			mTopoEdges.push_back(e);
		}

		// Point edge to this facet
		mTopoEdges[e_ind].AddIncFacet(facet_ind);

		// Point facet to this edge
		mTopoFacets[facet_ind].AddIncEdge(e_ind);
	}
	// --------------

	// Compute other connectivity
	for (i = 0; i < topofacet.GetNumberVertices(); i++) {
		// Add vertex-facet topology
		mTopoVerts[topofacet.GetVertexInd(i)].AddIncFacet(facet_ind);

		// Add vertex-vertex (edge) topology
		int prev = (i == 0) ? topofacet.GetNumberVertices() - 1 : i - 1;
		int next = (i == topofacet.GetNumberVertices() - 1) ? 0 : i + 1;

		mTopoVerts[topofacet.GetVertexInd(i)].AddIncVert(topofacet.GetVertexInd(prev));
		mTopoVerts[topofacet.GetVertexInd(i)].AddIncVert(topofacet.GetVertexInd(next));
	}

	// Facet-facet adjacency...
	for (i = 0; i < mTopoFacets[facet_ind].GetNumberEdges(); i++) {
		TopoEdge edge = mTopoEdges[mTopoFacets[facet_ind].GetIncEdge(i)];
		for (int j = 0; j < edge.GetNumberIncFacets(); j++) {
			if (edge.GetIncFacet(j) != facet_ind) {
				mTopoFacets[facet_ind].AddIncFacet(edge.GetIncFacet(j));
				mTopoFacets[edge.GetIncFacet(j)].AddIncFacet(facet_ind);
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
void Mesh::prepareRender()
{
	renderVerts.clear();
	renderNormals.clear();

	for (int fi = 0; fi < mTopoFacets.size(); ++fi)
	{
		int n_v = mTopoFacets[fi].GetNumberVertices();
		int vi[3];
		vi[0] = mTopoFacets[fi].GetVertexInd(0);
		// add positions
		for (int i = 1; i < n_v - 1; ++i)
		{
			vi[1] = mTopoFacets[fi].GetVertexInd(i);
			vi[2] = mTopoFacets[fi].GetVertexInd(i+1);
			for (int j = 0; j < 3; ++j)
			for (int k = 0; k < 3; ++k)
			{
				renderVerts.push_back(float(mGeomVerts[vi[j]].GetCo(k)));
			}
		}
		// add normals
		double a[3], b[3], x1,y1,z1,x2,y2,z2,x3,y3,z3;
		float c[3];
		vi[1] = mTopoFacets[fi].GetVertexInd(1);
		vi[2] = mTopoFacets[fi].GetVertexInd(2);
		x1 = mGeomVerts[vi[0]].GetCo(0);
		y1 = mGeomVerts[vi[0]].GetCo(1);
		z1 = mGeomVerts[vi[0]].GetCo(2);
		x2 = mGeomVerts[vi[1]].GetCo(0);
		y2 = mGeomVerts[vi[1]].GetCo(1);
		z2 = mGeomVerts[vi[1]].GetCo(2);
		x3 = mGeomVerts[vi[2]].GetCo(0);
		y3 = mGeomVerts[vi[2]].GetCo(1);
		z3 = mGeomVerts[vi[2]].GetCo(2);
		a[0] = x2 - x1; a[1] = y2 - y1; a[2] = z2 - z1;
		b[0] = x3 - x1; b[1] = y3 - y1; b[2] = z3 - z1;
		c[0] = float(a[1] * b[2] - a[2] * b[1]);
		c[1] = float(a[2] * b[0] - a[0] * b[2]);
		c[2] = float(a[0] * b[1] - a[1] * b[0]);
		Eigen::Vector3f cv(c[0], c[1], c[2]);
		cv.normalize();
		for (int i = 0; i < 3 * (n_v - 2); ++i)
		{
			renderNormals.push_back(cv[0]);
			renderNormals.push_back(cv[1]);
			renderNormals.push_back(cv[2]);
		}
	}

}

void Mesh::computeFaceEdgeCenters()
{
	vector<GeomVert> Ecenters;
	Ecenters.resize(mTopoEdges.size());
	for (int i = 0; i < mTopoEdges.size(); ++i)
	{
		int v0 = mTopoEdges[i].GetVertex(0);
		int v1 = mTopoEdges[i].GetVertex(1);
		Ecenters[i] = (mGeomVerts[v0] + mGeomVerts[v1]) / 2;
	}
	mFcenters.resize(mTopoFacets.size());
	mFEcenters.resize(mTopoFacets.size());
	for (int i = 0; i < mTopoFacets.size(); ++i)
	{
		GeomVert fcenter;
		int nv = mTopoFacets[i].GetNumberVertices();
		for (int j = 0; j < nv; ++j)
			fcenter += mGeomVerts[mTopoFacets[i].GetVertexInd(j)];
		mFcenters[i] = fcenter / nv;

		int ne = mTopoFacets[i].GetNumberEdges();
		mFEcenters[i].resize(ne);
		for (int j = 0; j < ne; ++j)
		{
			int ei = mTopoFacets[i].GetIncEdge(j);
			mFEcenters[i][j] = Ecenters[ei];
		}
	}
}

void Mesh::saveMesh()
{
	QFile file(QString("mesh.txt"));
	if (file.open(QIODevice::ReadWrite))
	{
		QTextStream stream(&file);
		int n_v = mGeomVerts.size();
		int n_f = mTopoFacets.size();
		stream << n_v << " " << n_f << endl;
		for (int i = 0; i < n_v; ++i)
		{
			stream << mGeomVerts[i].GetCo(0) << " "
				<< mGeomVerts[i].GetCo(1) << " "
				<< mGeomVerts[i].GetCo(2) << endl;
		}
		for (int i = 0; i < n_f; ++i)
		{
			int n = mTopoFacets[i].GetNumberVertices();
			stream << n;
			for (int j = 0; j < n; ++j)
			{
				stream << " " << mTopoFacets[i].GetVertexInd(j);
			}
			stream << endl;
		}
		file.close();
	}
}
// ------------------------------------------------------------

void Mesh::RevolveYaxis(const Eigen::MatrixXd & curve_pts, int n_curve_pts)
{
	Erase();
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

void Mesh::ExtrusionZaxis(const Eigen::MatrixXd & curve_pts, int n_curve_pts)
{
	Erase();
	Eigen::MatrixXd s0 = curve_pts.leftCols(n_curve_pts);
	Eigen::MatrixXd s1 = s0;
	Eigen::MatrixXd s2 = s0;
	double itv = m_depth / double (n_slice-1);
	for (int i = 1; i < n_slice; ++i)
	{
		for (int j = 0; j < n_curve_pts; ++j)
		{
			s2(2, j) += itv*i;
		}

		//add facet
		for (int k = 0; k < n_curve_pts - 1; ++k)
		{
			AddFacet(s1(0, k), s1(1, k), s1(2, k),
				s2(0, k), s2(1, k), s2(2, k),
				s2(0, k + 1), s2(1, k + 1), s2(2, k + 1),
				s1(0, k + 1), s1(1, k + 1), s1(2, k + 1));
		}
		s1 = s2;
	}
}

void Mesh::Sweep(const Eigen::MatrixXd & pts, int n_pts, const Eigen::MatrixXd & traj, int n_traj, bool closed)
{
	Erase();
	Eigen::MatrixXd s0 = pts.leftCols(n_pts);
	Eigen::MatrixXd s1, s2, s3;
	Eigen::Vector3d v1 = Eigen::Vector3d(0, 0, 1);
	Eigen::Vector3d v2 = traj.col(0) - traj.col(1);
	Eigen::Vector3d mean_p = pts.rowwise().mean();
	double theta = std::atan2((v1.cross(v2))[0],v1.dot(v2));
	double cos_theta = std::cos(theta);
	double sin_theta = std::sin(theta);

	Eigen::Matrix3d rot_mat;
	rot_mat << 1, 0, 0,
		0, cos_theta, -sin_theta,
		0, sin_theta, cos_theta;
	s1 = rot_mat * s0;
	Eigen::Vector3d offset = traj.col(0) - s1.rowwise().mean(); //rot_mat * mean_p;
	for (int i = 0; i < n_pts; ++i)
	{
		s1.col(i) += offset;
	}
	s3 = s1;

	for (int i = 1; i < n_traj; ++i)
	{
		v2 = traj.col(i-1) - traj.col(i);
		theta = std::atan2((v1.cross(v2))[0], v1.dot(v2));
		cos_theta = std::cos(theta);
		sin_theta = std::sin(theta);
		rot_mat(1, 1) = rot_mat(2, 2) = cos_theta;
		rot_mat(1, 2) = -sin_theta;
		rot_mat(2, 1) = sin_theta;
		s2 = rot_mat * s0;
		offset = traj.col(i) - s2.rowwise().mean();//rot_mat * mean_p;
		for (int ii = 0; ii < n_pts; ++ii)
		{
			s2.col(ii) += offset;
		}

		//add facet
		for (int k = 0; k < n_pts - 1; ++k)
		{
			AddFacet(s1(0, k), s1(1, k), s1(2, k),
				s1(0, k + 1), s1(1, k + 1), s1(2, k + 1),
				s2(0, k + 1), s2(1, k + 1), s2(2, k + 1),
				s2(0, k), s2(1, k), s2(2, k)
				);
		}
		s1 = s2;
	}


	if (closed)
	{
		s2 = s3;
		for (int k = 0; k < n_pts - 1; ++k)
		{
			AddFacet(s1(0, k), s1(1, k), s1(2, k),
				s1(0, k + 1), s1(1, k + 1), s1(2, k + 1),
				s2(0, k + 1), s2(1, k + 1), s2(2, k + 1),
				s2(0, k), s2(1, k), s2(2, k));
		}
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

