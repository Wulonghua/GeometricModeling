// Copyright (C) Geom Software e.U, Bernhard Kornberger, Graz/Austria
//
// This file is part of the Fade2D library. The student license is free
// of charge and covers personal non-commercial research. Licensees
// holding a commercial license may use this file in accordance with
// the Commercial License Agreement.
//
// This software is provided AS IS with NO WARRANTY OF ANY KIND,
// INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE.
//
// Please contact the author if any conditions of this licensing are
// not clear to you.
//
// Author: Bernhard Kornberger, bkorn (at) geom.at
// http://www.geom.at



#pragma once
#include "Point2.h"
#include "Segment2.h"
#include "Edge2.h"
#include <vector>

#if GEOM_PSEUDO3D==GEOM_TRUE
	namespace GEOM_FADE25D {
#elif GEOM_PSEUDO3D==GEOM_FALSE
	namespace GEOM_FADE2D {
#else
	#error GEOM_PSEUDO3D is not defined
#endif


/** \defgroup tools Tools
 *  @{
 */

#if GEOM_PSEUDO3D==GEOM_TRUE
/** \brief Get normal vector
 *
* Returns the normalized normal vector of the triangle defined by the
* three input points \p p0, \p p1, \p p2
*/
CLASS_DECLSPEC
Vector2 getNormalVector(const Point2& p0,const Point2& p1,const Point2& p2);
#endif



/** \brief Get directed edge
* The directed edges of \p vT are returned \p vDirectedEdgesOut.
* Directed means that each edge (a,b) with two adjacent triangles
* in vT is returned twice, as edge(a,b) and edge(b,a).
*/
CLASS_DECLSPEC
void getDirectedEdges(std::vector<Triangle2*>& vT,std::vector<Edge2>& vDirectedEdgesOut);

/** \brief Get undirected edges
 *
* A unique set of edges of \p vT is returned.
*/
CLASS_DECLSPEC
void getUndirectedEdges(std::vector<Triangle2*>& vT,std::vector<Edge2>& vUndirectedEdgesOut);






/** \brief Create polygons from a set of edges
 *
 *
 * A number of methods in Fade returns an unorganized set of edges that
 * delimit a certain area. But sometimes it is more beneficial to have
 * these edges organized as a set of one or more polygons. This is the
 * purpose of the present method.
 *
 * @param [in] vEdgesIn is a vector of oriented edges
 * @param [out] vvPolygonsOut contains one vector<Edge2> for each polygon
 * found in the input data.
 * @param [out] vRemainingOut is used to return the unused edges
 *
 * The present function adds for each polygon found in \p vEdgesIn one
 * vector<Edge2> to \p vvPolygonsOut. Edges that do not form a closed
 * polygon are collected in the last output vector.
 *
 * @note An Edge2 object represents an edge that is counterclockwise
 * (CCW) oriented with respect to its triangle. This orientation is
 * kept and therefore the output polygons can be both, CW or CCW-
 * oriented.
 *
 * \image html edges-to-polygons.png "Polygons created by edgesToPolygons"
 * \image latex edges-to-polygons.eps "Polygons created by edgesToPolygons" width=12cm
 *
 */
CLASS_DECLSPEC
void edgesToPolygons(
	std::vector<Edge2>& vEdgesIn,
	std::vector<std::vector<Edge2> >& vvPolygonsOut,
	std::vector<Edge2>& vRemainingOut
	);

/** \brief Get Borders
 *
 * Computes the border of the triangles in \p vT. The border consists
 * of all edges having only one adjacent triangle in vT.
 *
 * \param [in] vT are the input triangles
 * \param [out] vBorderSegmentsOut is used to return all border segments
*/
CLASS_DECLSPEC
void getBorders(const std::vector<Triangle2*>& vT,std::vector<Segment2>& vBorderSegmentsOut);
/** \brief Sort a vector of Segments
 *
 * The segments in vRing are reoriented and sorted such that subsequent
 * segments join at the endpoints.
*/
CLASS_DECLSPEC
bool sortRing(std::vector<Segment2>& vRing);

/** \brief Get the orientation of three points
 *
 * This function returns the \e exact orientation of the points \p p0, \p p1, \p p2
 * Possible values are \n
 * ORIENTATION2_COLLINEAR if \p p0, \p p1, \p p2 are located on a line, \n
 * ORIENTATION2_CCW if \p p0, \p p1, \p p2 are counterclockwise oriented \n
 * ORIENTATION2_CW if \p p0, \p p1, \p p2 are clockwise oriented \n
 *
 * Not thread-safe but a bit faster than the thread-safe version
*/

CLASS_DECLSPEC
Orientation2 getOrientation2(const Point2* p0,const Point2* p1,const Point2* p2);
/** \brief Get Orientation2 (MT)
 *
 * \see Orientation2 getOrientation2(const Point2* p0,const Point2* p1,const Point2* p2);
 *
 * This version is thread-safe.
*/

CLASS_DECLSPEC
Orientation2 getOrientation2_mt(const Point2* p0,const Point2* p1,const Point2* p2);

/// @private
CLASS_DECLSPEC
std::string getString(const Orientation2 ori);

/** @}*/


/** \defgroup codeInfo Version Information
 *  @{
 */

/** \brief Get the Fade2D version string
*/
CLASS_DECLSPEC
std::string getFade2DVersion();
/** \brief Get the major version number
*/
CLASS_DECLSPEC
int getMajorVersionNumber();
/** \brief Get the minor version number
*/
CLASS_DECLSPEC
int getMinorVersionNumber();
/** \brief Get the revision version number
*/
CLASS_DECLSPEC
int getRevisionNumber();
/** \brief Check if a RELEASE or a DEBUG version is used.
*/
CLASS_DECLSPEC
bool isRelease();


/// @private
CLASS_DECLSPEC
void setLic(
	const std::string& l1,
	const std::string& l2,
	const std::string& dt,
	const std::string& s1,
	const std::string& s2_
	);
/// @private
class Lic;
/// @private
Lic* getLic();

/** @}*/








/** \defgroup fileIO File I/O
 *  @{
 */


//////////////////////////
// READ AND WRITE, ASCII
//////////////////////////

/** \brief Write points to an ASCII file
 *
 * Writes points to an ASCII file,
 * \if SECTION_FADE25D
 * three coordinates (x y z) per line,
 * \else
 * two coordinates (x y) per line,
 * \endif
 * whitespace separated.
 *
 * \note Data exchange through ASCII files is easy and convenient but
 * floating point coordinates are not necessarily exact when represented
 * as decimal numbers. If the tiny rounding errors can't be accepted in
 * your setting you are advised to write binary files, (use
 * writePointsBIN(const char* filename,const std::vector<Point2>& vPointsIn))
 *
*/
CLASS_DECLSPEC
bool writePointsASCII(const char* filename,const std::vector<Point2*>& vPointsIn);

/** \brief Write points to an ASCII file
 *
 * \see bool writePointsASCII(const char* filename,const std::vector<Point2*>& vPointsIn);
*/
CLASS_DECLSPEC
bool writePointsASCII(const char* filename,const std::vector<Point2>& vPointsIn);





/** \brief Read (x y) points
 *
 * Reads points from an ASCII file. Expected file format: Two
 * coordinates (x y) per line, whitespace separated.
 *
 * \cond SECTION_FADE25D
 * The z coordinate is set to 0.
 * \endcond
*/
CLASS_DECLSPEC
bool readXY(const char* filename,std::vector<Point2>& vPointsOut);

#if GEOM_PSEUDO3D==GEOM_TRUE
// ONLY 2.5D
/** \brief Read (x y z) points
 *
 * Reads points from an ASCII file. Expected file format: Three
 * coordinates (x y z) per line, whitespace separated.
*/
CLASS_DECLSPEC
bool readXYZ(const char* filename,std::vector<Point2>& vPointsOut);
#endif


//////////////////////////
// READ AND WRITE, BINARY
//////////////////////////



/** \brief Write points to a binary file
 *
 * File format:\n
 * \if SECTION_FADE25D
 * int filetype (30)\n
 * \else
 * int filetype (20)\n
 * \endif
 * size_t numPoints (\p vPointsIn.size())\n
 * double x0\n
 * double y0\n
 * double z0\n
 * ...\n
 * double xn\n
 * double yn\n
 * double zn\n
*/
CLASS_DECLSPEC
bool writePointsBIN(const char* filename,std::vector<Point2>& vPointsIn);

/** \brief Write points to a binary file
 * \see writePointsBIN(const char* filename,const std::vector<Point2>& vPointsIn);
*/
CLASS_DECLSPEC
bool writePointsBIN(const char* filename,std::vector<Point2*>& vPointsIn);

/** \brief Read points from a binary file
 *
 * Reads points from a binary file of type 20 or 30
 * \see bool writePointsBIN(const char* filename,const std::vector<Point2>& vPointsIn)
*/
CLASS_DECLSPEC
bool readPointsBIN(const char* filename, std::vector<Point2>& vPointsIn);

/** \brief Write segments to a binary file
 *
 * Binary file format:\n
 * \if SECTION_FADE25D
 * int filetype (31) \n
 * \else
 * int filetype (21) \n
 * \endif
 * size_t numSegments (\p vSegmentsIn.size())  \n
 * double x0_source \n
 * double y0_source \n
 * \if SECTION_FADE25D
 * double z0_source \n
 * \endif
 * double x0_target \n
 * double y0_target \n
 * \if SECTION_FADE25D
 * double z0_target \n
 * \endif
 * ... \n
 * double xn_source \n
 * double yn_source \n
 * \if SECTION_FADE25D
 * double zn_source \n
 * \endif
 * double xn_target \n
 * double yn_target \n
 * \if SECTION_FADE25D
 * double zn_target \n
 * \endif
*/
CLASS_DECLSPEC
bool writeSegmentsBIN(const char* filename,std::vector<Segment2>& vSegmentsIn);

/** \brief Read segments from a binary file
 *
 * Reads segments from a binary file of type 21 or 31
 * \see bool writeSegmentsBIN(const char* filename,std::vector<Segment2>& vSegmentsIn);
*/
CLASS_DECLSPEC
bool readSegmentsBIN(const char* filename,std::vector<Segment2>& vSegmentsOut);






/** @}*/
} // NAMESPACE
