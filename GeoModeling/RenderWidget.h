#pragma once

#include "QGLViewer\qglviewer.h"
#include "Curve.h"
#include <iostream>
#include <memory>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>

class RenderWidget :
	public QGLViewer, protected QOpenGLFunctions
{
public:
	RenderWidget(QWidget *parent);
	~RenderWidget();

	void setCurve(std::shared_ptr<Curve>& curve) { m_curve = curve; }
	void updateRender();

protected:
	virtual void draw();
	virtual void init();
	virtual void postDraw();
	virtual void keyPressEvent(QKeyEvent *e);
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);

private:
	void transformPoint(GLdouble out[4], const GLdouble m[16], const GLdouble in[4]);
	GLint project(GLdouble objx, GLdouble objy, GLdouble objz,
		const GLdouble model[16], const GLdouble proj[16],
		const GLint viewport[4],
		GLdouble * winx, GLdouble * winy, GLdouble * winz);

	double m_ModelView[16];
	double m_Projection[16];
	int	   m_viewport[4];
	double unproj_p0[3];
	double unproj_p1[3];

	GLuint m_posAttr;

	QMatrix4x4 m_mvpMat;
	QOpenGLBuffer m_vertBuf;
	QOpenGLShaderProgram *m_program;


	std::shared_ptr<Curve> m_curve;
	int m_picked;
};

inline void RenderWidget::transformPoint(GLdouble out[4], const GLdouble m[16], const GLdouble in[4])
{
#define M(row,col)  m[col*4+row]
	out[0] =
		M(0, 0) * in[0] + M(0, 1) * in[1] + M(0, 2) * in[2] + M(0, 3) * in[3];
	out[1] =
		M(1, 0) * in[0] + M(1, 1) * in[1] + M(1, 2) * in[2] + M(1, 3) * in[3];
	out[2] =
		M(2, 0) * in[0] + M(2, 1) * in[1] + M(2, 2) * in[2] + M(2, 3) * in[3];
	out[3] =
		M(3, 0) * in[0] + M(3, 1) * in[1] + M(3, 2) * in[2] + M(3, 3) * in[3];
#undef M
}

inline GLint RenderWidget::project(GLdouble objx, GLdouble objy, GLdouble objz,
	const GLdouble model[16], const GLdouble proj[16],
	const GLint viewport[4],
	GLdouble * winx, GLdouble * winy, GLdouble * winz)
{
	GLdouble in[4], out[4];

	in[0] = objx;
	in[1] = objy;
	in[2] = objz;
	in[3] = 1.0;
	transformPoint(out, model, in);
	transformPoint(in, proj, out);

	if (in[3] == 0.0)
		return GL_FALSE;

	in[0] /= in[3];
	in[1] /= in[3];
	in[2] /= in[3];

	*winx = viewport[0] + (1 + in[0]) * viewport[2] / 2;
	*winy = viewport[1] + (1 + in[1]) * viewport[3] / 2;

	*winz = (1 + in[2]) / 2;
	return GL_TRUE;
}