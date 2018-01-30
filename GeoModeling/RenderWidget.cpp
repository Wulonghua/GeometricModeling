#include "RenderWidget.h"



RenderWidget::RenderWidget(QWidget *parent):
	QGLViewer(parent)
{
}


RenderWidget::~RenderWidget()
{
}

void RenderWidget::keyPressEvent(QKeyEvent *e) 
{
	QGLViewer::keyPressEvent(e);
}

void RenderWidget::mousePressEvent(QMouseEvent * e)
{
	camera()->getModelViewMatrix(m_ModelView);
	camera()->getProjectionMatrix(m_Projection);
	camera()->getViewport(m_viewport);

	float winX = e->pos().x();
	float winY = e->pos().y();

	gluUnProject(winX, winY, 0.0, m_ModelView, m_Projection, m_viewport, &unproj_p0[0], &unproj_p0[1], &unproj_p0[2]);
	gluUnProject(winX, winY, 1.0, m_ModelView, m_Projection, m_viewport, &unproj_p1[0], &unproj_p1[1], &unproj_p1[2]);
	std::cout << unproj_p0[0] << " " << unproj_p0[1] << " " << unproj_p0[2] << std::endl;
	std::cout << unproj_p1[0] << " " << unproj_p1[1] << " " << unproj_p1[2] << std::endl;
	//QGLViewer::mousePressEvent(e);
}

void RenderWidget::mouseMoveEvent(QMouseEvent * e)
{
	//QGLViewer::mouseMoveEvent(e);
}

void RenderWidget::mouseReleaseEvent(QMouseEvent * e)
{
	//QGLViewer::mouseReleaseEvent(e);
}

void RenderWidget::init()
{
	restoreStateFromFile();
}

void RenderWidget::draw()
{
	//const float nbSteps = 200.0;

	//glBegin(GL_QUAD_STRIP);
	//for (int i = 0; i < nbSteps; ++i) {
	//	const float ratio = i / nbSteps;
	//	const float angle = 21.0 * ratio;
	//	const float c = cos(angle);
	//	const float s = sin(angle);
	//	const float r1 = 1.0 - 0.8f * ratio;
	//	const float r2 = 0.8f - 0.8f * ratio;
	//	const float alt = ratio - 0.5f;
	//	const float nor = 0.5f;
	//	const float up = sqrt(1.0 - nor * nor);
	//	glColor3f(1.0 - ratio, 0.2f, ratio);
	//	glNormal3f(nor * c, up, nor * s);
	//	glVertex3f(r1 * c, alt, r1 * s);
	//	glVertex3f(r2 * c, alt + 0.05f, r2 * s);
	//}
	//glEnd();
}

void RenderWidget::postDraw()
{
	QGLViewer::postDraw();
}