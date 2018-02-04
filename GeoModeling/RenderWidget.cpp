#include "RenderWidget.h"

static const char *vertexShaderSource =
"attribute highp vec3 posAttr;\n"
"uniform highp mat4 mvp;\n"
"void main() {\n"
"   gl_Position = mvp * vec4(posAttr,1);\n"
"}\n";

static const char *fragmentShaderSource =
"uniform highp vec3 col;\n"
"void main() {\n"
"   gl_FragColor = vec4(col,1);\n"
"}\n";

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
	if (m_curve->m_contrlType == Curve::ADD)
	{
		camera()->getModelViewMatrix(m_ModelView);
		camera()->getProjectionMatrix(m_Projection);
		camera()->getViewport(m_viewport);

		float winX = e->pos().x();
		float winY = e->pos().y();

		gluUnProject(winX, winY, 0.0, m_ModelView, m_Projection, m_viewport, &unproj_p0[0], &unproj_p0[1], &unproj_p0[2]);
		gluUnProject(winX, winY, 1.0, m_ModelView, m_Projection, m_viewport, &unproj_p1[0], &unproj_p1[1], &unproj_p1[2]);

		m_curve->addControlPoints(unproj_p0, unproj_p1);
		m_ctlBuf.bind();
		Eigen::Matrix3Xf  ctls = m_curve->m_ctls.cast<float>();
		m_ctlBuf.allocate(ctls.data(), sizeof(GLfloat) * 3 * m_curve->n_ctls);
		update();
		//std::cout << m_curve->m_ctls.col(m_curve->n_ctls - 1) << std::endl;
		//std::cout << unproj_p0[0] << " " << unproj_p0[1] << " " << unproj_p0[2] << std::endl;
		//std::cout << unproj_p1[0] << " " << unproj_p1[1] << " " << unproj_p1[2] << std::endl;
	}
	else if(m_curve->m_contrlType == Curve::VIEW)
	{ 
		QGLViewer::mousePressEvent(e);
	}
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
	initializeOpenGLFunctions();
	m_program = new QOpenGLShaderProgram(this);
	m_program->addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShaderSource);
	m_program->addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShaderSource);
	m_program->link();
	m_program->bind();

	m_posAttr = m_program->attributeLocation("posAttr");
	//m_colUniform = m_program->attributeLocation("col");
	//m_mvpUniform = m_program->uniformLocation("mvp");

	m_ctlBuf.create();
	//m_ctlBuf.bind();
	//if (m_curve && m_curve->n_ctls > 0) 
	//{
	//	Eigen::Matrix3Xf  ctls = m_curve->m_ctls.cast<float>();
	//	m_ctlBuf.allocate(ctls.data(), sizeof(GLfloat) * 3 * m_curve->n_ctls);
	//}
		
	
}

void RenderWidget::draw()
{
	m_program->bind();
	camera()->getModelViewProjectionMatrix(m_mvpMat.data());
	m_program->setUniformValue("mvp", m_mvpMat);
	m_program->setUniformValue("col", QVector3D(1.0, 0.0, 0.0));

	m_ctlBuf.bind();
	m_program->setAttributeBuffer(m_posAttr, GL_FLOAT, 0, 3);
	m_program->enableAttributeArray(m_posAttr);

	glPointSize(5);
	glDrawArrays(GL_POINTS, 0, m_curve->n_ctls);
	

	m_program->setUniformValue("col", QVector3D(1.0, 1.0, 0.0));
	glLineWidth(2);
	glDrawArrays(GL_LINE_STRIP, 0, m_curve->n_ctls);

	m_program->release();

}

void RenderWidget::postDraw()
{
	QGLViewer::postDraw();
}