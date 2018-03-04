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
	QGLViewer(parent),m_picked(-1)
{
}


RenderWidget::~RenderWidget()
{
}

void RenderWidget::updateMesh()
{
	if (m_mesh->m_buildType == Mesh::REVOLUTION)
		m_mesh->RevolveYaxis(m_curve->m_points, m_curve->n_points);
	else if (m_mesh->m_buildType == Mesh::EXTRUSION)
		m_mesh->ExtrusionZaxis(m_curve->m_points, m_curve->n_points);
}

void RenderWidget::updateRender()
{
	m_program->bind();
	m_vertBuf.bind();
	m_vertBuf.allocate(sizeof(GLfloat) * 3 * (m_curve->n_ctls + m_curve->n_points));
	Eigen::Matrix3Xf  ctls = m_curve->m_ctls.leftCols(m_curve->n_ctls).cast<float>();
	Eigen::Matrix3Xf points = m_curve->m_points.leftCols(m_curve->n_points).cast<float>();
	m_vertBuf.write(0, ctls.data(), 3 * m_curve->n_ctls * sizeof(float));
	m_vertBuf.write(3 * m_curve->n_ctls * sizeof(float), points.data(), 3 * m_curve->n_points * sizeof(float));
	m_program->release();

	if (m_mesh->renderVerts.size() > 0)
	{
		m_s_program->bind();
		m_s_vertBuf.bind();
		m_s_vertBuf.allocate(sizeof(GLfloat) * (m_mesh->renderVerts.size() + m_mesh->renderNormals.size()));
		m_s_vertBuf.write(0, m_mesh->renderVerts.data(), m_mesh->renderVerts.size() * sizeof(float));
		m_s_vertBuf.write(m_mesh->renderVerts.size() * sizeof(float), m_mesh->renderNormals.data(),m_mesh->renderNormals.size()*sizeof(float));
		m_s_program->release();
	}
	update();
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
		updateMesh();
		updateRender();
	}
	else if(m_curve->m_contrlType == Curve::MOVE)
	{ 
		float winX = e->pos().x();
		float winY = e->pos().y();
		
		camera()->getModelViewMatrix(m_ModelView);
		camera()->getProjectionMatrix(m_Projection);
		camera()->getViewport(m_viewport);

		gluUnProject(winX, winY, 0.0, m_ModelView, m_Projection, m_viewport, &unproj_p0[0], &unproj_p0[1], &unproj_p0[2]);
		gluUnProject(winX, winY, 1.0, m_ModelView, m_Projection, m_viewport, &unproj_p1[0], &unproj_p1[1], &unproj_p1[2]);
		//auto tmp = this->camera()->position();
		//unproj_p0[0] = tmp.x;
		//unproj_p0[1] = tmp.y;
		//unproj_p0[2] = tmp.z;
		m_picked =  m_curve->pickControlPoint(unproj_p0, unproj_p1);
	}
	else if (m_curve->m_contrlType == Curve::VIEW)
	{
		QGLViewer::mousePressEvent(e);
	}
}

void RenderWidget::mouseMoveEvent(QMouseEvent * e)
{
	if (m_curve->m_contrlType == Curve::MOVE && m_picked > -1)
	{
		float winX = e->pos().x();
		float winY = e->pos().y();
		gluUnProject(winX, winY, 0.0, m_ModelView, m_Projection, m_viewport, &unproj_p0[0], &unproj_p0[1], &unproj_p0[2]);
		gluUnProject(winX, winY, 1.0, m_ModelView, m_Projection, m_viewport, &unproj_p1[0], &unproj_p1[1], &unproj_p1[2]);
		m_curve->setControlPoint(unproj_p0, unproj_p1, m_picked);
		//m_mesh->RevolveYaxis(m_curve->m_points, m_curve->n_points);
		updateMesh();
		updateRender();
	}
	else if (m_curve->m_contrlType == Curve::VIEW)
	{
		QGLViewer::mouseMoveEvent(e);
	}
}

void RenderWidget::mouseReleaseEvent(QMouseEvent * e)
{
	if(m_curve->m_contrlType == Curve::MOVE)
		m_picked = -1;
	else if (m_curve->m_contrlType == Curve::VIEW)
	{
		QGLViewer::mouseReleaseEvent(e);
	}

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
	m_vertBuf.create();	
	m_program->release();

	m_s_program = new QOpenGLShaderProgram(this);
	m_s_program->addShaderFromSourceFile(QOpenGLShader::Vertex, "./vshader_phong.glsl");
	m_s_program->addShaderFromSourceFile(QOpenGLShader::Fragment, "./fshader_phong.glsl");
	m_s_program->link();
	m_s_program->bind();
	m_s_vertBuf.create();
	m_s_program->release();
}

void RenderWidget::draw()
{
	camera()->getModelViewProjectionMatrix(m_mvpMat.data());
	camera()->getModelViewMatrix(m_mvMat.data());
	// draw curve with simple shading
	m_program->bind();
	m_program->setUniformValue("mvp", m_mvpMat);
	m_program->setUniformValue("col", QVector3D(1.0, 0.0, 0.0));
	
	m_vertBuf.bind();
	m_program->setAttributeBuffer(m_posAttr, GL_FLOAT, 0, 3);
	m_program->enableAttributeArray(m_posAttr);
	
	glPointSize(5);
	glDrawArrays(GL_POINTS, 0, m_curve->n_ctls);

	m_program->setUniformValue("col", QVector3D(0.0, 1.0, 0.0));
	glLineWidth(1);
	glDrawArrays(GL_LINE_STRIP, 0, m_curve->n_ctls);

	m_program->setUniformValue("col", QVector3D(1.0, 1.0, 0.0));
	glLineWidth(2);
	glDrawArrays(GL_LINE_STRIP, m_curve->n_ctls, m_curve->n_points);

	m_program->release();

	// draw surface with phong shading
	if (m_mesh->renderVerts.size() > 0) 
	{
		m_s_program->bind();
		m_s_program->setUniformValue("light_pos", m_light_pos);
		m_s_program->setUniformValue("ambient_coef", 0.8, 0.8, 0.8, 1.0);
		m_s_program->setUniformValue("diffuse_coef", 0.8, 0.8, 0.8, 1.0);
		m_s_program->setUniformValue("specular_coef", 0.8, 0.8, 0.8, 1.0);
		m_s_program->setUniformValue("light_ambient", 1.0, 1.0, 1.0, 1.0);
		m_s_program->setUniformValue("light_diffuse", 0.3, 0.3, 0.3, 1.0);
		m_s_program->setUniformValue("light_specular", 0.6, 0.6, 0.6, 1.0);
		m_s_program->setUniformValue("uVertexColor", 0.258, 0.525, 0.957, 1.0);
		
		m_s_program->setUniformValue("uMVPMatrix", m_mvpMat);
		m_s_program->setUniformValue("uMVMatrix", m_mvMat);
		m_nMat = m_mvMat.inverted().transposed();
		m_s_program->setUniformValue("uNMatrix", m_nMat);

		int vertexLocation = m_s_program->attributeLocation("aVertexPosition");
		int normalLocation = m_s_program->attributeLocation("aVertexNormal");

		m_s_vertBuf.bind();
		m_s_program->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3);
		m_s_program->enableAttributeArray(vertexLocation);

		m_s_program->setAttributeBuffer(normalLocation, GL_FLOAT, m_mesh->renderVerts.size()*sizeof(float), 3);
		m_s_program->enableAttributeArray(normalLocation);

		//glEnable(GL_BLEND);
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDrawArrays(GL_TRIANGLES, 0, m_mesh->renderVerts.size()/3);
		m_s_program->release();
	}

	//if (m_curve->m_contrlType==Curve::VIEW || m_curve->m_contrlType==Curve::MOVE)
	//{
	//	glBegin(GL_LINES);
	//	glColor3f(0, 0, 1.0);
	//	glVertex3d(unproj_p0[0], unproj_p0[1], unproj_p0[2]);
	//	glVertex3d(unproj_p1[0], unproj_p1[1], unproj_p1[2]);
	//	glEnd();
	//}
}

void RenderWidget::postDraw()
{
	QGLViewer::postDraw();
}