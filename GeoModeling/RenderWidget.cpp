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
	QGLViewer(parent),m_picked(-1), m_plane(0)
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
	else if(m_mesh->m_buildType == Mesh::SWEEP)
		m_mesh->Sweep(m_curve->m_points, m_curve->n_points, m_traj->m_points, m_traj->n_points, m_traj->m_closed);
}

void RenderWidget::setLookatPlane()
{
	if (m_plane == 0) //xy-plane
	{
		camera()->setPosition(qglviewer::Vec(0, 0, 3));
		camera()->setUpVector(qglviewer::Vec(0, 1, 0));
	}
	else if (m_plane == 1) //yz-plane
	{
		camera()->setPosition(qglviewer::Vec(3, 0, 0));
		camera()->setUpVector(qglviewer::Vec(0, 0, 1));
	}
	
	camera()->lookAt(sceneCenter());
	update();
}

void RenderWidget::updateRender()
{
	m_program->bind();
	if (m_curve->n_ctlsb==0)
	{
		m_vertBuf.bind();
		m_vertBuf.allocate(sizeof(GLfloat) * 3 * (m_curve->n_ctls + m_curve->n_points));
		Eigen::Matrix3Xf  ctls = m_curve->m_ctls.leftCols(m_curve->n_ctls).cast<float>();
		Eigen::Matrix3Xf points = m_curve->m_points.leftCols(m_curve->n_points).cast<float>();
		m_vertBuf.write(0, ctls.data(), 3 * m_curve->n_ctls * sizeof(float));
		m_vertBuf.write(3 * m_curve->n_ctls * sizeof(float), points.data(), 3 * m_curve->n_points * sizeof(float));

		m_t_vertBuf.bind();
		m_t_vertBuf.allocate(sizeof(GLfloat) * 3 * (m_traj->n_ctls + m_traj->n_points));
		ctls = m_traj->m_ctls.leftCols(m_traj->n_ctls).cast<float>();
		points = m_traj->m_points.leftCols(m_traj->n_points).cast<float>();
		m_t_vertBuf.write(0, ctls.data(), 3 * m_traj->n_ctls * sizeof(float));
		m_t_vertBuf.write(3 * m_traj->n_ctls * sizeof(float), points.data(), 3 * m_traj->n_points * sizeof(float));
	}
	else
	{
		m_vertBuf.bind();
		int n = m_curve->n_ctls*m_curve->n_ctlsb;
		m_vertBuf.allocate(sizeof(GLfloat) * 3 * n);
		Eigen::Matrix3Xf  ctls = m_curve->m_ctls.leftCols(n).cast<float>();
		m_vertBuf.write(0, ctls.data(), 3 * n * sizeof(float));
	}

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
	if (m_plane==0) // xy-plane, generator curve
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
		if (m_curve->m_contrlType == Curve::DEL)
		{
			camera()->getModelViewMatrix(m_ModelView);
			camera()->getProjectionMatrix(m_Projection);
			camera()->getViewport(m_viewport);

			float winX = e->pos().x();
			float winY = e->pos().y();

			gluUnProject(winX, winY, 0.0, m_ModelView, m_Projection, m_viewport, &unproj_p0[0], &unproj_p0[1], &unproj_p0[2]);
			gluUnProject(winX, winY, 1.0, m_ModelView, m_Projection, m_viewport, &unproj_p1[0], &unproj_p1[1], &unproj_p1[2]);
			
			m_picked = m_curve->pickControlPoint(unproj_p0, unproj_p1);
			if (m_picked > -1 && QMessageBox::Yes == QMessageBox::question(this,
				tr("Delete Point"),
				tr("Confirm to delete this chosen point?")))
			{
				m_curve->deleteControlPoint(m_picked);
				m_curve->generateCurves();
				updateMesh();
				updateRender();
			}
		}
		else if (m_curve->m_contrlType == Curve::MOVE)
		{
			float winX = e->pos().x();
			float winY = e->pos().y();

			camera()->getModelViewMatrix(m_ModelView);
			camera()->getProjectionMatrix(m_Projection);
			camera()->getViewport(m_viewport);

			gluUnProject(winX, winY, 0.0, m_ModelView, m_Projection, m_viewport, &unproj_p0[0], &unproj_p0[1], &unproj_p0[2]);
			gluUnProject(winX, winY, 1.0, m_ModelView, m_Projection, m_viewport, &unproj_p1[0], &unproj_p1[1], &unproj_p1[2]);
			m_picked = m_curve->pickControlPoint(unproj_p0, unproj_p1);
		}
		else if (m_curve->m_contrlType == Curve::VIEW)
		{
			QGLViewer::mousePressEvent(e);
		}
	}
	else if (m_plane == 1)
	{
		if (m_traj->m_contrlType == Curve::ADD)
		{
			camera()->getModelViewMatrix(m_ModelView);
			camera()->getProjectionMatrix(m_Projection);
			camera()->getViewport(m_viewport);

			float winX = e->pos().x();
			float winY = e->pos().y();

			gluUnProject(winX, winY, 0.0, m_ModelView, m_Projection, m_viewport, &unproj_p0[0], &unproj_p0[1], &unproj_p0[2]);
			gluUnProject(winX, winY, 1.0, m_ModelView, m_Projection, m_viewport, &unproj_p1[0], &unproj_p1[1], &unproj_p1[2]);

			m_traj->addControlPoints(unproj_p0, unproj_p1,m_plane);
			updateMesh();
			updateRender();
		}
		else if (m_traj->m_contrlType == Curve::MOVE)
		{
			float winX = e->pos().x();
			float winY = e->pos().y();

			camera()->getModelViewMatrix(m_ModelView);
			camera()->getProjectionMatrix(m_Projection);
			camera()->getViewport(m_viewport);

			gluUnProject(winX, winY, 0.0, m_ModelView, m_Projection, m_viewport, &unproj_p0[0], &unproj_p0[1], &unproj_p0[2]);
			gluUnProject(winX, winY, 1.0, m_ModelView, m_Projection, m_viewport, &unproj_p1[0], &unproj_p1[1], &unproj_p1[2]);
			m_picked = m_traj->pickControlPoint(unproj_p0, unproj_p1);
		}
		else if (m_traj->m_contrlType == Curve::VIEW)
		{
			QGLViewer::mousePressEvent(e);
		}

	}
}

void RenderWidget::mouseMoveEvent(QMouseEvent * e)
{
	if (m_plane==0)
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
	else if (m_plane == 1)
	{
		if (m_traj->m_contrlType == Curve::MOVE && m_picked > -1)
		{
			float winX = e->pos().x();
			float winY = e->pos().y();
			gluUnProject(winX, winY, 0.0, m_ModelView, m_Projection, m_viewport, &unproj_p0[0], &unproj_p0[1], &unproj_p0[2]);
			gluUnProject(winX, winY, 1.0, m_ModelView, m_Projection, m_viewport, &unproj_p1[0], &unproj_p1[1], &unproj_p1[2]);
			m_traj->setControlPoint(unproj_p0, unproj_p1, m_picked,m_plane);
			//m_mesh->RevolveYaxis(m_curve->m_points, m_curve->n_points);
			updateMesh();
			updateRender();
		}
		else if (m_traj->m_contrlType == Curve::VIEW)
		{
			QGLViewer::mouseMoveEvent(e);
		}
	}
}

void RenderWidget::mouseReleaseEvent(QMouseEvent * e)
{
	if (m_plane==0)
	{
		if (m_curve->m_contrlType == Curve::MOVE)
			m_picked = -1;
		else if (m_curve->m_contrlType == Curve::VIEW)
		{
			QGLViewer::mouseReleaseEvent(e);
		}
	}
	else if (m_plane == 1)
	{
		if (m_traj->m_contrlType == Curve::MOVE)
			m_picked = -1;
		else if (m_traj->m_contrlType == Curve::VIEW)
		{
			QGLViewer::mouseReleaseEvent(e);
		}
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
	m_t_vertBuf.create();
	m_program->release();

	m_s_program = new QOpenGLShaderProgram(this);
	m_s_program->addShaderFromSourceFile(QOpenGLShader::Vertex, "./vshader_phong.glsl");
	m_s_program->addShaderFromSourceFile(QOpenGLShader::Fragment, "./fshader_phong.glsl");
	m_s_program->link();
	m_s_program->bind();
	m_s_vertBuf.create();
	m_s_program->release();

	camera()->setPosition(qglviewer::Vec(0, 0, 3));
	camera()->setUpVector(qglviewer::Vec(0, 1, 0));
	camera()->lookAt(sceneCenter());
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
	if (m_curve->n_ctlsb == 0)
	{
		glPointSize(5);
		glDrawArrays(GL_POINTS, 0, m_curve->n_ctls);

		m_program->setUniformValue("col", QVector3D(0.0, 1.0, 0.0));
		glLineWidth(1);
		if (m_curve->m_closed)
			glDrawArrays(GL_LINE_LOOP, 0, m_curve->n_ctls);
		else
			glDrawArrays(GL_LINE_STRIP, 0, m_curve->n_ctls);

		m_program->setUniformValue("col", QVector3D(1.0, 1.0, 0.0));
		glLineWidth(2);
		glDrawArrays(GL_LINE_STRIP, m_curve->n_ctls, m_curve->n_points);

		/***************************************************************/
	}
	else // contol polygon
	{
		glPointSize(5);
		glDrawArrays(GL_POINTS, 0, m_curve->n_ctls*m_curve->n_ctlsb);
	}

	/**********************trajactory cuve**************************/
	m_program->setUniformValue("col", QVector3D(1.0, 0.0, 0.0));
	m_t_vertBuf.bind();
	m_program->setAttributeBuffer(m_posAttr, GL_FLOAT, 0, 3);
	m_program->enableAttributeArray(m_posAttr);

	glPointSize(5);
	glDrawArrays(GL_POINTS, 0, m_traj->n_ctls);

	m_program->setUniformValue("col", QVector3D(0.0, 1.0, 0.0));
	glLineWidth(1);
	if (m_traj->m_closed)
		glDrawArrays(GL_LINE_LOOP, 0, m_traj->n_ctls);
	else
		glDrawArrays(GL_LINE_STRIP, 0, m_traj->n_ctls);

	m_program->setUniformValue("col", QVector3D(1.0, 1.0, 0.0));
	glLineWidth(2);
	glDrawArrays(GL_LINE_STRIP, m_traj->n_ctls, m_traj->n_points);
	
	/***************************************************************/
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
		//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawArrays(GL_TRIANGLES, 0, m_mesh->renderVerts.size()/3);
		//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
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