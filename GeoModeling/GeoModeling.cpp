#include "GeoModeling.h"

GeoModeling::GeoModeling(QWidget *parent)
	: QMainWindow(parent)
{
	QWidget *dockw = new QWidget();
	ui_control.setupUi(dockw);
	QDockWidget *dock = new QDockWidget();
	dock->setWidget(dockw);
	dock->setMinimumWidth(250);
	this->addDockWidget(Qt::LeftDockWidgetArea, dock);
	ui.setupUi(this);

	m_curve = std::make_shared<Curve>();
	m_traj = std::make_shared<Curve>();
	m_mesh = std::make_shared<Mesh>();
	ui.glWidget->setCurve(m_curve, m_traj);
	ui.glWidget->setMesh(m_mesh);
	m_curve->reset();
	m_traj->reset();
	m_mesh->reset();
	initConnections();
}

void GeoModeling::changeControlState()
{
	if (ui_control.radioButton_addPoints->isChecked())
	{
		m_curve->m_contrlType = Curve::ADD;
		m_traj->m_contrlType = Curve::ADD;
	}
	if (ui_control.radioButton_delPoints->isChecked())
	{
		m_curve->m_contrlType = Curve::DEL;
		m_traj->m_contrlType = Curve::DEL;
	}
	else if (ui_control.radioButton_movePoints->isChecked())
	{
		m_curve->m_contrlType = Curve::MOVE;
		m_traj->m_contrlType = Curve::MOVE;
	}
	else if (ui_control.radioButton_View->isChecked())
	{
		m_curve->m_contrlType = Curve::VIEW;
		m_traj->m_contrlType = Curve::VIEW;
	}
}

void GeoModeling::changeCurveType()
{
	std::shared_ptr<Curve> curve;
	if(ui.glWidget->m_plane==0)
		curve = m_curve;
	else if(ui.glWidget->m_plane == 1)
		curve = m_traj;
	curve->changeCloseStatus(ui_control.checkBox_closed->isChecked());
	if (ui_control.radioButton_Sampling->isChecked())
	{
		curve->m_genType = Curve::SAMPLE;
		ui_control.label_precision->setText(QString("Precisions: "));
		curve->n_precision = ui_control.spinBox_precision->value();
	}
	else if (ui_control.radioButton_Subdivision->isChecked())
	{
		curve->m_genType = Curve::SUBDIVISION;
		ui_control.label_precision->setText(QString("Iterations: "));
		curve->n_precision = ui_control.spinBox_precision->value();
	}

	if (ui_control.radioButton_Bezier->isChecked())
	{
		if (curve->m_closed && curve->m_curveType != Curve::Bezier && curve->n_ctls>0)
		{
			curve->m_ctls.col(curve->n_ctls) = curve->m_ctls.col(0);
			curve->n_ctls += 1;
		}
		curve->m_curveType = Curve::Bezier;
		curve->generateBezierPoints();
	}
	else if (ui_control.radioButton_Quadric_Bspline->isChecked())
	{
		if (curve->m_closed && curve->m_curveType == Curve::Bezier && curve->n_ctls>0)
			curve->n_ctls -= 1;
		curve->m_curveType = Curve::Quadric_B_spline;
		curve->generateQuadBspline();
	}
	else if (ui_control.radioButton_Cubic_Bspline->isChecked())
	{
		if (curve->m_closed && curve->m_curveType == Curve::Bezier && curve->n_ctls>0)
			curve->n_ctls -= 1;
		curve->m_curveType = Curve::Cubic_B_spline;
		curve->generateCubicBspline();
	}
	ui.glWidget->updateMesh();
	ui.glWidget->updateRender();
}

void GeoModeling::changePlane()
{
	if (ui_control.radioButton_generator->isChecked())
	{
		ui.glWidget->m_plane = 0; // xy-plane
	}
	else if (ui_control.radioButton_trajactory->isChecked())
	{
		ui.glWidget->m_plane = 1; // yz-plane
	}
	ui.glWidget->setLookatPlane();

	std::shared_ptr<Curve> curve;
	if (ui.glWidget->m_plane == 0)
		curve = m_curve;
	else if (ui.glWidget->m_plane == 1)
		curve = m_traj;

	ui_control.checkBox_closed->setChecked(curve->m_closed);
	if (curve->m_curveType == Curve::Bezier)
	{
		ui_control.radioButton_Quadric_Bspline->setChecked(false);
		ui_control.radioButton_Cubic_Bspline->setChecked(false);
		ui_control.radioButton_Bezier->setChecked(true);
	}
	else if (curve->m_curveType == Curve::Quadric_B_spline)
	{
		ui_control.radioButton_Bezier->setChecked(false);
		ui_control.radioButton_Cubic_Bspline->setChecked(false);
		ui_control.radioButton_Quadric_Bspline->setChecked(true);
	}	
	else if (curve->m_curveType == Curve::Cubic_B_spline)
	{
		ui_control.radioButton_Bezier->setChecked(false);
		ui_control.radioButton_Quadric_Bspline->setChecked(false);
		ui_control.radioButton_Cubic_Bspline->setChecked(true);
	}
		

	ui_control.spinBox_precision->setValue(curve->n_precision);

	if (curve->m_genType == Curve::SAMPLE)
		ui_control.radioButton_Sampling->setChecked(true);
	else if (curve->m_genType == Curve::SUBDIVISION)
		ui_control.radioButton_Subdivision->setChecked(true);
}

void GeoModeling::changePrecision(int x)
{
	std::shared_ptr<Curve> curve;
	if (ui.glWidget->m_plane == 0)
		curve = m_curve;
	else if (ui.glWidget->m_plane == 1)
		curve = m_traj;

	curve->n_precision = x;
	if (ui_control.radioButton_Bezier->isChecked())
	{
		curve->m_curveType = Curve::Bezier;
		curve->generateBezierPoints();
	}
	else if (ui_control.radioButton_Quadric_Bspline->isChecked())
	{
		curve->m_curveType = Curve::Quadric_B_spline;
		curve->generateQuadBspline();
	}
	else if (ui_control.radioButton_Cubic_Bspline->isChecked())
	{
		curve->m_curveType = Curve::Cubic_B_spline;
		curve->generateCubicBspline();
	}
	if (m_mesh->GetNumberFacets() > 0)
	{
		ui.glWidget->updateMesh();
	}
	ui.glWidget->updateRender();
}

void GeoModeling::changeSlices(int x)
{
	m_mesh->n_slice = x;
	//m_mesh->reset();
	//m_mesh->RevolveYaxis(m_curve->m_points, m_curve->n_points);
	ui.glWidget->updateMesh();
	ui.glWidget->updateRender();
}

void GeoModeling::changeZdepth(double x)
{
	m_mesh->m_depth = x;
	if (m_mesh->m_buildType == Mesh::EXTRUSION)
	{
		m_mesh->ExtrusionZaxis(m_curve->m_points, m_curve->n_points);
		ui.glWidget->updateRender();
	}
}

void GeoModeling::clearState()
{
	m_curve->reset();
	m_traj->reset();
	m_mesh->reset();
	ui.glWidget->updateRender();
}

void GeoModeling::doYRevolution()
{
	// just for tri-mesh test
	//m_mesh->AddFacet(0, 0, 0, 1, 0, 0, 0, 1, 0);
	//m_mesh->AddFacet(1, 0, 0, 1, 1, 0, 0, 1, 0);
	//m_mesh->AddFacet(1, 0, 0, 2, 0, 0, 1, 1, 0);
	//m_mesh->AddFacet(0, 1, 0, 1, 1, 0, 0, 2, 0);

	// for quad mesh test
	//m_mesh->AddFacet(0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0);
	//m_mesh->AddFacet(1, 0, 0, 2, 0, 0, 2, 1, 0, 1, 1, 0);
	//m_mesh->AddFacet(0, 1, 0, 1, 1, 0, 1, 2, 0, 0, 2, 0);
	//m_mesh->AddFacet(1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0);

	//std::cout << std::endl;
	if (m_curve->n_points < 1)
	{
		QMessageBox msgBox;
		msgBox.setText("There is no curve to generate revolution surface.");
		msgBox.exec();
		return;
	}

	if (!m_curve->isValidYrevolve())
	{
		QMessageBox msgBox;
		msgBox.setText("Warning: curve overlaps with Y axis");
		msgBox.exec();
		return;
	}
	m_mesh->m_buildType = Mesh::REVOLUTION;
	m_mesh->RevolveYaxis(m_curve->m_points, m_curve->n_points);
	m_curve->m_contrlType = Curve::VIEW;
	ui_control.radioButton_View->setChecked(true);
	ui.glWidget->updateRender();
	//std::cout << m_mesh->GetNumberVertices() << " " << m_mesh->GetNumberEdges() << " " << m_mesh->GetNumberFacets() << std::endl;

}

void GeoModeling::doZExtrusion()
{
	if (m_curve->n_points < 1)
	{
		QMessageBox msgBox;
		msgBox.setText("There is no curve to generate extrusion surface.");
		msgBox.exec();
		return;
	}
	m_mesh->m_buildType = Mesh::EXTRUSION;
	m_mesh->ExtrusionZaxis(m_curve->m_points, m_curve->n_points);
	m_curve->m_contrlType = Curve::VIEW;
	ui_control.radioButton_View->setChecked(true);
	ui.glWidget->updateRender();
}

void GeoModeling::doSweep()
{
	if (m_curve->n_points < 1)
	{
		QMessageBox msgBox;
		msgBox.setText("There is no generator curve to generate sweeping surface.");
		msgBox.exec();
		return;
	}
	if (m_traj->n_points < 1)
	{
		QMessageBox msgBox;
		msgBox.setText("There is no trajactory curve to generate sweeping surface.");
		msgBox.exec();
		return;
	}

	m_mesh->m_buildType = Mesh::SWEEP;
	m_mesh->Sweep(m_curve->m_points, m_curve->n_points,m_traj->m_points,m_traj->n_points, m_traj->m_closed);
	m_curve->m_contrlType = Curve::VIEW;
	m_traj->m_contrlType = Curve::VIEW;
	ui_control.radioButton_View->setChecked(true);
	ui.glWidget->updateRender();
}

void GeoModeling::doSaveMesh()
{
	if (m_mesh->GetNumberFacets() < 1)
	{
		QMessageBox msgBox;
		msgBox.setText("There is no Mesh.");
		msgBox.exec();
		return;
	}
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Mesh"),
		"./untitled.off",
		tr("Mesh (*.off)"));

	m_mesh->saveMesh(filename);
}

void GeoModeling::doLoadMesh()
{
	QString filename = QFileDialog::getOpenFileName(this,
		tr("Load Mesh"), "./", tr("Model File (*.off)"));
	m_mesh->LoadModel(filename);
	ui.glWidget->updateRender();
	m_curve->m_contrlType = Curve::VIEW;
	m_traj->m_contrlType = Curve::VIEW;
	ui_control.radioButton_View->setChecked(true);
}

void GeoModeling::doSubDooSabin()
{
	std::shared_ptr<Mesh> mesh=std::make_shared<Mesh>();
	m_mesh->SubDooSabin(mesh);
	m_mesh = mesh;
	ui.glWidget->setMesh(m_mesh);
	m_mesh->prepareRender();
	ui.glWidget->updateRender();
}

void GeoModeling::doCatmullClark()
{
	std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>();
	m_mesh->SubCatmullClark(mesh);
	m_mesh = mesh;
	ui.glWidget->setMesh(m_mesh);
	m_mesh->prepareRender();
	ui.glWidget->updateRender();
}

void GeoModeling::doLoop()
{
	std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>();
	m_mesh->SubLoop(mesh);
	m_mesh = mesh;
	ui.glWidget->setMesh(m_mesh);
	m_mesh->prepareRender();
	ui.glWidget->updateRender();
}

void GeoModeling::doGenCtlPoly()
{
	int n_slices = ui_control.spinBox_slices->value();
	m_curve->generateControlPolygon(n_slices);
	ui.glWidget->updateRender();
	m_curve->m_contrlType = Curve::VIEW;
	m_traj->m_contrlType = Curve::VIEW;
	ui_control.radioButton_View->setChecked(true);
}

void GeoModeling::doGenBezierSurface()
{
	if (m_curve->n_ctlsb < 3)
	{
		std::cerr << "Does not have valid control polygon" << std::endl;
		return;
	}
	m_mesh->reset();
	m_curve->generateBezierSurface(m_mesh);
	m_mesh->prepareRender();
	ui.glWidget->updateRender();
	std::cout << "Done Generating Bezier Surface" << std::endl;
}

void GeoModeling::doGenCubicSplineSurface()
{
	if (m_curve->n_ctlsb < 4)
	{
		std::cerr << "Does not have valid control polygon" << std::endl;
		return;
	}
	m_mesh->reset();
	m_curve->generateCubicSplineSurface(m_mesh);
	m_mesh->prepareRender();
	ui.glWidget->updateRender();
	std::cout << "Done Generating Bezier Surface" << std::endl;
}

void GeoModeling::doNNCrust()
{
	m_curve->NNCrust();
	m_curve->generateCurves();
	ui.glWidget->updateRender();
}

void GeoModeling::doCrust()
{
	m_curve->Crust();
	m_curve->generateCurves();
	ui.glWidget->updateRender();
}

void GeoModeling::initConnections()
{
	connect(ui_control.radioButton_addPoints, &QRadioButton::clicked, this, &GeoModeling::changeControlState);
	connect(ui_control.radioButton_delPoints, &QRadioButton::clicked, this, &GeoModeling::changeControlState);
	connect(ui_control.radioButton_movePoints, &QRadioButton::clicked, this, &GeoModeling::changeControlState);
	connect(ui_control.radioButton_View, &QRadioButton::clicked, this, &GeoModeling::changeControlState);
	connect(ui_control.radioButton_Bezier, &QRadioButton::clicked, this, &GeoModeling::changeCurveType);
	connect(ui_control.radioButton_Quadric_Bspline, &QRadioButton::clicked, this, &GeoModeling::changeCurveType);
	connect(ui_control.radioButton_Cubic_Bspline, &QRadioButton::clicked, this, &GeoModeling::changeCurveType);
	connect(ui_control.radioButton_Sampling, &QRadioButton::clicked, this, &GeoModeling::changeCurveType);
	connect(ui_control.radioButton_Subdivision, &QRadioButton::clicked, this, &GeoModeling::changeCurveType);
	connect(ui_control.radioButton_generator, &QRadioButton::clicked, this, &GeoModeling::changePlane);
	connect(ui_control.radioButton_trajactory, &QRadioButton::clicked, this, &GeoModeling::changePlane);
	connect(ui_control.checkBox_closed, &QCheckBox::stateChanged, this, &GeoModeling::changeCurveType);
	connect(ui_control.spinBox_precision, QOverload<int>::of(&QSpinBox::valueChanged), this, &GeoModeling::changePrecision);
	connect(ui_control.pushButton_clear, &QPushButton::clicked, this, &GeoModeling::clearState);
	connect(ui_control.pushButton_revolution, &QPushButton::clicked, this, &GeoModeling::doYRevolution);
	connect(ui_control.pushButton_sweep, &QPushButton::clicked, this, &GeoModeling::doSweep);
	connect(ui_control.spinBox_slices, QOverload<int>::of(&QSpinBox::valueChanged), this, &GeoModeling::changeSlices);
	connect(ui_control.pushButton_extrusion, &QPushButton::clicked, this, &GeoModeling::doZExtrusion);
	connect(ui_control.doubleSpinBox_zdepth, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &GeoModeling::changeZdepth);
	connect(ui_control.pushButton_saveMesh, &QPushButton::clicked, this, &GeoModeling::doSaveMesh);
	connect(ui_control.pushButton_loadMesh, &QPushButton::clicked, this, &GeoModeling::doLoadMesh);
	connect(ui_control.pushButton_DooSabin, &QPushButton::clicked, this, &GeoModeling::doSubDooSabin);
	connect(ui_control.pushButton_CatmullClark, &QPushButton::clicked, this, &GeoModeling::doCatmullClark);
	connect(ui_control.pushButton_Loop, &QPushButton::clicked, this, &GeoModeling::doLoop);
	connect(ui_control.pushButton_GenCtlPoly, &QPushButton::clicked, this, &GeoModeling::doGenCtlPoly);
	connect(ui_control.pushButton_BezierSurface, &QPushButton::clicked, this, &GeoModeling::doGenBezierSurface);
	connect(ui_control.pushButton_CubicBsplineSurface, &QPushButton::clicked, this, &GeoModeling::doGenCubicSplineSurface);
	connect(ui_control.pushButton_NNCrust, &QPushButton::clicked, this, &GeoModeling::doNNCrust);
	connect(ui_control.pushButton_Crust, &QPushButton::clicked, this, &GeoModeling::doCrust);
}
