#include "GeoModeling.h"

GeoModeling::GeoModeling(QWidget *parent)
	: QMainWindow(parent)
{
	QWidget *dockw = new QWidget();
	ui_control.setupUi(dockw);
	QDockWidget *dock = new QDockWidget();
	dock->setWidget(dockw);
	dock->setMinimumWidth(230);
	this->addDockWidget(Qt::LeftDockWidgetArea, dock);
	ui.setupUi(this);

	m_curve = std::make_shared<Curve>();
	m_mesh = std::make_shared<Mesh>();
	ui.glWidget->setCurve(m_curve);
	ui.glWidget->setMesh(m_mesh);

	initConnections();
}

void GeoModeling::changeControlState()
{
	if (ui_control.radioButton_addPoints->isChecked())
	{
		m_curve->m_contrlType = Curve::ADD;
	}
	else if (ui_control.radioButton_movePoints->isChecked())
	{
		m_curve->m_contrlType = Curve::MOVE;
	}
	else if (ui_control.radioButton_View->isChecked())
	{
		m_curve->m_contrlType = Curve::VIEW;
	}
}

void GeoModeling::changeCurveType()
{
	if (ui_control.radioButton_Sampling->isChecked())
	{
		m_curve->m_genType = Curve::SAMPLE;
		ui_control.label_precision->setText(QString("Precisions: "));
		m_curve->n_precision = ui_control.spinBox_precision->value();
	}
	else if (ui_control.radioButton_Subdivision->isChecked())
	{
		m_curve->m_genType = Curve::SUBDIVISION;
		ui_control.label_precision->setText(QString("Iterations: "));
		m_curve->n_precision = ui_control.spinBox_precision->value();
	}

	if (ui_control.radioButton_Bezier->isChecked())
	{
		m_curve->m_curveType = Curve::Bezier;
		m_curve->generateBezierPoints();
	}
	else if (ui_control.radioButton_Quadric_Bspline->isChecked())
	{
		m_curve->m_curveType = Curve::Quadric_B_spline;
		m_curve->generateQuadBspline();
	}
	else if (ui_control.radioButton_Cubic_Bspline->isChecked())
	{
		m_curve->m_curveType = Curve::Cubic_B_spline;
		m_curve->generateCubicBspline();
	}
	ui.glWidget->updateRender();
}

void GeoModeling::changePrecision(int x)
{
	m_curve->n_precision = x;
	if (ui_control.radioButton_Bezier->isChecked())
	{
		m_curve->m_curveType = Curve::Bezier;
		m_curve->generateBezierPoints();
	}
	else if (ui_control.radioButton_Quadric_Bspline->isChecked())
	{
		m_curve->m_curveType = Curve::Quadric_B_spline;
		m_curve->generateQuadBspline();
	}
	else if (ui_control.radioButton_Cubic_Bspline->isChecked())
	{
		m_curve->m_curveType = Curve::Cubic_B_spline;
		m_curve->generateCubicBspline();
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

void GeoModeling::initConnections()
{
	connect(ui_control.radioButton_addPoints, &QRadioButton::clicked, this, &GeoModeling::changeControlState);
	connect(ui_control.radioButton_movePoints, &QRadioButton::clicked, this, &GeoModeling::changeControlState);
	connect(ui_control.radioButton_View, &QRadioButton::clicked, this, &GeoModeling::changeControlState);
	connect(ui_control.radioButton_Bezier, &QRadioButton::clicked, this, &GeoModeling::changeCurveType);
	connect(ui_control.radioButton_Quadric_Bspline, &QRadioButton::clicked, this, &GeoModeling::changeCurveType);
	connect(ui_control.radioButton_Cubic_Bspline, &QRadioButton::clicked, this, &GeoModeling::changeCurveType);
	connect(ui_control.radioButton_Sampling, &QRadioButton::clicked, this, &GeoModeling::changeCurveType);
	connect(ui_control.radioButton_Subdivision, &QRadioButton::clicked, this, &GeoModeling::changeCurveType);
	connect(ui_control.spinBox_precision, QOverload<int>::of(&QSpinBox::valueChanged), this, &GeoModeling::changePrecision);
	connect(ui_control.pushButton_clear, &QPushButton::clicked, this, &GeoModeling::clearState);
	connect(ui_control.pushButton_revolution, &QPushButton::clicked, this, &GeoModeling::doYRevolution);
	connect(ui_control.spinBox_slices, QOverload<int>::of(&QSpinBox::valueChanged), this, &GeoModeling::changeSlices);
	connect(ui_control.pushButton_extrusion, &QPushButton::clicked, this, &GeoModeling::doZExtrusion);
	connect(ui_control.doubleSpinBox_zdepth, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &GeoModeling::changeZdepth);
}
