#include "GeoModeling.h"

GeoModeling::GeoModeling(QWidget *parent)
	: QMainWindow(parent)
{
	QWidget *dockw = new QWidget();
	ui_control.setupUi(dockw);
	QDockWidget *dock = new QDockWidget();
	dock->setWidget(dockw);
	dock->setMinimumWidth(240);
	this->addDockWidget(Qt::LeftDockWidgetArea, dock);
	ui.setupUi(this);

	m_curve = std::make_shared<Curve>();
	ui.glWidget->setCurve(m_curve);

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
	ui.glWidget->updateRender();
}

void GeoModeling::clearState()
{
	m_curve->reset();
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
}
