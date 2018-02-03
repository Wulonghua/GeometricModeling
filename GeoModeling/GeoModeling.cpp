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

void GeoModeling::initConnections()
{
	connect(ui_control.radioButton_addPoints, &QRadioButton::clicked, this, &GeoModeling::changeControlState);
	connect(ui_control.radioButton_movePoints, &QRadioButton::clicked, this, &GeoModeling::changeControlState);
	connect(ui_control.radioButton_View, &QRadioButton::clicked, this, &GeoModeling::changeControlState);
}
