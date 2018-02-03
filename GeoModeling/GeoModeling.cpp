#include "GeoModeling.h"

GeoModeling::GeoModeling(QWidget *parent)
	: QMainWindow(parent)
{
	QWidget *dockw = new QWidget();
	ui_control.setupUi(dockw);
	QDockWidget *dock = new QDockWidget();
	dock->setWidget(dockw);
	dock->setMinimumWidth(300);
	this->addDockWidget(Qt::LeftDockWidgetArea, dock);
	ui.setupUi(this);

	m_curve = std::make_shared<Curve>();
}
