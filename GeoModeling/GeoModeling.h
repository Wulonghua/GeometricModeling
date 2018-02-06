#pragma once

#include <QtWidgets/QMainWindow>
#include <QDockWidget>
#include <memory>

#include "Curve.h"
#include "ui_GeoModeling.h"
#include "ui_Control.h"

class GeoModeling : public QMainWindow
{
	Q_OBJECT

public:
	GeoModeling(QWidget *parent = Q_NULLPTR);
	
public slots:
	void changeControlState();
	void changeCurveType();


private:
	void initConnections();
	Ui::GeoModelingClass ui;
	Ui::ControlPanel ui_control;

	std::shared_ptr<Curve> m_curve;
};
