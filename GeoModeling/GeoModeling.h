#pragma once

#include <QtWidgets/QMainWindow>
#include <QDockWidget>
#include "ui_GeoModeling.h"
#include "ui_Control.h"

class GeoModeling : public QMainWindow
{
	Q_OBJECT

public:
	GeoModeling(QWidget *parent = Q_NULLPTR);

private:
	Ui::GeoModelingClass ui;
	Ui::ControlPanel ui_control;
};
