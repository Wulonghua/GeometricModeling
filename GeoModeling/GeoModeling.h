#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_GeoModeling.h"

class GeoModeling : public QMainWindow
{
	Q_OBJECT

public:
	GeoModeling(QWidget *parent = Q_NULLPTR);

private:
	Ui::GeoModelingClass ui;
};
