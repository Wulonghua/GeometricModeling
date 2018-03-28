#pragma once

#include <QtWidgets/QMainWindow>
#include <QDockWidget>
#include <QMessagebox>
#include <QFileDialog>
#include <memory>

#include "Curve.h"
#include "ui_GeoModeling.h"
#include "ui_Control.h"
#include "mesh.h"

class GeoModeling : public QMainWindow
{
	Q_OBJECT

public:
	GeoModeling(QWidget *parent = Q_NULLPTR);
	
public slots:
	void changeControlState();
	void changeCurveType();
	void changePlane();
	void changePrecision(int x);
	void changeSlices(int x);
	void changeZdepth(double x);
	void clearState();
	void doYRevolution();
	void doZExtrusion();
	void doSweep();
	void doSaveMesh();
	void doLoadMesh();
	void doSubDooSabin();

private:
	void initConnections();
	Ui::GeoModelingClass ui;
	Ui::ControlPanel ui_control;

	std::shared_ptr<Curve> m_curve;
	std::shared_ptr<Curve> m_traj;  // trajactory curve in yz plane
	std::shared_ptr<Mesh> m_mesh;

};
