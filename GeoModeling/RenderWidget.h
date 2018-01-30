#pragma once

#include "QGLViewer\qglviewer.h"
#include <QKeyEvent>

class RenderWidget :
	public QGLViewer
{
public:
	RenderWidget(QWidget *parent);
	~RenderWidget();

protected:
	virtual void draw();
	virtual void init();
	virtual void postDraw();
	virtual void keyPressEvent(QKeyEvent *e);
};

