#include "GeoModeling.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	a.setAttribute(Qt::AA_UseDesktopOpenGL);
	GeoModeling w;
	w.show();
	return a.exec();
}
