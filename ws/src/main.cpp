#include "MainWindow.hpp"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    
    // Set application properties
    app.setApplicationName("Fish Simulation");
    app.setApplicationVersion("1.0");
    app.setOrganizationName("BioMorf");
    
    MainWindow window;
    window.show();
    
    return app.exec();
}