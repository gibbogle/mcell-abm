// myvtk.h
#ifndef MYVTK_H
#define MYVTK_H

#include <QtGui>
#include <QtCore>
#include <QIODevice>
#include <QVTKWidget.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include "vtkSphereSource.h"
#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
//#include <vtkMPEG2Writer.h>
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkBMPWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSmartPointer.h>
#include <vtkImageCast.h>

#include <vtkAppendPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>

#include <vtkSuperquadricSource.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>

//#include <vtkConfigure.h>

#include <QInputDialog>
#include <QFileDialog>
#ifdef _WIN32
#include "windows.h"
#define sleep(n) Sleep(n)   // milliseconds
#else
#define sleep(n) usleep(1000*n)
#endif

using namespace std;

struct cell_pos {
	int tag;
//	int x, y, z;
    double x, y, z;
    double diameter;
//	double state;
	int state;
    int highlight;
};
typedef cell_pos CELL_POS;

struct bond_pos {
	int BCtag;
	int DCtag;
};
typedef bond_pos BOND_POS;

struct actor_str {
    bool active;
    vtkActor *actor;
    vtkSmartPointer<vtkMatrix4x4> vtkM;
    vtkSmartPointer<vtkTransform> transform;
};
typedef actor_str ACTOR_TYPE;

struct colour_str {
    double r, g, b;
};
typedef colour_str COLOUR_TYPE;

#define USE_CELLTYPE_COLOUR false
#define DISPLAY_SQUADS 0
#define DISPLAY_BLOCKS 1
#define DISPLAY_SPHERES 2

class MyVTK
{
private slots:
//    void on_checkBox_CELLDISPLAY_1_toggled(bool display);
//    void on_checkBox_CELLDISPLAY_2_toggled(bool);

public:
    MyVTK(QWidget *, QWidget *);
	~MyVTK();

    void key_canvas(QWidget *);
    void createMappers();
	void read_cell_positions(QString, QString, bool);
    void get_cell_positions();
	void init();
	void cleanup();
	void unpack(int x, double *, double *, double *);
    void renderCells();
    void setPoints(vtkSmartPointer<vtkPoints> p);
    void process_Mcells();
    void process_Tcells();
    void process_squads();
	bool startPlayer(QString, QTimer *, bool);
	bool nextFrame();
	void pause();
	void playon();
	void saveSnapshot(QString, QString);
    void startRecorder(QString basename, int nframes);
    void stopRecorder();
    void recorder();
    void stop();
    void set_celltype_colour(COLOUR_TYPE *, QString str);
    void makeColors(vtkSmartPointer<vtkUnsignedCharArray> colors, int nV);
    void set_diameter(double);
    void toggle_display_mode(int button);

    QList<CELL_POS > TCpos_list;
//	QList<CELL_POS > DCpos_list;
//	QList<BOND_POS > bondpos_list;
//	QList<vtkActor *> B_Actor_list;
    QList<ACTOR_TYPE> T_Actor_list;
    QList<ACTOR_TYPE> S_Actor_list;
//  QList<vtkActor *> D_Actor_list;
//  QList<ACTOR_TYPE> D_Actor_list;
//    QList<vtkActor *> Bnd_Actor_list;

    QWidget *page_VTK;
	QVTKWidget* qvtkWidget;
	vtkRenderWindow *renWin;	
	vtkRenderer* ren;
	vtkRenderWindowInteractor * iren;
    vtkPolyDataMapper *TcellMapper;
    vtkSmartPointer<vtkPolyDataMapper> squadMapper;
//	vtkPolyDataMapper *DcellMapper;
//	vtkPolyDataMapper *bondMapper;
//	vtkPolyDataMapper *FDcellMapper;
//	vtkMPEG2Writer *mpg;
//	vtkSmartPointer<vtkPNGWriter> writer;
//	vtkSmartPointer<vtkBMPWriter> writer;
    vtkSmartPointer<vtkJPEGWriter> jpgwriter;
    vtkSmartPointer<vtkPNGWriter> pngwriter;
//	vtkSmartPointer<vtkTIFFWriter> writer;
//	vtkSmartPointer<vtkImageCast> castFilter;
//	vtkWindowToImageFilter *w2img;
//	vtkSmartPointer<vtkPNGWriter> pngwriter;
//	vtkSmartPointer<vtkJPEGWriter> jpgwriter;
    vtkSmartPointer<vtkWindowToImageFilter> w2i;

// For Mcells
    vtkSmartPointer<vtkPoints> points;
//    vtkSmartPointer<vtkPolygon> polygon1, polygon2;
    vtkSmartPointer<vtkCellArray> polygons;
    vtkSmartPointer<vtkPolyData> polygonPolyData;
    vtkSmartPointer<vtkPolyDataMapper> hexmapper;
    vtkSmartPointer<vtkActor> hexactor;

	char msg[2048];
	double zoomlevel;
    double opacity;
	double Pi;
	bool DCmotion;
	bool DCfade;
	bool first_VTK;
	bool playing;
	bool paused;
	bool save_image;
    bool display_celltype[10];
    bool reset;
    bool display_spheres;
    double diameter;
    int display_mode;
//    QString celltype_colour[10];
    QColor celltype_colour[10];
    QString casename;
	int framenum;
	QTimer *timer;
	QString infile;
	QFile *playerData;
	QTextStream *playerStream;

    // Image recorder
//    bool record;
//    QString record_basename;
//    int record_nframes;
//    int record_it;
//    QTemporaryFile * tempFile;

public slots:
//    void on_checkBox_CELLDISPLAY_1_toggled(bool display);
//    void on_checkBox_CELLDISPLAY_2_toggled(bool);

};

#endif
