#include <string>
#include <fstream>
#ifdef _WIN32
#include <windows.h>
#endif
#include <QTcpServer>
#include <QTcpSocket>
#include <QtGui>
#include <QTcpServer>
#include <QMessageBox>

#include "misc.h"
#include "log.h"

#include "libmcell.h"

#include "global.h"

LOG_USE();
char msg[2048];

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
SocketHandler::SocketHandler(int newport, QObject *parent)
	: QThread(parent)
{
    exiting = false;
    port = newport;
	sprintf(msg,"SocketHandler: port: %d",port);
	LOG_MSG(msg);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
SocketHandler::~SocketHandler() // make sure the worker object is destroyed
{
    exiting = true;
    wait();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::stop()
{
	LOG_MSG("SocketHandler::stop: set stopped");
	stopped = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::run()
{
//	QObject::moveToThread(this);
	sprintf(msg,"run: port: %d", port);
	LOG_MSG(msg);
	quint16 qport = port;
	QString addressStr = "127.0.0.1";
	QHostAddress hostAddress;
	hostAddress.setAddress(addressStr);
    tcpServer = new QTcpServer(this);
	stopped = false;
	connect(tcpServer, SIGNAL(newConnection()), this, SLOT(processor()), Qt::DirectConnection);
    if (!tcpServer->listen(hostAddress,qport)) {
 //       QMessageBox::critical(this, tr("Fortune Server"),
 //                              tr("Unable to start the server: %1.")
 //                              .arg(tcpServer->errorString()));
		sprintf(msg,"Unable to start the server: port: %d", port);
		LOG_MSG(msg);
        return;
    }
	sprintf(msg,"Listening on port: %d",tcpServer->serverPort());
	LOG_MSG(msg);
	LOG_MSG("serverAddress:");
	LOG_QMSG((tcpServer->serverAddress()).toString());
	bool timedOut = false;
	tcpServer->waitForNewConnection(-1,&timedOut);
	exec();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void SocketHandler::processor()
{
//	LOG_MSG("In processor");
    socket = tcpServer->nextPendingConnection();
	sprintf(msg,"got server connection: %p",socket);	
	LOG_MSG(msg);
    emit sh_connected();
	QString qdata;
	QByteArray ba;
	ba.resize(1024);
	while (true) {
		if (stopped) {
			LOG_MSG("SocketHandler::processor: stopped!");
			break;
		}
		socket->waitForReadyRead(100);
		int nb = socket->bytesAvailable();
		if (nb > 0) {
			ba = socket->readLine(1024);
			qdata = QString(ba);
			QStringList s = qdata.split("^",QString::SkipEmptyParts);
			for (int k=0; k<s.length(); k++) {
				emit sh_output(s[k]); // Emit signal to update GUI
				if (port == CPORT0) {
					LOG_QMSG(s[k]);
				}
			}
			if (quitMessage(qdata)) {
				sprintf(msg,"Closing connection: port: %d", port);
				LOG_MSG(msg);
		        break;
			} else {
//				LOG_MSG("No bytes yet");
			}
		}
	}
	socket->close();
	tcpServer->close();
	if (port == CPORT0) {
		emit sh_disconnected();		// Is it right that both threads emit this?
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
ExecThread::ExecThread(QString infile)
{
	inputFile = infile;
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::run()
{
	LOG_MSG("Invoking DLL...");
	int res=0;
    int nsumm_interval, hour;
	const char *infile, *outfile;
	QString infile_path, outfile_path;
	int len_infile, len_outfile;
    bool cused[16];

	infile_path = inputFile;
	QString casename = QFileInfo(inputFile).baseName();
	len_infile = infile_path.length();
	std::string std_infile = infile_path.toStdString();
	infile = std_infile.c_str();
	outfile_path = casename.append(".res");
	len_outfile = outfile_path.length();
	std::string std_outfile = outfile_path.toStdString();
    outfile = std_outfile.c_str();

	paused = false;
    execute(&ncpu,const_cast<char *>(infile),&len_infile,const_cast<char *>(outfile),&len_outfile,&res);
    if (res != 0) {
        terminate_run(&res);
        return;
    }
    get_dimensions(&nsteps, &Global::Ncirc, &Global::Nlong, &Global::DELTA_T);
    if (Global::hex_list) free(Global::hex_list);
    Global::hex_list = (HEXAHEDRON *)malloc(Global::Ncirc*Global::Nlong*sizeof(HEXAHEDRON));
    emit setupC();
    nsumm_interval = (60)/Global::DELTA_T;   // number of time steps per hour
    sprintf(msg,"exthread: nsteps: %d Ncirc: %d Nlong: %d DELTA_T: %f nsumm_interval: %d nt_vtk: %d",
            nsteps,Global::Ncirc,Global::Nlong,Global::DELTA_T,nsumm_interval,Global::nt_vtk);
    LOG_MSG(msg);

    Global::conc_nc = 0;
    hour = 0;

//    mutex1.lock();
//    get_summary(Global::summaryData);
//    mutex1.unlock();
//    emit summary(hour);		// Emit signal to initialise summary plots
//    summary_done.wait(&mutex3);

    for (int i=1; i <= nsteps+1; i++) {
		bool updated = false;
		if (paused && !updated) {
			snapshot();
            sprintf(msg,"got snapshot: i: %d",i);
            LOG_MSG(msg);
            updated = true;
		}
        while(paused || Global::leftb) {
            sleep(100);
		}
        if (stopped) {
            res = -1;
            break;
        }

        simulate_step(&res);
        if (res != 0) {
            LOG_MSG("simulate_step: error: res != 0");
            break;
        } else {
//            LOG_MSG("did simulate_step OK");
        }

        if (i%nsumm_interval == 0) {
//          if (i%nsumm_interval == -1) {
//            mutex1.lock();
//            get_summary(Global::summaryData);

//            conc_nc = 0;
//            vol_nv = 0;
//            oxy_nv = 0;
//            get_concdata(&conc_nc, &conc_dx, concData);
//            get_volprob(&vol_nv, &vol_v0, &vol_dv, volProb);
//            get_oxyprob(&oxy_nv, &oxy_dv, oxyProb);
//            mutex1.unlock();
            hour++;
            emit summary(hour);		// Emit signal to update summary plots, at hourly intervals
//            summary_done.wait(&mutex3);
        }

        if (stopped) {
            res = -1;
            break;
        }
        if (i%Global::nt_vtk == 0) {
            if (Global::showingVTK) {
				snapshot();
                Global::istep = i;
                sleep(10);
			}
		}
        if (stopped) {
            res = -1;
            break;
        }
    }
    LOG_MSG("ExecThread::run: stopped or completed");
    snapshot();
    LOG_MSG("got snapshot:");
    sleep(100);
	LOG_MSG("ExecThread::run: call terminate_run");
	terminate_run(&res);

	return;
}

/*
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::wait_to_go()
{
    for (;;) {
        sprintf(msg,"waiting");
        if (goflag || stopped) break;
    }
}
*/

//-----------------------------------------------------------------------------------------
// Note that storage for BC_list, DC_list, bond_list is provided in the GUI code
// (see mainwindow.cpp and transfer.h).  The DLL fills in the data, and the number of elements
// in the lists is returned in nBC_list, nDC_list, nbond_list.
//-----------------------------------------------------------------------------------------
void ExecThread::snapshot()
{
    get_scene(&Global::Nhex, Global::hex_list);
//    emit displayF(); // Emit signal to update Field display
//    LOG_MSG("did get_scene, emit display");
    emit display(); // Emit signal to update VTK display
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::saveGradient2D(int i)
{
    LOG_QMSG("saveGradient2D");
    paused = true;
    SimpleView2D *sv2D = new SimpleView2D();
    sv2D->makeFrame(i);
    paused = false;
    delete sv2D;
}
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::stop()
{
	stopped = true;
	LOG_MSG("ExecThread::stop: stopped");
}

	//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::pause()
{
	paused = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void ExecThread::unpause()
{
	paused = false;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool quitMessage(QString msg)
{
	if (msg.contains("__EXIT__",Qt::CaseSensitive))
		return true;
	else
		return false;
}
