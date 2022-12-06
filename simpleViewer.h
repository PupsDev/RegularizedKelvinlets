/****************************************************************************

 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of the QGLViewer library version 2.8.0.

 http://www.libqglviewer.com - contact@libqglviewer.com

 This file may be used under the terms of the GNU General Public License 
 versions 2.0 or 3.0 as published by the Free Software Foundation and
 appearing in the LICENSE file included in the packaging of this file.
 In addition, as a special exception, Gilles Debunne gives you certain 
 additional rights, described in the file GPL_EXCEPTION in this package.

 libQGLViewer uses dual licensing. Commercial/proprietary software must
 purchase a libQGLViewer Commercial License.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*****************************************************************************/

#include <QGLViewer/qglviewer.h>
#include <QMatrix4x4>
#include <Mesh.h>
#include <QGLViewer/camera.h>

#include <QOpenGLFramebufferObject>
#include <QOpenGLContext>
#include <QOpenGLExtraFunctions>
#include <QGLViewer/manipulatedFrame.h>
using namespace qglviewer;
typedef struct KelvinParameters
{
    Eigen::Vector3d x0;
    Eigen::Vector3d x;
    double a,b,c;
    double epsilon;
    Eigen::Vector3d force;
    Eigen::Vector3d axis;
    double dt;

}KelvinParameters;
typedef struct KelvinResult
{
    Eigen::Matrix3d force;
    double density;
}KelvinResult;

class Viewer : public QGLViewer {
	public:
        Viewer() : wireframe_(false), flatShading_(false), hideFlag(true), grabFlag(false){};
	protected:
		virtual void draw();
		void mousePressEvent(QMouseEvent* e);
		void mouseMoveEvent(QMouseEvent *e);
		void mouseReleaseEvent(QMouseEvent* e);
        void wheelEvent(QWheelEvent *event);
		virtual void keyPressEvent(QKeyEvent *e);
		virtual void init();
		void postSelection(const QPoint &point);
		void move(Eigen::Vector3d force);

         Eigen::Vector3d rungeKutta( Eigen::Vector3d x,int index,const KelvinParameters&  parameters,  float h);
	private :
		Mesh mesh;
		
		ManipulatedFrame	mf;

		float size;
		bool wireframe_, flatShading_;
		bool selection;
        bool mySelection;
        bool startSelection;
		bool grabFlag;
        bool hideFlag;

        double radiusBall = 1.;

		QPoint startClick;
		QPoint currentMovement;
        double distanceToCamera;

        GLdouble thisMatrix[16];
		
		float offsetX=0.;
		float offsetY=0.;
        Eigen::Vector3d test;

		// Runge-Kutta vectors
		// a vector of point per point
		std::vector<std::vector<Eigen::Vector3d>> rungePoints;

		double mini;
		double maxi;
		int indiceTomove;


		Vec startHandler;
		Vec currentHandler;
		Vec endHandler;

};
