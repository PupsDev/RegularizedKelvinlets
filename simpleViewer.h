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
#include "Eigen/Dense"
#include <QOpenGLFramebufferObject>
#include <QOpenGLContext>
#include <QOpenGLExtraFunctions>
using namespace qglviewer;

class Viewer : public QGLViewer {
public:
  Viewer() : wireframe_(false), flatShading_(false), grabFlag(false){};
protected:
  virtual void draw();
  virtual void keyPressEvent(QKeyEvent *e);
   //void mouseMoveEvent(QMouseEvent* const e);
  virtual void init();
  void drawSelection(double radius);
  void Grab(double radius);
  void move();
private :
QOpenGLFramebufferObject* mFBO = nullptr;
      Mesh mesh;
      Camera* mycamera;
      float size;
      bool wireframe_, flatShading_;
      bool selection;
      bool grabFlag;

      float offsetX=0.;
      float offsetY=0.;
     Eigen::Vector3d test;
    double mini;
    double maxi;
    int indiceTomove;
      QMatrix4x4 modelViewMatrix;
      QMatrix4x4 projectionMatrix;
       QPoint pixelMouse;

        QPoint pos, prevPos;
        bool moved;

};
