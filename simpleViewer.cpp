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

#include "simpleViewer.h"

using namespace std;

void drawBigAxis(qreal length)
{
     glLineWidth(2.);
     glBegin(GL_LINES);
     glColor3d(1.,0.,0.);
     glVertex3d(-length, 0., 0);
     glVertex3d(length, 0., 0.);
     glEnd();
      glLineWidth(1.);
}

// Draws a spiral
void Viewer::draw() {

  glColor3d(1.,1.,1.);
  drawGrid(size,2*size);
  //drawBigAxis(20.);
  drawAxis(1.);


    glBegin (GL_TRIANGLES);
    for (unsigned int i = 0; i < mesh.T.size (); i++)
        for (unsigned int j = 0; j < 3; j++) {
            const MeshVertex & v = mesh.V[mesh.T[i].v[j]];
            glNormal3f (v.n[0], v.n[1], v.n[2]);
            glVertex3f (v.p[0], v.p[1], v.p[2]);
        }
    glEnd ();

}


void Viewer::init() {
    setSceneRadius(50.);
    mesh.loadOFF (std::string("C:\\Users\\pups\\libqgl\\libQGLViewer-2.8.0\\examples\\simpleViewer\\models\\arma.off"));
    size = mesh.size;
}

