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
#include <QCursor>
#include <QKeyEvent>
#include <QMap>
#include <QMenu>
#include <QPoint>
#include <QMouseEvent>
#include <QCursor>
#include <limits>
#include <cfloat>
#include <QtMath>
#include <QGLViewer/camera.h>
#ifdef Q_OS_WIN
#include "gl/GLU.h"
#else
#include <glu.h>
#endif
using namespace std;

const float dy = 3.;

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
void Viewer::keyPressEvent(QKeyEvent *e) {
  // Get event modifiers key
  const Qt::KeyboardModifiers modifiers = e->modifiers();

  bool handled = false;
  if((modifiers == Qt::NoButton))
  {
      switch(e->key())
      {

          case  Qt::Key_W:
            wireframe_ = !wireframe_;
            handled = true;

          break;
          case  Qt::Key_B:
            selection = !selection;
            handled = true;

          break;
          case  Qt::Key_G:
            grabFlag = !grabFlag;
            handled = true;
                //Grab(0.5);

          break;
          case  Qt::Key_F:
            flatShading_ = !flatShading_;
            handled = true;
          break;
      case Qt::Key_Up:
          //std::cout<<"ok up"<<std::endl;
          if(offsetY<0)
            offsetY*=-1;
          //if(offsetY>0)offsetY+=dy;
              test[1]+=dy;
          handled = true;
          break;
        case Qt::Key_Down:
          //std::cout<<"ok down"<<std::endl;
          /*
          if(offsetY>0)
            offsetY*=-1;
           //if(offsetY<0)offsetY-=dy;
           */
          //this->camera()->setPosition (Vec(0.,0.,0.));
          handled = true;
          break;
      }
      if(handled)
      {
        move();
        update();
      }
  }


  if (!handled)
    QGLViewer::keyPressEvent(e);
}
/*
void Viewer::mouseMoveEvent(QMouseEvent *event)
{
    glEnable(GL_DEPTH_TEST);
    if(!selection)
    {

        QPointF wpos = event->windowPos();
      

        QColor color;
        glReadPixels(wpos.x(), this->height() - wpos.y(), 1, 1, GL_RGBA, GL_FLOAT, &color);

        QVector3D screenCoordinates = QVector3D(wpos.x(),  this->height() - wpos.y(), 0.);
        modelViewMatrix.setToIdentity();
        projectionMatrix.setToIdentity();

        //GLfloat mvmatrix[16];
        //GLfloat pmatrix[16];
        //this->camera()->getModelViewMatrix (mvmatrix);
        //this->camera()->getProjectionMatrix (pmatrix);

        QVector3D mouseIn3D = screenCoordinates.unproject(modelViewMatrix,
                                                          projectionMatrix,
                                                          QRect(0, 0, width(), height()));
        //std::cout << "Clicked on pixel " << event->x() << ", " << event->y() <<endl;
        bool found;
        auto point = QPoint(event->x(),event->y());
        auto cam = this->camera();
          pixelMouse = QPoint(wpos.x(),wpos.y());
        Vec t = cam->pointUnderPixel(QPoint(wpos.x(),wpos.y()),found);
        std::cout <<t[0]<<" "<<t[1]<<" "<<t[2]<<std::endl;

        auto pp = this->camera()->unprojectedCoordinatesOf(Vec(wpos.x(),wpos.y(),0.));
        //std::cout<<found<<std::endl;
        //std::cout <<t[0]<<" "<<t[1]<<" "<<t[2]<<std::endl;
        //std::cout <<pp[0]<<" "<<pp[1]<<" "<<pp[2]<<std::endl;
        /*
        gluUnProject(	GLdouble winX,
            GLdouble winY,
            GLdouble winZ,
            const GLdouble * model,
            const GLdouble * proj,
            const GLint * view,
            GLdouble* objX,
            GLdouble* objY,
            GLdouble* objZ);
        */
/*
        GLdouble Bx,By,Bz;
        GLint viewport[4];
        GLdouble mvmatrix[16], projmatrix[16];
        //glGetDoublev(GL_MODELVIEW_MATRIX,mvmatrix );
        //glGetDoublev(GL_PROJECTION_MATRIX,projmatrix );
        this->camera()->getModelViewMatrix (mvmatrix);
        this->camera()->getProjectionMatrix (projmatrix);
        glGetIntegerv(GL_VIEWPORT,viewport );
        GLfloat*  depths;
        std::cout <<"test "<<std::endl;
        glReadPixels (wpos.x(),this->height()-wpos.y(), this->width(), this->height(), GL_DEPTH_COMPONENT, GL_FLOAT, depths);
        glGetError();
        std::cout <<"depths "<<depths<<std::endl;
        //gluUnProject(wpos.x(),400-wpos.y(),depths, mvmatrix, projmatrix,viewport,&Bx,&By,&Bz  );
        // std::cout <<Bx<<" "<<By<<" "<<Bz<<std::endl;


    }
}

*/
/*
void Viewer::mouseMoveEvent(QMouseEvent* const e)
 {

   if (moved)
   {
       // Add position delta to current pos
       pos += e->pos() - prevPos;
       prevPos = e->pos();
   }
   if(selection)
   {
     // std::cout << "Clicked on pixel " << e->x() << ", " << e->y() <<endl;
      pixelMouse = QPoint(e->x(),e->y());
      update();
   }

 }
*/
// Draws a spiral
void Viewer::drawSelection(double radius) {
    glColor3d(1.,0.,0.);
    glLineWidth(4.);
    glBegin (GL_LINES);
    double theta = 0.;
    int segments = 64;
    QPoint p = mapFromGlobal(QCursor::pos());
    //p =  pixelMouse;
    //std::cout << p.x()<<" "<<p.y()<<std::endl;
    std::cout << this->width()<<" "<<this->height()<<std::endl;

    for(int i = 0 ; i < segments;i++)
    {
        theta = ((2*3.1415)*i)/segments;
        double x = radius* cos(theta);
        double y = radius* sin(theta);

        //double x2 = (double)(p.x() -this->width()/2. )/250. +x;
        //double y2 = (double)((this->height()-p.y()) -this->height()/2. )/250. +y;

        double xscreen = (double)(2.*p.x()/ (double)this->width());
        double yscreen = (double)(2.*(this->height()-p.y()) /(double)this->height());
        std::cout << xscreen<<" "<<yscreen<<std::endl;

        /*glVertex3fv(camera()->worldCoordinatesOf(
            qglviewer::Vec( x2,   y2, -1.0)));
        */
        glVertex3fv(qglviewer::Vec(x+xscreen-1.,1.25*(y+yscreen-1.),0.));
    }

    glEnd ();
    update();
}
void Viewer::Grab(double radius) {

    QPoint mousePos = QCursor::pos();
    bool found=false;
    auto p = mycamera->pointUnderPixel(mousePos,found);
    std::cout<<found<<std::endl;
    std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;	
}
double density_fonction(double epsilon, double repsilon)
{
  
  return (15*pow(epsilon,4)/(8*M_PI))/(pow(repsilon,7));
}
double r_epsilon(double radius,double epsilon)
{
  return  sqrt(radius*radius + epsilon*epsilon);
}
Eigen::Vector3d force_u()
{

    //return test;
}
void Viewer::move() {
    double mu = 1.;
    double v = 0.5;
    double a = 1./(4*M_PI*mu);
    double b = a / (4*(1.-v));
    double epsilon = 2.;
    //double radius = 1.;
    //double repsilon =r_epsilon(radius,epsilon);
    std::vector<Eigen::Vector3d> points;
    for(auto vertex : mesh.V)
    {
        points.push_back(Eigen::Vector3d(vertex.p[0],vertex.p[1],vertex.p[2]));
    }
    auto identity = Eigen::Matrix3d::Identity();
    Eigen::Vector3d x0 = points[indiceTomove];
    
    for(int i = 0 ; i < points.size(); i++ )
    {
        Eigen::Vector3d rvector = x0 - points[i];
        double radius = rvector.norm();
        double repsilon =r_epsilon(radius,epsilon);

        auto first = ((a-b)*1./repsilon)*identity;
        auto second = b/(repsilon*repsilon*repsilon);
        auto last = ((a*epsilon*epsilon)/(2*repsilon*repsilon*repsilon))*identity;
        
        Eigen::Matrix3d rrt = rvector * rvector.transpose();
        
        auto Kef = first + second*rrt + last;
 
        Eigen::Vector3d f = density_fonction(epsilon, repsilon) * test;
        Eigen::Vector3d uForce = Kef * f;

        points[i]+=uForce;
    }
    int k =0;
    for(auto p : points)
    {
       mesh.V[k++].p = Vec3(p[0],p[1],p[2]);
    }

}
void Viewer::draw() {

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (float)600 / 400, 0.1, 100.);
    glMatrixMode(GL_MODELVIEW);
    glColor3d(1.,1.,1.);
    //glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDepthMask(GL_TRUE);

    std::cout<<this->camera()->position()[0]<<" "<<this->camera()->position()[1]<<" "<<this->camera()->position()[2]<<std::endl;



    drawGrid(size,2*size);
    //drawBigAxis(20.);
    drawAxis(1.);
    wireframe_ ?  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) :  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    flatShading_ ? glShadeModel(GL_FLAT) :  glShadeModel(GL_SMOOTH);
    glBegin (GL_TRIANGLES);
    for (unsigned int i = 0; i < mesh.T.size (); i++)
        for (unsigned int j = 0; j < 3; j++) {
            const MeshVertex & v = mesh.V[mesh.T[i].v[j]];
            glNormal3f (v.n[0], v.n[1], v.n[2]);
            glVertex3f (v.p[0], v.p[1], v.p[2]);
        }
    glEnd ();
    glPointSize(8.);
    glBegin (GL_POINTS);
    int k=0;
    for(auto v : mesh.V)
    {
        if(k!=indiceTomove)
        {
            glColor3f(	0.,0.0,1.);
        }
        else
        {
            glColor3f(	0.,1.0,0.);
        }
        glVertex3f (v.p[0]+0.005*v.p[0], v.p[1]+0.005*v.p[1], v.p[2]+0.005*v.p[2]);
        k++;
    }
    glEnd();
    if(selection)
    {
        glDisable(GL_LIGHTING);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        glPushMatrix();

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glPushMatrix();

        drawSelection(0.1);
        glMatrixMode(GL_MODELVIEW);
    }
    update();

}
void Viewer::init() {

   
    modelViewMatrix.lookAt(
        QVector3D(0.0, 0.0, 10.0), // Eye
        QVector3D(0.0, 0.0, 0.0), // Focal Point
        QVector3D(0.0, 1.0, 0.0)); // Up vector

    // Window size is fixed at 800.0 by 600.0
    projectionMatrix.perspective(45.0, 800.0 / 600.0, 1.0, 100.0);

    QMatrix4x4 mvp = (projectionMatrix*modelViewMatrix );




    setSceneRadius(50.);

    setMouseTracking(true);

    auto c =this->camera();
    auto p = c->position();
    std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;




    this->camera()->setPosition(Vec(0.,0.,-20.));
    this->camera()->setPivotPoint(Vec(0.,0.,0.));
    this->camera()->lookAt(Vec(0.,0.,0));

    p = this->camera()->position();
    std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;

    

    mesh.loadOFF (std::string("C:\\Users\\pups\\libqgl\\libQGLViewer-2.8.0\\examples\\simpleViewer\\models\\arma.off"));
    size = mesh.size;
    mini = FLT_MAX;
    maxi = -FLT_MAX;

    for(auto vertex : mesh.V)
    { 
    mini = std::min(mini,vertex.p[1]);    
    maxi = std::max(maxi,vertex.p[1]);

    }
    offsetY = dy;
    int k =0;
    for(auto & vertex : mesh.V)
    { 
        if( vertex.p[1] == maxi )
        {
            indiceTomove = k; 
        }
        k++;

    }
    std::cout<<"pts->"<<mesh.V[indiceTomove].p[0]<<" "<< mesh.V[indiceTomove].p[1]<<" "<<mesh.V[indiceTomove].p[2]<<std::endl;
    // Add custom key description (see keyPressEvent).
    setKeyDescription(Qt::Key_W, "Toggles wire frame display");
    setKeyDescription(Qt::Key_F, "Toggles flat shading display");
        test = Eigen::Vector3d(0.,15.,0.);
    // glViewport(0,0,1600,600);
}

