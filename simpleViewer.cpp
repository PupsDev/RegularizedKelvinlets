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
#include <qmessagebox.h>
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
Vec computeMovement(Vec startHandler, QPoint start, QPoint end)
{
    double dx=0,dy=0;
    if( (end.x() - start.x()) > 0)
    {
        dx += 0.1;
    }
    else
    {
        dx -= 0.1;
    }
    if( (end.y() - start.y()) > 0)
    {
        dy += 0.1;
    }
    else
    {
        dy -= 0.1;
    }
    //Vec ray = startHandler - this->camera->position();
    Vec ray(0.,1.,0.);
    Vec x = cross	(ray, startHandler.unit());
    return startHandler + dx*x  - dy*ray;
}
void printMatrix(GLdouble m[16])
{
    std::cout<<"matrix"<<std::endl;
    for(int i = 0 ; i < 4 ; i++)
    {
        for(int j = 0 ; j < 4 ; j++)
        {
            std::cout<<m[4*i+j]<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"matrix"<<std::endl;
}
void printMatrix(GLfloat m[16])
{
    std::cout<<"matrix"<<std::endl;
    for(int i = 0 ; i < 4 ; i++)
    {
        for(int j = 0 ; j < 4 ; j++)
        {
            std::cout<<m[4*i+j]<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"matrix"<<std::endl;
}
void Viewer::mousePressEvent(QMouseEvent* e)
{
    if ((e->modifiers() == Qt::ControlModifier))
    {
        selection = true;
        std::cout<<"Begin1"<<std::endl;
        QPoint point(e->pos());
        beginSelection(point);
        bool found;
        startHandler = this->camera()->pointUnderPixel(point,found);

    }
    else
        QGLViewer::mousePressEvent(e);
}
void Viewer::mouseMoveEvent(QMouseEvent *e)
{
    if (selection)
    {
        /*std::cout<<"move1"<<std::endl;
        QPoint point(e->pos());
        bool found;
        startHandler = this->camera()->pointUnderPixel(point,found);
        */

    }
    else if(mySelection)
    {
        //std::cout<<"move"<<std::endl;
        QPoint point(e->pos());
        //std::cout<<point.x()<<" "<<point.y()<<std::endl;
        currentMovement = point;
        bool found;
        startClick = currentMovement;
        auto v = this->camera()->pointUnderPixel(point,found);
        auto vec = v-this->camera()->position();
        vec = vec.unit();

        endHandler = this->camera()->position() + distanceToCamera*vec;
        Vec delta = endHandler - currentHandler;

        move(Eigen::Vector3d(delta[0],delta[1],delta[2]));
        currentHandler = endHandler;
        auto vec2 = currentHandler - this->camera()->position();
        distanceToCamera = vec2.norm();
        


    }
    else
        QGLViewer::mouseMoveEvent(e);
}
void Viewer::mouseReleaseEvent(QMouseEvent* e)
{
    if (selection)
    {
         selection = false;
         std::cout<<"End1"<<std::endl;
         QPoint point(e->pos());
         std::cout<<point.x()<<" "<<point.y()<<std::endl;
         bool found;
         startHandler = this->camera()->pointUnderPixel(point,found);

    }
    else if(mySelection)
    {
         std::cout<<"End"<<std::endl;
         QPoint point(e->pos());
         //std::cout<<point.x()<<" "<<point.y()<<std::endl;
         //bool found;
         //auto v = this->camera()->pointUnderPixel(point,found);
         //std::cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<std::endl;

        //glFlush();


         mySelection = false;
    }
    else
        QGLViewer::mouseReleaseEvent(e);
}
void Viewer::wheelEvent(QWheelEvent *event)
{
    double angle = event->angleDelta().y();

    radiusBall += angle/120;
    std::cout<<"radiusBall "<<radiusBall<<std::endl;
}
void Viewer::postSelection(const QPoint &point) {
    bool found=false;
   	
    startHandler = this->camera()->pointUnderPixel(point,found);

    startClick = point;

    auto vec = startHandler - this->camera()->position();
    currentHandler = startHandler;
    distanceToCamera = vec.norm();

    mySelection = true;

    startSelection = !startSelection;
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
            break;
            case  Qt::Key_F:
            flatShading_ = !flatShading_;
            handled = true;
            break;
        case Qt::Key_Up:
            if(offsetY<0)
            offsetY*=-1;
                test[1]+=dy;
            handled = true;
            break;
        case Qt::Key_Down:
            handled = true;
            break;
        }
        if(handled)
        {
        //move();
        update();
        }
    }


    if (!handled)
        QGLViewer::keyPressEvent(e);
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
void Viewer::move(Eigen::Vector3d force) {
    double mu = 1.;
    double v = -1000000000;
    double a = 1./(4*M_PI*mu);
    double b = a / (4*(1.-v));
    double epsilon = radiusBall;
    //double radius = 1.;
    //double repsilon =r_epsilon(radius,epsilon);
    std::cout<<"force"<<force[0]<<" "<<force[1]<<" "<<force[2]<<std::endl;
    force = 50.*radiusBall*force;
    std::cout<<"force"<<force[0]<<" "<<force[1]<<" "<<force[2]<<std::endl;
    std::vector<Eigen::Vector3d> points;
    for(auto vertex : mesh.V)
    {
        points.push_back(Eigen::Vector3d(vertex.p[0],vertex.p[1],vertex.p[2]));
    }
    auto identity = Eigen::Matrix3d::Identity();
    Eigen::Vector3d x0 = Eigen::Vector3d(startHandler[0],startHandler[1],startHandler[2]);//points[indiceTomove];

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

        Eigen::Vector3d f = density_fonction(epsilon, repsilon) * force;
        Eigen::Vector3d uForce = Kef * f;

        points[i]+=uForce;
    }
    int k =0;
    for(auto p : points)
    {
        mesh.V[k++].p = Vec3(p[0],p[1],p[2]);
    }

}
void drawSphere(Vec origin, double r, int lats, int longs) {
    int i, j;
    for(i = 0; i <= lats; i++) {
        double lat0 = M_PI * (-0.5 + (double) (i - 1) / lats);
        double z0  = sin(lat0);
        double zr0 =  cos(lat0);

        double lat1 = M_PI * (-0.5 + (double) i / lats);
        double z1 = sin(lat1);
        double zr1 = cos(lat1);

        glBegin(GL_QUAD_STRIP);
        for(j = 0; j <= longs; j++) {
            double lng = 2 * M_PI * (double) (j - 1) / longs;
            double x = cos(lng);
            double y = sin(lng);

            glNormal3f(origin[0]+x * zr0,origin[1]+ y * zr0, origin[2]+z0);
            glVertex3f(origin[0]+x * zr0*r,origin[1]+ y * zr0*r, origin[2]+z0*r);
            glNormal3f(origin[0]+x * zr1,origin[1]+ y * zr1, origin[2]+z1);
            glVertex3f(origin[0]+x * zr1*r,origin[1]+ y * zr1*r, origin[2]+z1*r);
        }
        glEnd();
    }
}
void Viewer::draw() {

    glColor3d(1.,1.,1.);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    glEnable(GL_LIGHTING);

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
    /*if(!startSelection)
        glColor4f (0.0, 1.0, 0.0, 0.2);
    else
        glColor4f (1.0, 0.0, 0.0, 0.2);
        */
    //glVertex3f (startHandler[0], startHandler[1], startHandler[2]);

     glDisable(GL_LIGHTING);
    glPointSize(8.);
    glBegin (GL_POINTS);
    int k=0;

    //glDisable(GL_LIGHTING);
    for(auto v : mesh.V)
    {
        Vec vToSphere = startHandler-Vec(v.p[0],v.p[1],v.p[2]);
        if(vToSphere.norm()>radiusBall)
        {
            glColor3f(	0.,0.0,1.);
        }
        else
        {
            glColor3f(	0.,1.0,0.);
        }
        double offset =0.005;
        glVertex3f (v.p[0]+offset*v.n[0], v.p[1]+offset*v.n[1], v.p[2]+offset*v.n[2]);
        k++;
    }
    glEnd();
    glColor4f (0.0, 1.0, 0.0, 0.2);
    drawSphere(startHandler, radiusBall, 40, 40);
     glColor4f (1.0, 0.0, 0.0, 0.2);
    drawSphere(currentHandler, radiusBall, 40, 40);
    //glEnable(GL_LIGHTING);

    update();

}
void Viewer::init() {

    setSceneRadius(50.);
    setMouseTracking(true);
     mf = ManipulatedFrame();
     startSelection = true;

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
    GLfloat light_position[] = { 0.0, 0.0, 10.0, 0.0 };
    GLfloat light_position2[] = { 0.0, 0.0, 10.0, 0.0 };

    GLfloat light_diffuse[] = { 1.0, 0.0, 0.0, 0.0 };

    glLightfv(GL_LIGHT1, GL_POSITION, light_position);
    glLightfv(GL_LIGHT2, GL_POSITION, light_position2);

    glEnable(GL_LIGHTING);

    GLfloat cyan[] = {0.f, .8f, .8f, 1.f};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, cyan);
    //cyan[] = {0.f, .8f, .8f, 1.f};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, cyan);

    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, light_diffuse);

    this->camera()->setPosition(Vec(0.,0.,-20.));
     this->camera()->setViewDirection(Vec(0.,0.,1.));

    startHandler  = Vec(0.,0.,0.);
    currentHandler  = Vec(0.,0.,0.);
    //glEnable(GL_LIGHT1);
    //glEnable(GL_LIGHT2);

}

