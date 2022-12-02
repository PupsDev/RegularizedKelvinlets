#include "Mesh.h"
#include <iostream>
#include <fstream>
#include <limits>
#include <cfloat>
#include <tuple>

bool compVec3(Vec3 a, Vec3 b)
{
    return std::make_tuple(a[0], a[1], a[2]) < std::make_tuple(b[0], b[1], b[2]);
}
void Mesh::loadOFF (const std::string & filename) {
    std::ifstream in (filename.c_str ());
    if (!in)
    {
        std::cout<<"can't read file"<<std::endl;
        exit (EXIT_FAILURE);
    }
    std::string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    V.resize (sizeV);
    T.resize (sizeT);

    for (unsigned int i = 0; i < sizeV; i++) {
        in >> V[i].p;
        V[i].pInit = V[i].p;
    }
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
        in >> s;
        for (unsigned int j = 0; j < 3; j++)
            in >> T[i].v[j];
    }
    Vec3 mini(FLT_MAX,FLT_MAX,FLT_MAX);
    Vec3 maxi(-FLT_MAX,-FLT_MAX,-FLT_MAX);
    for(auto & v : V)
    {
        mini = std::min(mini,v.p, compVec3);
        maxi = std::max(maxi,v.p,compVec3);
    }
    Vec3 d = maxi-mini;
    size = 1.;//d.length();
    in.close ();
    centerAndScaleToUnit ();
    //catmullClark();
    recomputeNormals ();

    for(auto v : V)
    {
        vertex.push_back(Eigen::Vector3d(v.p[0],v.p[1],v.p[2]));
        //std::cout<<"force"<<v[0]<<" "<<v[1]<<" "<<v[2]<<std::endl;
    }
}

void Mesh::catmullClark()
{
    for (unsigned int i = 0; i < T.size (); i++) {
        Vec3 centroid = 1./3. * (V[T[i].v[0]].p+V[T[i].v[1]].p+V[T[i].v[2]].p);
    
    }

}
void Mesh::recomputeNormals () {
    for (unsigned int i = 0; i < V.size (); i++)
        V[i].n = Vec3 (0.0, 0.0, 0.0);
    for (unsigned int i = 0; i < T.size (); i++) {
        Vec3 e01 = V[T[i].v[1]].p -  V[T[i].v[0]].p;
        Vec3 e02 = V[T[i].v[2]].p -  V[T[i].v[0]].p;
        Vec3 n = Vec3::cross (e01, e02);
        n.normalize ();
        for (unsigned int j = 0; j < 3; j++)
            V[T[i].v[j]].n += n;
    }
    for (unsigned int i = 0; i < V.size (); i++)
        V[i].n.normalize ();
}

void Mesh::centerAndScaleToUnit () {
    Vec3 c(0,0,0);
    for  (unsigned int i = 0; i < V.size (); i++)
        c += V[i].p;
    c /= V.size ();
    float maxD = (V[0].p - c).length()/10.;
     size = 10.*maxD;
    /*for (unsigned int i = 0; i < V.size (); i++){
        float m = (V[i].p - c).length();
        if (m > maxD)
            maxD = m;
    }*/
    for  (unsigned int i = 0; i < V.size (); i++) {
        V[i].p = (V[i].p - c) / maxD;
        V[i].pInit = (V[i].pInit - c) / maxD;
    }
}
