#ifndef MESH_H
#define MESH_H


#include <vector>
#include <string>
#include "Vec3.h"



// -------------------------------------------
// Basic Mesh class
// -------------------------------------------

struct MeshVertex {
    inline MeshVertex () {}
    inline MeshVertex (const Vec3 & _p, const Vec3 & _n) : p (_p), pInit (_p), n (_n) {}
    inline MeshVertex (const MeshVertex & v) : p (v.p), pInit (v.pInit), n (v.n) {}
    inline virtual ~MeshVertex () {}
    inline MeshVertex & operator = (const MeshVertex & v) {
        p = v.p;
        pInit = v.pInit;
        n = v.n;
        return (*this);
    }
    double & operator [] (unsigned int c) {
        return p[c];
    }
    double operator [] (unsigned int c) const {
        return p[c];
    }

    // membres :
    Vec3 p; // une position
    Vec3 pInit; // une position
    Vec3 n; // une normale
};

struct MeshTriangle {
    inline MeshTriangle () {
        v[0] = v[1] = v[2] = 0;
    }
    inline MeshTriangle (const MeshTriangle & t) {
        v[0] = t.v[0];   v[1] = t.v[1];   v[2] = t.v[2];
    }
    inline MeshTriangle (unsigned int v0, unsigned int v1, unsigned int v2) {
        v[0] = v0;   v[1] = v1;   v[2] = v2;
    }
    inline virtual ~MeshTriangle () {}
    inline MeshTriangle & operator = (const MeshTriangle & t) {
        v[0] = t.v[0];   v[1] = t.v[1];   v[2] = t.v[2];
        return (*this);
    }
    unsigned int & operator [] (unsigned int c) {
        return v[c];
    }
    unsigned int operator [] (unsigned int c) const {
        return v[c];
    }

    // membres :
    unsigned int v[3];
};




class Mesh {
public:
    std::vector<MeshVertex> V;
    std::vector<MeshTriangle> T;
    float size;
    void loadOFF (const std::string & filename);
    void recomputeNormals ();
    void centerAndScaleToUnit ();
    void scaleUnit ();

    void draw() const {




    }
};



#endif
