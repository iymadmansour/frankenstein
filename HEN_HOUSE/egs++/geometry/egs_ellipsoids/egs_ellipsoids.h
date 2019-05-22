#ifndef EGS_CELLIPSOIDS_
#define EGS_CELLIPSOIDS_

#include "egs_base_geometry.h"

#ifdef WIN32

#ifdef BUILD_ELLIPSOIDS_DLL
#define EGS_ELLIPSOIDS_EXPORT __declspec(dllexport)
#else
#define EGS_ELLIPSOIDS_EXPORT __declspec(dllimport)
#endif
#define EGS_ELLIPSOIDS_LOCAL 

#else

#ifdef HAVE_VISIBILITY
#define EGS_ELLIPSOIDS_EXPORT __attribute__ ((visibility ("default")))
#define EGS_ELLIPSOIDS_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define EGS_ELLIPSOIDS_EXPORT
#define EGS_ELLIPSOIDS_LOCAL
#endif

#endif

class EGS_ELLIPSOIDS_EXPORT EGS_cEllipsoids : public EGS_BaseGeometry
{     
public:

    // construct some ellipsoid
    EGS_cEllipsoids(const EGS_Float &a, const EGS_Float &b, const EGS_Float &c, const EGS_Vector &position, const string &Name = "");

    // destruct ellipsoid from memory
    ~EGS_cEllipsoids()
	{}

    // method to determine which ellipsoids we are in(between)
    int inside(const EGS_Vector &x);
    bool isInside(const EGS_Vector &x);
    int isWhere(const EGS_Vector &x);

    EGS_Float howfarToOutside(int ireg, const EGS_Vector &x,
                    const EGS_Vector &u);
    // howfar is particle trajectory from sphere boundry
    int howfar(int ireg, const EGS_Vector &x, const EGS_Vector &u, 
            EGS_Float &t, int *newmed=0, EGS_Vector *normal=0);

    // hownear - closest perpendicular distance to sphere surface
    EGS_Float hownear(int ireg, const EGS_Vector &x);

    int getMaxStep() const { return 2*nreg; };

    const string &getType() const { return type; };

    void printInfo() const;

private:

	EGS_Float A, B, C; // x, y and z semi-axis
	EGS_Float max, min; // largest and smallest semi-axis
	EGS_Float A2, B2, C2; // x, y and z semi-axis squared
	EGS_Float max2, min2; // largest and smallest semi-axis squared
	// such that, (x/A)^2+(y/B)^2+(z/C)^2=1
	
	EGS_Vector xo; // midpoint
    static string type;
};

#endif
