#include "egs_ellipsoids.h"
#include "egs_input.h"
#include "egs_functions.h"

#include <vector>
using std::vector;

string EGS_cEllipsoids::type = "EGS_cEllipsoids";

// generate the concentric spheres
EGS_cEllipsoids::EGS_cEllipsoids(const EGS_Float &a, const EGS_Float &b, const EGS_Float &c, const EGS_Vector &position, const string &Name) : EGS_BaseGeometry(Name), xo(position)
{
    A=a;
	A2=A*A;
	B=b;
	B2=B*B;
	C=c;
	C2=C*C;
	min=A<B?(A<C?A:C):(B<C?B:C);
	min2=min*min;
	max=A>B?(A>C?A:C):(B>C?B:C);
	max2=max*max;
    nreg=1;
}

bool EGS_cEllipsoids::isInside(const EGS_Vector &x)
{
    EGS_Vector tmp(x-xo);
	EGS_Float rad = tmp.x*tmp.x/A2+tmp.y*tmp.y/B2+tmp.z*tmp.z/C2;
    if(rad > 1)
		return false;
    return true;
}

int EGS_cEllipsoids::isWhere(const EGS_Vector &x)
{
    EGS_Vector tmp(x-xo);
	EGS_Float rad = tmp.x*tmp.x/A2+tmp.y*tmp.y/B2+tmp.z*tmp.z/C2;
    if(rad > 1)
		return -1;
    return 0;
}

// method to determine which spheres we are in(between)
int EGS_cEllipsoids::inside(const EGS_Vector &x)
{
    EGS_Vector tmp(x-xo);
	EGS_Float rad = tmp.x*tmp.x/A2+tmp.y*tmp.y/B2+tmp.z*tmp.z/C2;
    if(rad > 1)
		return -1;
    return 0;
}

// howfar is particle trajectory from sphere boundary
/* note that in general we will be between two spheres (if inside a sphere at
 * all... so we need to check if the flight path will intersect the inner or
 * outer sphere
 */

#ifdef ELLIPSOIDS_DEBUG
EGS_Vector last_x, last_u;
int last_ireg;
EGS_Float last_d,last_t,last_aa,last_bb2,last_R2b2,last_tmp;
#endif

EGS_Float EGS_cEllipsoids::howfarToOutside(int ireg, const EGS_Vector &x, const EGS_Vector &u)
{
	EGS_Vector nx = x-xo, nu = u;
	EGS_Float Qa =    nu.x*nu.x/A2 + nu.y*nu.y/B2 + nu.z*nu.z/C2;     // The a in the quadratic formula
	EGS_Float Qb = 2*(nx.x*nu.x/A2 + nx.y*nu.y/B2 + nx.z*nu.z/C2);    // The b in the quadratic formula
	EGS_Float Qc =    nx.x*nx.x/A2 + nx.y*nx.y/B2 + nx.z*nx.z/C2 - 1; // The c in the quadratic formula
	
	EGS_Float QR = Qb*Qb - 4.0*Qa*Qc; // The radical in the quadratic formula
	
	if (!(QR >= 0)) // Imaginary roots, ie, no intersection (this should never be true if ireg >= 0)
		return -1;
	
	// m1 and m2 are the constants by which you must multiply u to get a vector that goes from x to the surface of the ellipsoid
	EGS_Float m1 = (-Qb - sqrt(QR))/(2.0*Qa);
	EGS_Float m2 = (-Qb + sqrt(QR))/(2.0*Qa);
	if (m1 < 0.0000001) // Round it off to zero
		m1 = 0;
	if (m2 < 0.0000001) // Round it off to zero
		m2 = 0;
	
	EGS_Float m; // The m yielding the closer intersection
	EGS_Float d; // The distance to the closer intersection
	EGS_Vector temp (nu); // Vector from point to intersection
	
	if (ireg < 0) // Outside of the ellipsoid
	{
		if (m1 < 0) // We are past it (m1 and m2 should have the same sign)
			return -1;
		
		m = m1*m1<m2*m2?m1:m2;
		temp *= m;
		
		return temp.length();
	}
	else // Inside of the ellipsoid
	{
		m = m1>0?m1:m2; // One scalar must be negative, therefore the other one is positive
		temp *= m;
		
		return temp.length();
	}
}

int EGS_cEllipsoids::howfar(int ireg, const EGS_Vector &x, const EGS_Vector &u, EGS_Float &t, int *newmed, EGS_Vector *normal)
{
	EGS_Vector nx = x-xo, nu = u;
	//cout << x.x << "," << x.y << "," << x.z << ";" << u.x << "," << u.y << "," << u.z << "\n"; cout.flush();
	EGS_Float Qa =    nu.x*nu.x/A2 + nu.y*nu.y/B2 + nu.z*nu.z/C2;     // The a in the quadratic formula
	EGS_Float Qb = 2*(nx.x*nu.x/A2 + nx.y*nu.y/B2 + nx.z*nu.z/C2);    // The b in the quadratic formula
	EGS_Float Qc =    nx.x*nx.x/A2 + nx.y*nx.y/B2 + nx.z*nx.z/C2 - 1; // The c in the quadratic formula
	
	EGS_Float QR = Qb*Qb - 4.0*Qa*Qc; // The radical in the quadratic formula
	
	//cout << Qa << "," << Qb << "," << Qc << ";" << QR << "\n"; cout.flush();
	if (!(QR >= 0)) // Imaginary roots, ie, no intersection (this should never be true if ireg >= 0)
	{
		//cout << "Miss\n\n"; cout.flush();
		return -1;
	}
	
	// m1 and m2 are the constants by which you must multiply u to get a vector that goes from x to the surface of the ellipsoid
	QR=sqrt(QR);
	EGS_Float m1 = (-Qb - QR)/(2.0*Qa);
	EGS_Float m2 = (-Qb + QR)/(2.0*Qa);
	if (fabs(m1) < 0.0000001) // Round it off to zero
		m1 = 0;
	if (fabs(m2) < 0.0000001) // Round it off to zero
		m2 = 0;
	//cout << m1 << "," << m2 << "\n"; cout.flush();
	
	EGS_Float m; // The m yielding the closer intersection
	EGS_Float d; // The distance to the closer intersection
	EGS_Vector temp (nu); // Vector from point to intersection
	
	if (ireg < 0) // Outside of the ellipsoid
	{
		if (m1 < 0 || m2 < 0) // We are past it (m1 and m2 should have the same sign unless one is zero)
		{
			//cout << "Past it\n\n"; cout.flush(); sleep(1);
			return -1;
		}
		
		m = m1*m1<m2*m2?m1:m2;
		temp *= m;
		//cout << temp.x << "," << temp.y << "," << temp.z << "\n"; cout.flush();
		
		d = temp.length();
		//cout << d << "," << t << "\n\n"; cout.flush(); sleep(1);
		
		if (d <= t) // Does it make it?
		{
			t = d;
			if (newmed)
				*newmed = medium(0);
			if (normal)
			{
				*normal = EGS_Vector(2/A2*(nx.x+temp.x),2/B2*(nx.y+temp.y),2/C2*(nx.z+temp.z));
				normal->normalize();
			}
			return 0;
		}
		else
			return -1;
	}
	else // Inside of the ellipsoid
	{
		m = m1>0?m1:m2; // One scalar must be negative, therefore the other one is positive
		temp *= m;
		//cout << temp.x << "," << temp.y << "," << temp.z << "\n"; cout.flush();
		
		d = temp.length();
		//cout << d << "\n\n"; cout.flush();
		
		if (d <= t) // Does it make it?
		{
			t = d;
			if (newmed)
				*newmed = medium(-1);
			if (normal)
			{
				*normal = EGS_Vector(2/A2*(nx.x+temp.x),2/B2*(nx.y+temp.y),2/C2*(nx.z+temp.z));
				normal->normalize();
			}
			return -1;
		}
		else
			return 0;
	}
}

// hownear - closest perpendicular distance to sphere surface
EGS_Float EGS_cEllipsoids::hownear(int ireg, const EGS_Vector &x)
{
	EGS_Vector d = (x-xo);
	EGS_Float d2 = d.length2();
	EGS_Vector n = d; n.normalize();
	// Far and near conditions for easy outs
	if (d2 >= max2)
		return sqrt(d2)-max;
	else if (d2 <= min2)
		return min-sqrt(d2);
	else
		return (d-EGS_Vector(n.x*A2,n.y*B2,n.z*C2)).length()/(A2*B2*C2);
	// The above is a huge underestimate of the closest point on the ellipse,
	// but the increased step count should still take less time than finding
	// a more exact solution.
}

void EGS_cEllipsoids::printInfo() const
{
    EGS_BaseGeometry::printInfo();
    egsInformation(" midpoint of ellipsoid = (%g,%g,%g)\n",xo.x,xo.y,xo.z);
    egsInformation(" ellipsoid equation = x^2/%g^2 + y^2/%g^2 + z^2/%g^2 = 1",A,B,C);
    egsInformation("\n=======================================================\n");
}

extern "C"
{
EGS_ELLIPSOIDS_EXPORT EGS_BaseGeometry* createGeometry(EGS_Input *input)
{
    if( !input ) {
        egsWarning("createGeometry(ellipsoids): null input?\n");
        return 0;
    }
    EGS_Vector xo; vector<EGS_Float> Xo;
    int err = input->getInput("midpoint",Xo);
    if( !err && Xo.size() == 3 )
		xo = EGS_Vector(Xo[0],Xo[1],Xo[2]);
		
	EGS_Float xs;
    err = input->getInput("x semi-axis",xs);
    if( err || xs <= 0)
	{
        egsWarning("createGeometry(ellipsoids): wrong/missing 'x semi-axis' input\n");
        return 0;
    }
	
	EGS_Float ys;
    err = input->getInput("y semi-axis",ys);
    if(err  || ys <= 0)
	{
        egsWarning("createGeometry(ellipsoids): wrong/missing 'y semi-axis' input\n");
        return 0;
    }
	
	EGS_Float zs;
    err = input->getInput("z semi-axis",zs);
    if(err || zs <= 0)
	{
        egsWarning("createGeometry(ellipsoids): wrong/missing 'z semi-axis' input\n");
        return 0;
    }
	
	EGS_BaseGeometry *result = new EGS_cEllipsoids(xs,ys,zs,xo);
    result->setName(input); result->setMedia(input);
    return result;
}
}
