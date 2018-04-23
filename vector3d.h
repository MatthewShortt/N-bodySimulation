#ifndef vec3_h
#define vec3_h

#include <cmath>

class vec3
{
public:
	// Data
	double x, y, z;
	int r, g, b;
	double mass, velocity;
    int direction;
		

	// Ctors
	vec3( double InX, double InY, double InZ, int InR, int InG, int InB, double InMass, double InVelocity, int InDirection ) : x( InX ), y( InY ), z( InZ )
		{
			r=InR;
			g=InG;
			b=InB;
			mass=InMass;
			velocity=InVelocity;
            direction = InDirection;
		}

	void SetDoublePoint( const double *v ) {  x=v[0]; y=v[1]; z=v[2]; }

	vec3( ) : x(0), y(0), z(0)
		{
		}

	// Operator Overloads
	inline bool operator== (const vec3& V2) const 
		{
		return (x == V2.x && y == V2.y && z == V2.z);
		}

	inline vec3 operator+ (const vec3& V2) const 
		{
		return vec3( x + V2.x,  y + V2.y,  z + V2.z, r, g, b, mass, velocity, direction);
		}

	inline vec3 operator- (const vec3& V2) const
		{
		return vec3( x - V2.x,  y - V2.y,  z - V2.z, r, g, b, mass, velocity, direction);
		}
	inline vec3 SubP(const double *v) const
		{
		  return vec3( x - v[0],  y - v[1],  z - v[2], r, g, b, mass, velocity, direction);
		}

	inline vec3 operator- ( ) const
		{
		return vec3(-x, -y, -z, r, g, b, mass, velocity, direction);
		}

	inline vec3 operator/ (double S ) const
		{
		double fInv = 1.0 / S;
		return vec3 (x * fInv , y * fInv, z * fInv, r, g, b, mass, velocity, direction);
		}
	inline vec3 operator/ (const vec3& V2) const
		{
		return vec3 (x / V2.x,  y / V2.y,  z / V2.z, r, g, b, mass, velocity, direction);
		}
	inline vec3 operator* (const vec3& V2) const
		{
		return vec3 (x * V2.x,  y * V2.y,  z * V2.z, r, g, b, mass, velocity, direction);
		}
	inline vec3 operator* (double S) const
		{
		return vec3 (x * S,  y * S,  z * S, r, g, b, mass, velocity, direction);
		}
	inline vec3 operator+ (double S) const
		{
		return vec3 (x + S,  y + S,  z + S, r, g, b, mass, velocity, direction);
		}
	inline vec3 operator- (double S) const
		{
		return vec3 (x - S,  y - S,  z - S, r, g, b, mass, velocity, direction);
		}

	inline void operator+= ( const vec3& V2 )
		{
		x += V2.x;
		y += V2.y;
		z += V2.z;
		}
	inline void operator-= ( const vec3& V2 )
		{
		x -= V2.x;
		y -= V2.y;
		z -= V2.z;
		}

	inline double operator[] ( int i )
		{
		if ( i == 0 ) return x;
		else if ( i == 1 ) return y;
		else return z;
		}

	// Functions
	inline int getR()
		{
		return r;
		}
	inline int getG()
		{
		return g;
		}
	inline int getB()
		{
		return b;
		}
	inline double getX()
		{
		return x;
		}
	inline double getY()
		{
		return y;
		}
	inline double getZ()
		{
		return z;
        }
    inline void setX(double thetaX)
    {
        x = thetaX;
    }
    inline void setY(double thetaY)
    {
        y = thetaY;
    }
    inline void setZ(double thetaZ)
    {
        z = thetaZ;
    }
	inline double getMass()
		{
		return mass;
		}
	inline int getDirection()
		{
		return direction;
        }
    
    inline void setPosition(vec3& part, double timeSubStep)
        {
            switch(part.getDirection()) {
                case 0 :
                    //changeX = part.getX() + timeSubStep*part.getVelocity();
                    //printf("Original X: %f \tchangeX: %f\t\t", part.getX(), changeX);
                    part.setX(part.getX() + timeSubStep*part.getVelocity());
                    //printf("NEW X: %f \t", part.getX());
                    break;
                case 1 :
                    part.setY(part.getY() + timeSubStep*part.getVelocity());
                    break;
                case 2 :
                    part.setZ(part.getZ() + timeSubStep*part.getVelocity());
                    break;
                case 3 :
                    part.setX(part.getX() + (timeSubStep*part.getVelocity())/sqrt(2));
                    part.setY(part.getY() + (timeSubStep*part.getVelocity())/sqrt(2));
                    break;
                case 4 :
                    part.setX(part.getX() + (timeSubStep*part.getVelocity())/sqrt(2));
                    part.setZ(part.getZ() + (timeSubStep*part.getVelocity())/sqrt(2));
                    break;
                case 5 :
                    part.setY(part.getY() + (timeSubStep*part.getVelocity())/sqrt(2));
                    part.setZ(part.getZ() + (timeSubStep*part.getVelocity())/sqrt(2));
                    break;
                case 6 :
                    part.setX(part.getX() + (timeSubStep*part.getVelocity())/sqrt(2));
                    part.setY(part.getY() + (timeSubStep*part.getVelocity())/sqrt(2));
                    part.setZ(part.getZ() + (timeSubStep*part.getVelocity())/sqrt(2));
                    break;
                default :
                    part.setX(part.getX() + timeSubStep*part.getVelocity());
                    
            }
        }
    inline double getVelocity()
        {
        return velocity;
        }
	inline void setVelocity(vec3& part, double timeSubStep, double force)
		{
		velocity = part.getVelocity() + ((timeSubStep*force)/part.getMass()) ;
		}
	
	inline double Dot( const vec3 &V1 ) const
		{
		return V1.x*x + V1.y*y + V1.z*z;
		}

	// These require math.h for the sqrt function
	double Magnitude( ) const
		{
		return sqrt( x*x + y*y + z*z );
		}
    double MagnitudeSquared( ) const
        {
        return x*x + y*y + z*z;
        }

	inline void Normalize()
		{
		double fMag = ( x*x + y*y + z*z );
		if (fMag == 0) {return;}

		double fMult = 1.0/sqrt(fMag);            
		x *= fMult;
		y *= fMult;
		z *= fMult;
		return;
		}

};


inline vec3 SubtractDoubleDouble(const double *d1, const double *d2)
{
  return vec3(d1[0]-d2[0], d1[1]-d2[1], d1[2]-d2[2], '0', '0', '0', 0, 0, 0);
}

inline double clamp(double d, double min, double max)
{
  if (d < min)
    return min;
  if (d > max)
    return max;
  return d;
}

inline void clamp(vec3 &v, double min, double max) 
{
  v.x = clamp(v.x,min,max);
  v.y = clamp(v.y,min,max);
  v.z = clamp(v.z,min,max);
}



#endif
