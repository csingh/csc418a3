/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		   light source classes

***********************************************************/

#include "util.h"

// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray 
// differently.
class LightSource {
public:
	virtual void shade( Ray3D&, int shadingType ) = 0;
	virtual Point3D get_position() = 0; 
};

// A point light is defined by its position in world space and its
// colour.
class PointLight : public LightSource {
public:
	PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col), 
	_col_diffuse(col), _col_specular(col) {}
	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ) 
	: _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
	_col_specular(specular) {}
	void shade( Ray3D& ray, int shadingType );
	Point3D get_position() { return _pos; }
	
private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular;
};

// An area light is defined as a plane, given its center point, two vectors
// that define its plane in the world, and its colour
class AreaLight : public LightSource {
public:
	AreaLight( Point3D pos, Vector3D u, Vector3D v, Colour col ) : _pos(pos), _u(u), _v(v), _col_ambient(col), 
	_col_diffuse(col), _col_specular(col) {}
	AreaLight( Point3D pos, Vector3D u, Vector3D v, Colour ambient, Colour diffuse, Colour specular ) 
	: _pos(pos), _u(u), _v(v), _col_ambient(ambient), _col_diffuse(diffuse), 
	_col_specular(specular) {}
	void shade( Ray3D& ray, int shadingType );
	Point3D get_position();
	Vector3D get_u() const { return _u; }
	Vector3D get_v() const { return _v; }
private:
	Point3D _pos;
	Vector3D _u;
	Vector3D _v;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular;
};
