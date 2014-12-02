/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		classes defining primitives in the scene

***********************************************************/

#include "util.h"

// All primitives should provide a intersection function.  
// To create more primitives, inherit from SceneObject.
// Namely, you can create, Sphere, Cylinder, etc... classes
// here.
class SceneObject {
public:
	// Returns true if an intersection occured, false otherwise.
	virtual bool intersect( Ray3D&, const Matrix4x4&, const Matrix4x4& ) = 0;

};

// Example primitive you can create, this is a unit square on 
// the xy-plane.
class UnitSquare : public SceneObject {
public:
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
};

class UnitSphere : public SceneObject {
public:
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
};

class TexturedUnitSphere : public SceneObject {
public:
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );
};

// Cylinder aligned along the z-axis with height 1, centered at origin
class Cylinder : public SceneObject { 
	public:
	const static double MIN_Z = -0.5; 
	const static double MAX_Z = 0.5; 
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );

	// test if intersect cylinder wall (takes in transformed origin and direction)
	bool intersectWall(Point3D origin, Vector3D direction, double& t); 

	// test if intersect top or bottom (takes in transformed origin and direction)
	bool intersectCap(Point3D origin, Vector3D direction, double& t, Vector3D& normal); 
};

class Cone : public SceneObject { 
	public:
	const static double MIN_Z = 0; 
	const static double MAX_Z = 1.0; 
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld );

	bool intersectWall(Point3D origin, Vector3D direction, double& t_value, Vector3D& normal);
	
	// test if intersect the base of the cone (takes in transformed origin and direction)
	bool intersectCap(Point3D origin, Vector3D direction, double& t_value, Vector3D& normal); 
};
