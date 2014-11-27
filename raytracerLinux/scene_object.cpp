/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	// Followed pseudo-code from old tutorial notes here:
	// https://csc.cdf.toronto.edu/mybb/showthread.php?tid=8668

	// square vars
	Vector3D normal(0, 0, 1);

	// transform ray to model space
	ray.origin = worldToModel * ray.origin;
	ray.dir = worldToModel * ray.dir;

	// intersection position
	double t = -ray.origin[2] / ray.dir[2];

	// intersection must be in front of camera
	if (t <= 0) {
		return false;
	}

	// x/y values for intersection on xy-plane
	double x = ray.origin[0] + t*ray.dir[0];
	double y = ray.origin[1] + t*ray.dir[1];
	Point3D p(x, y, 0);

	// make sure intersection on xy-plane is in the defined shape
	if (x >= -0.5 && x <= 0.5 && y >= -0.5 && y <= 0.5) {
		// update if ray has no intersection, or
		// this intersection is closer to camera
		if (ray.intersection.none || t < ray.intersection.t_value) {
			ray.intersection.t_value = t;
			ray.intersection.point = modelToWorld * p;
			normal = worldToModel.transpose() * normal;
			normal.normalize();
			ray.intersection.normal = normal;
			ray.intersection.none = false;
			return true;
		}
	}

	return false;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	// transform ray to model space
	Point3D o = worldToModel * ray.origin;
	Vector3D origin(o[0], o[1], o[2]);
	Vector3D dir = worldToModel * ray.dir;
	dir.normalize();

	// algorithm from:
	// http://stackoverflow.com/questions/6533856/ray-sphere-intersection
	// https://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm

	double a = dir.dot(dir);
	double b = 2 * ( dir.dot(origin) );
	double c = origin.dot(origin) - 1;

	double delta = pow(b,2) - 4*a*c;

	if (delta < 0) {
		// no intersection
		return false;
	} else { 
		double t;
		if (delta == 0) {
			// single intersection (ray just touches the sphere)
			double t = -b / (2*a);
		} else {
			// two intersection points (ray passes through sphere)
			double t1 = (-b - sqrt(delta)) / (2*a);
			double t2 = (-b + sqrt(delta)) / (2*a);
			t = fmin(t1, t2); // pick the closer one non-negative one

		}

		if (t<= 0) {
			return false; 
		}

		// intersection pt
		Point3D p = o + (t * dir);

		// surface normal at intersection
		Vector3D normal(p[0], p[1], p[2]);
		normal.normalize();

		if (ray.intersection.none || t < ray.intersection.t_value) {
			ray.intersection.t_value = t;
			ray.intersection.point = modelToWorld * p;
			normal = worldToModel.transpose() * normal;
			normal.normalize();
			ray.intersection.normal = normal;
			ray.intersection.none = false;
			return true;
		}
	}
	
	return false;
}


bool Cylinder::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) { 

	//algo from https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html

	// transform ray to model space
	Point3D o = worldToModel * ray.origin;
	Vector3D origin(o[0], o[1], o[2]);
	Vector3D dir = worldToModel * ray.dir;
	dir.normalize();

	double min_z = -0.5; double max_z = 0.5; 

	double a = dir[0]*dir[0] + dir[1]*dir[1];
	double b = 2*(o[0]*dir[0]+o[1]*dir[1]);
	double c = o[0]*o[0] + o[1]*o[1] - 1; 

	double d = pow(b,2) - 4*a*c;

	if (d<0) {
		return false; 
	}	

	double t1 = (-b - sqrt(d)) / (2*a);
	double t2 = (-b + sqrt(d)) / (2*a);
	double t; 

	if (t1 <= 0 && t2 <= 0) {
		return false; 
	} else if (t1 <= 0) {
		t = t2; 
	} else if (t2 <= 0) {
		t = t1; 
	} else {
		t = fmin(t1, t2); 
	}

	//check if within range of z
	double x = o[0] + (t*dir[0]); 
	double y = o[1] + (t*dir[1]); 
	double z = o[2] + (t*dir[2]); 
	if (z>max_z || z<min_z) {
		return false; 
	}

	// intersection pt
	Point3D p(x, y, z);

	// surface normal at intersection
	Vector3D normal(p[0], p[1], 0);
	normal.normalize();

	if (ray.intersection.none || t < ray.intersection.t_value) {
		ray.intersection.t_value = t;
		ray.intersection.point = modelToWorld * p;
		normal = worldToModel.transpose() * normal;
		normal.normalize();
		ray.intersection.normal = normal;
		ray.intersection.none = false;
		return true;
	}

	return false; 
}
