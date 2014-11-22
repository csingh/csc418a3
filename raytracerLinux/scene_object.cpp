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
			normal = modelToWorld.transpose() * normal;
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

	// sphere vars
	Point3D sphere_center(0, 0, 0);
	double sphere_r = 1; //radius

	// transform ray to model space
	ray.origin = worldToModel * ray.origin;
	ray.dir = worldToModel * ray.dir;

	// algorithm from:
	// http://stackoverflow.com/questions/6533856/ray-sphere-intersection

	// 2 points on ray
	double xa = ray.origin[0];
	double ya = ray.origin[1];
	double za = ray.origin[2];
	double xb = ray.origin[0] + ray.dir[0];
	double yb = ray.origin[1] + ray.dir[1];
	double zb = ray.origin[2] + ray.dir[2];

	double a = pow(xb-xa,2) + pow(yb-ya,2) + pow(zb-za,2);
	double b = 2 * ( (xb-xa)*(xa) + (yb-ya)*(ya) + (zb-za)*(za) );
	double c = pow(xa,2) + pow(ya,2) + pow(za,2) - pow(sphere_r, 2);

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
			t = (t1 < t2) ? t1 : t2; // pick the closer one
		}

		double x = xa + t*(xb-xa);
		double y = ya + t*(yb-ya);
		double z = za + t*(zb-za);
		Point3D p(x,y,z);
		if (ray.intersection.none || t < ray.intersection.t_value) {
			ray.intersection.t_value = t;
			ray.intersection.point = modelToWorld * p;
			//normal = modelToWorld.transpose() * normal;
			//normal.normalize();
			//ray.intersection.normal = normal;
			ray.intersection.none = false;
			return true;
		}
	}
	
	return false;
}

