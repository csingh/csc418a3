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
			Vector3D normal(0, 0, 1);
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
	
	return false;
}

