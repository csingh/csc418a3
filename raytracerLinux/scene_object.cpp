/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include <limits>
#include "scene_object.h"

		#include <stdio.h>

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
	Ray3D r;
	r.origin = worldToModel * ray.origin;
	r.dir = worldToModel * ray.dir;

	// intersection position
	double t = -r.origin[2] / r.dir[2];

	// intersection must be in front of camera
	if (t <= 0) {
		return false;
	}

	// x/y values for intersection on xy-plane
	double x = r.origin[0] + t*r.dir[0];
	double y = r.origin[1] + t*r.dir[1];
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

	// check for cap intersection 
	double t_cap = std::numeric_limits<double>::max();  
	Vector3D cap_normal;
	bool intersectedCap = intersectCap(o, dir, t_cap, cap_normal);

	// checking for intersection of the sides of the cylinder

	double t_wall = std::numeric_limits<double>::max();  
	bool intersectedWall = intersectWall(o, dir, t_wall); 
	
	double t;
	if (!intersectedCap && !intersectedWall) {
		return false; 
	} else if (intersectedWall && intersectedCap) {
		t = fmin(t_wall, t_cap); 
	} else if (intersectedWall) {
		t = t_wall;
	} else {
		t = t_cap; 
	}

	// intersection pt
	Point3D p = o + t*dir;

	// surface normal at intersection
	Vector3D wall_normal(p[0], p[1], 0);
	
	Vector3D normal; 
	if (t == t_wall) {
		normal = wall_normal; 
	} else {
		normal = cap_normal; 
	}

	normal.normalize();

	if (ray.intersection.none || t < ray.intersection.t_value) {
		ray.intersection.t_value = t;
		ray.intersection.point = modelToWorld * p;
		normal = worldToModel.transpose() * normal;

		ray.intersection.normal = normal;
		ray.intersection.none = false;
		return true;
	}

	return false; 
}

bool Cylinder::intersectWall(Point3D origin, Vector3D direction, double& t) {

	double a = direction[0]*direction[0] + direction[1]*direction[1];
	double b = 2*(origin[0]*direction[0] + origin[1]*direction[1]);
	double c = origin[0]*origin[0] + origin[1]*origin[1] - 1; 

	double d = pow(b,2) - 4*a*c;

	if (d<=0) {
		return false; 
	}

	double t1 = (-b - sqrt(d)) / (2*a);
	double t2 = (-b + sqrt(d)) / (2*a);
	if (t1 <= 0 && t2 <= 0) {
		return false; 
	} else if (t1 <= 0 && t2 > 0) {
		t = t2; 
	} else if (t2 <=0 && t1 > 0) {
		t = t1; 
	} else {
		t = fmin(t1, t2); 
	} 

	double x = origin[0] + (t*direction[0]); 
	double y = origin[1] + (t*direction[1]); 
	double z = origin[2] + (t*direction[2]); 
	if (z > MAX_Z || z < MIN_Z) {
		return false; 
	}
	return true; 
}

bool Cylinder::intersectCap(Point3D origin, Vector3D direction, double& t, Vector3D &normal) {
		//check for interesecting caps
	
	Vector3D min_normal(0, 0, -1);
	Vector3D max_normal(0, 0, 1);
	double t1 = (MIN_Z - origin[2])/direction[2];
	double t2 = (MAX_Z - origin[2])/direction[2];

	if (t1<=0 && t2 <= 0) {
		return false; 
	}

	if (t1 <= 0 && t2 <= 0) {
		return false; 
	} else if (t1 <= 0 && t2 > 0) {
		t = t2; 
	} else if (t2 <=0 && t1 > 0) {
		t = t1; 
	} else if (t1 > 0 && t2 > 0) {
		t = fmin(t1, t2); 
	} 

	double x = origin[0] + t*direction[0]; 
	double y = origin[1] + t*direction[1]; 

	if (!(x*x + y*y <= 1)) {
		return false; 
	}
	
	if (t == t1) {
		normal = min_normal;
	} else {
		normal = max_normal;
	}

	return true; 
}

bool Cone::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld ) {
	
	// transform ray to model space
	Point3D o = worldToModel * ray.origin;
	Vector3D origin(o[0], o[1], o[2]);
	Vector3D d = worldToModel * ray.dir;
	d.normalize();

	double t_wall = std::numeric_limits<double>::max();
	Vector3D wall_normal;
	bool intersectedWall = intersectWall(o, d, t_wall, wall_normal);

	double t_cap = std::numeric_limits<double>::max();  
	Vector3D cap_normal; 
	bool intersectedCap = intersectCap(o, d, t_cap, cap_normal);

	Vector3D normal; 
	double t;
	if (!intersectedCap && !intersectedWall) {
		return false; 
	} else if (intersectedWall && intersectedCap) {
		if (t_wall < t_cap) {

			t= t_wall; 
			normal = wall_normal; 
		} else {

			t=t_cap; 
			normal = cap_normal; 
		}
		// t = fmin(t_wall, t_cap); 
	} else if (intersectedWall) {

		t = t_wall;
		normal = wall_normal;
	} else {
		t = t_cap; 
		normal = cap_normal;
	}

	// intersection pt
	Point3D p = o + t*d;

	if (t == t_wall) {
		normal = wall_normal; 
	} else {
		normal = cap_normal; 
	}
	
	normal.normalize();

	if (ray.intersection.none || t < ray.intersection.t_value) {
		ray.intersection.t_value = t;
		ray.intersection.point = modelToWorld * p;
		normal = worldToModel.transpose() * normal;

		ray.intersection.normal = normal;
		ray.intersection.none = false;
		return true;
	}

	return false; 

}

bool Cone::intersectWall(Point3D o, Vector3D d, double& t, Vector3D& n){

	double a = d[0]*d[0] + d[1]*d[1] - d[2]*d[2]; 
	double b = 2*(o[0]*d[0] + o[1]*d[1] - o[2]*d[2]);
	double c = o[0]*o[0] + o[1]*o[1] - o[2]*o[2]; 

	double delta = pow(b,2) - 4*a*c;

	if (delta<0) {
		return false; 
	}

	double t1 = (-b - sqrt(delta)) / (2*a);
	double t2 = (-b + sqrt(delta)) / (2*a);

	if (t1 <= 0 && t2 <= 0) {
		return false; 
	} 

	// check the point for the two roots returned
	double z1 = o[2]+t1*d[2];
	double z2 = o[2]+t2*d[2]; 
	Point3D p; 

	// check if z value is in range
	bool p1_inrange = (z1 < MAX_Z && z1 >= MIN_Z);
	bool p2_inrange = (z2 < MAX_Z && z2 >= MIN_Z);

	if(t1>0  && t2 > 0 && p1_inrange && p2_inrange){
		t = fmin(t1, t2); 
	} else if (t2 > 0 && p2_inrange) {
		t = t2; 
	} else if (t1 > 0 && p1_inrange) {
		t = t1;
	} else {
		return false; 
	}

	p = o+t*d; 
	n = Vector3D(2*p[0], 2*p[1], -2*p[2]);

	return true; 

}

bool Cone::intersectCap(Point3D o, Vector3D d, double& t, Vector3D& n) {
		
	Vector3D max_normal(0, 0, 1);
	t = (MAX_Z - o[2])/d[2];

	if (t < 0) {
		return false; 
	}

	Point3D	p = o + t*d; 

	if (!(p[0]*p[0] + p[1]*p[1] <= p[2]*p[2])) {
		return false; 
	}
	n = max_normal;

	return true; 
}
