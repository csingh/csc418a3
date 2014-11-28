/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray, int shadingType ) { 
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  


	// normalized vector from intersection point to light
	Vector3D light_ray = _pos - ray.intersection.point; 
	light_ray.normalize(); 

	Vector3D normal = ray.intersection.normal; 
	normal.normalize(); 

	// normalized vector in the reflection direction
	Vector3D reflect_dir = 2*(normal.dot(light_ray))*normal - light_ray; 
	reflect_dir.normalize(); 

	Material* mat = ray.intersection.mat; 

	Colour amb_ = mat->ambient*_col_ambient; 

	Vector3D camera = -ray.dir; 
	camera.normalize(); 

	double dif_comp = fmax(0.0, (light_ray.dot(normal))); 
	Colour dif_ = dif_comp * (mat->diffuse*_col_diffuse);
	
	double spec_comp = fmax(0.0, pow((reflect_dir.dot(camera)), mat->specular_exp));
	Colour spec_= spec_comp * (mat->specular*_col_specular); 

	if (shadingType == 1) {
		ray.col = amb_;
	} else {
		ray.col = amb_ + dif_ + spec_; 
	}
	ray.col.clamp(); 

}

