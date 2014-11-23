/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) { 
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	Vector3D light_dir = _pos - ray.intersection.point; 
	light_dir.normalize(); 
	Vector3D normal = ray.intersection.normal; 
	normal.normalize(); 

	Vector3D reflect_dir = 2*(light_dir.dot(normal))*normal - light_dir; 
	reflect_dir.normalize(); 

	Material* mat = ray.intersection.mat; 

	double a_r, a_g, a_b, d_r, d_g,  d_b, s_r, s_g, s_b; 

	a_r = mat->ambient[0]*_col_ambient[0];
	a_g = mat->ambient[1]*_col_ambient[1];
	a_b = mat->ambient[2]*_col_ambient[2];

	Colour amb_(a_r, a_g, a_b);

	d_r = mat->diffuse[0]*(light_dir.dot(normal))*_col_diffuse[0];
	d_g = mat->diffuse[1]*(light_dir.dot(normal))*_col_diffuse[1];
	d_b = mat->diffuse[2]*(light_dir.dot(normal))*_col_diffuse[2];

	Colour dif_(d_r, d_g, d_b); 

	s_r = mat->specular[0]*pow((reflect_dir.dot(-ray.dir)),mat->specular_exp)*_col_specular[0];
	s_g = mat->specular[1]*pow((reflect_dir.dot(-ray.dir)),mat->specular_exp)*_col_specular[1];
	s_b = mat->specular[2]*pow((reflect_dir.dot(-ray.dir)),mat->specular_exp)*_col_specular[2];

	Colour spec_(s_r, s_g, s_b);

	ray.col = amb_ + dif_ + spec_; 

}

