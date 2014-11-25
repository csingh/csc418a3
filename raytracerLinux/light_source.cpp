/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include <algorithm>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) { 
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	Vector3D light_ray = _pos - ray.intersection.point; 
	light_ray.normalize(); 
	
	Vector3D normal = ray.intersection.normal;
	normal.normalize();

	Vector3D reflect_dir = 2*(light_ray.dot(normal))*normal - light_ray; 
	reflect_dir.normalize();

	Vector3D camera = -ray.dir;
	camera.normalize();

	Material* mat = ray.intersection.mat;

	Colour c;

	Colour amb = mat->ambient + _col_ambient;
	Colour diff = fmax(0, light_ray.dot(normal)) * mat->diffuse * _col_diffuse;
	Colour spec = fmax( 0, pow(camera.dot(reflect_dir), mat->specular_exp) ) * mat->specular * _col_specular;

	c = c + amb + diff + spec;

	// double a_r, a_g, a_b, d_r, d_g,  d_b, s_r, s_g, s_b; 

	// a_r = mat->ambient[0]*_col_ambient[0];
	// a_g = mat->ambient[1]*_col_ambient[1];
	// a_b = mat->ambient[2]*_col_ambient[2];

	// Colour amb_(a_r, a_g, a_b);

	// d_r = mat->diffuse[0]*_col_diffuse[0]*std::max(0.0, (light_ray.dot(normal)));
	// d_g = mat->diffuse[1]*_col_diffuse[1]*std::max(0.0, (light_ray.dot(normal)));
	// d_b = mat->diffuse[2]*_col_diffuse[2]*std::max(0.0, (light_ray.dot(normal)));

	// Colour dif_(d_r, d_g, d_b); 

	// s_r = mat->specular[0]*_col_specular[0]*std::max(0.0, pow((reflect_dir.dot(ray.dir)),mat->specular_exp));
	// s_g = mat->specular[1]*_col_specular[1]*std::max(0.0, pow((reflect_dir.dot(ray.dir)),mat->specular_exp));
	// s_b = mat->specular[2]*_col_specular[2]*std::max(0.0, pow((reflect_dir.dot(ray.dir)),mat->specular_exp));

	// Colour spec_(s_r, s_g, s_b);

	// ray.col = amb_ + dif_+ spec_; 

	ray.col = ray.col + c;
	// ray.col.clamp();
}

