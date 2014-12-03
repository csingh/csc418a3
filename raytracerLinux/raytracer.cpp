/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel; 
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
	for (;;) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.

		int i = 0;
		for (i = 0; i < NUM_SHADOW_RAYS; i++) {
			// Implement shadows here if needed.
			Point3D int_point = ray.intersection.point;
			Point3D light_p = curLight->light->get_position();
			Vector3D dir = light_p - int_point;
			dir.normalize();
			// offset point slighty towards light, otherwise
			// all rays will intersect object at t = 0
			Point3D p = int_point + (0.01 * dir);

			Ray3D shadowRay(p, dir);

			traverseScene(_root, shadowRay);

			// if intersection happened between intersection point and the light source
			// TODO: clean up (somewhat of a mess)
			bool do_as_normal = true;
			if ( !shadowRay.intersection.none ) {
				// shadow ray hit something
				double dist_to_light = (light_p - int_point).length();
				double dist_between_points = (shadowRay.intersection.point - ray.intersection.point).length();
				if (dist_between_points < dist_to_light) {
					// shadow ray hit something, and its between point and light, so light is being blocked
					curLight->light->shade(ray, 1);
					do_as_normal = false;
				}
			}

			if (do_as_normal) {
				// shadow ray didnt hit anything, compute pixel color as normal
				curLight->light->shade(ray, 0);				
			}
		}

		if (i == 0) {
			curLight->light->shade(ray, 0);
		} else {
			ray.col = (1.0 / i) * ray.col;
		}

		ray.col.clamp();

		curLight = curLight->next;
	}
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray ) {
	Colour col(0.0, 0.0, 0.0);
	traverseScene(_root, ray);

	// Don't bother shading if the ray didn't hit 
	// anything.
	if (ray.intersection.none) {
		return col;
	}

	// --- TEXTURE CODE ---
	// TODO: move to computeShading? Yes if we want to add lighting to textures
	if (ray.intersection.mat->is_texture) {
		int w = ray.intersection.mat->texture_width;
		int h = ray.intersection.mat->texture_height;
		double u = ray.texture_u;
		double v = ray.texture_v;
		int x = u * w; // x coordinate of pixel
		int y = v * h; // y coordinate of pixel
		// printf("w, h: %d, %d\n", w, h);
		// printf("u, v: %f, %f\n", u, v);
		// printf("x, y: %d, %d\n", x, y);
		int i = (x * w) + y; // index into rbgarrays
		// printf("i: %d\n", i);
		double r = (1.0/256) * ray.intersection.mat->rarray[i];
		double g = (1.0/256) * ray.intersection.mat->garray[i];
		double b = (1.0/256) * ray.intersection.mat->barray[i];
		// printf("(r,g,b): (%f, %f, %f)\n", r, g, b);
		col = Colour( r, g, b );
		return col;
	}

	// You'll want to call shadeRay recursively (with a different ray, 
	// of course) here to implement reflection/refraction effects.

	computeShading(ray); 
	col = ray.col;

	if ( ray.num_reflections < NUM_REFLECTION_RAYS) {

		Colour reflectCol(0, 0, 0); 
		Colour refractCol(0, 0, 0); 

		// --- REFLECTION CODE ----
		// Set up incident ray
		Point3D p = ray.intersection.point;
		Vector3D normal(ray.intersection.normal);
		Vector3D dir = ray.dir - 2*( normal.dot(ray.dir) * normal );
		dir.normalize();
		p = p + (0.01 * dir);

		// Find vectors u,v that make a plane perpendular to incident ray
		Vector3D u = get_orthonormal_vector(dir);
		Vector3D v = dir.cross(u);
		v.normalize();

		if (ray.intersection.mat->reflectance > 0) {
			int i = 0;
			Ray3D gloss_ray;
			for (i = 0; i < NUM_GLOSSY_REFLECTION_RAYS; i++) {
				double u_mult = rand_range(-GLOSS_AMOUNT, GLOSS_AMOUNT);  
				double v_mult = rand_range(-GLOSS_AMOUNT, GLOSS_AMOUNT);

				Point3D glossy_ray_point = p + (u_mult * u) + (v_mult * v);
				// Vector3D glossy_ray_dir = dir + (u_mult * u) + (v_mult * v);
				Vector3D glossy_ray_dir = dir;
				glossy_ray_dir.normalize();

				gloss_ray = Ray3D(glossy_ray_point, glossy_ray_dir);
				// TODO: check (don't think we need this mat line)
				gloss_ray.intersection.mat = ray.intersection.mat;
				// computeShading(gloss_ray);
				// col = col + gloss_ray.col;

				gloss_ray.num_reflections = ray.num_reflections + 1;
				reflectCol = reflectCol + ray.intersection.mat->reflectance * shadeRay(gloss_ray);

			}

			if (i > 0) {
				reflectCol = ( 1.0 / (i*ray.intersection.mat->reflectance) ) * reflectCol;
			}

			reflectCol.clamp();
		}

		// --- REFRACTION CODE ---
		if (ray.intersection.mat->refractance > 0 ) {
			Vector3D normal(ray.intersection.normal);
		
			normal.normalize();
			
			Vector3D rdir = ray.dir;
			rdir.normalize(); 
			
		
			// check if leaving or entering material by checking normal
			double n1 = ray.refrac_ind;
			double n2 = ray.intersection.mat->refracive_ind; 
			double cosI = rdir.dot(normal); 
			if  (cosI > 0) {
				// ray is inside material
				n1 = ray.intersection.mat->refracive_ind; 
				n2 = ray.refrac_ind; 

				normal = -normal; 
			
			} else {
				// ray is outside material
				cosI = - cosI;
			}
		
			// https://www.cs.unc.edu/~rademach/xroads-RT/RTarticle.html
			double n = n1/n2; 
			double cosT = 1.0-pow(n,2)*(1.0-pow(cosI,2)); 
			
			if (cosT >= 0.0)
			{

				Vector3D refractDir = (n*rdir) + (n*cosI -sqrt(cosT))*normal; 
				refractDir.normalize(); 
				
				Point3D newOrigin = ray.intersection.point  + 0.01*refractDir;
				Ray3D refractRay(newOrigin, refractDir);
				
				// get the exit point of the ray and
				refractCol= Colour(ray.intersection.mat->refractance*shadeRay(refractRay));
				refractCol.clamp();				
			}
		}
		col = col + refractCol + reflectCol; 
		col.clamp(); 
	}

	return col;

}	

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);

	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);

	// float min_x = 0, max_x = 0, min_y = 0, max_y = 0;

	bool firstCall = true;
	// Construct a ray for each pixel.
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {

			// dividing the pixel into n*n subpixels for anti-aliasing
			int n = NUM_AA_SAMPLES; 
			double scale = 1.0/(double)n; 
			Colour col(0, 0, 0); 
			for (int m = 0; m < n*n; m++) {
				int y = (int) m / n; 
				int x = m % n; 

				// to sample a random x and y in this sub pixel 
				double rand_x = rand_range(x*scale, (x+1)*scale);  
				double rand_y = rand_range(y*scale, (y+1)*scale); 

				double x_ = (-double(width)/2 + rand_x + j)/factor;
				double y_ = (-double(height)/2 + rand_y + i)/factor; 
				
				// Sets up ray origin and direction in view space, 
				// image plane is at z = -1.
				Point3D origin(0, 0, 0);
				Point3D imagePlane;
				imagePlane[0] = x_;
				imagePlane[1] = y_;
				imagePlane[2] = -1;
				
				origin = viewToWorld * origin;
				Vector3D dir(imagePlane[0], imagePlane[1], imagePlane[2]);
				dir = viewToWorld * dir;
				dir.normalize();

				Ray3D ray(origin, dir);

				Colour subcol = shadeRay(ray);
				
				col = col + subcol; 
			}

			_rbuffer[i*width+j] = int((col[0]/pow(n,2))*255);
			_gbuffer[i*width+j] = int((col[1]/pow(n,2))*255);
			_bbuffer[i*width+j] = int((col[2]/pow(n,2))*255);
			
		}

		loadBar(i, _scrHeight, 50, 50, firstCall);

		if (i == 0) firstCall = false;

	}

	flushPixelBuffer(fileName);
}

int main(int argc, char* argv[])
{	
	//setting up the randomizer
	srand((unsigned)time(NULL));

	int width = 320; 
	int height = 240; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// original_scene(width, height); 
	// scene_part_b_cylinder_cone(width, height);
	// scene_single_textured_sphere(width, height); 
	refraction_reflection_scene(width, height); 
	return 0;
}

// from: https://www.ross.click/2011/02/creating-a-progress-bar-in-c-or-any-other-console-app/
// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
static inline void loadBar(int x, int n, int r, int w, bool firstCall)
{

	//Make sure that load bar is displayed on next line and don’t delete current console line during first run
	if (firstCall)
	{
	    std::cout << std::endl;
	    firstCall = false;
	}

    // Only update r times.
    if ( x % (n/r) != 0 ) return;
 
    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    for (int x=0; x<c; x++)
       printf("=");
 
    for (int x=c; x<w; x++)
       printf(" ");
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
    // printf("]\n33[F33[J");
printf("\r"); // Move to the first column
fflush(stdout);
}

void original_scene(int width, int height) { 
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it frm two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;

	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2, 1.0 , 0.0, 0.0);
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.0 , 0.0, 0.0);

	Material silver( Colour(0.2, 0.2, 0.2), Colour(0.50754, 0.50754, 0.50754),
			Colour(0.508273, 0.508273, 0.508273), 
			51.2, 1.0 , 0.0, 0.0);

	Material glass( Colour(0.0, 0.0, 0.0), Colour(0.1, 0.1, 0.1), 
			Colour(0.5, 0.5, 0.5),
			0.25, 0.0 , 1.0, .9); 


	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(10, 10, 3), 
				Colour(0.9, 0.9, 0.9) ) );

	// Add a unit square into the scene with material mat.
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &jade );
	SceneDagNode* plane2 = raytracer.addObject( new UnitSquare(), &jade );
	SceneDagNode* cylinder = raytracer.addObject(new Cone(), &glass);
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &glass );

	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	double factor3[3] = { 2.0, 1.5, 2.0};

	double factor4[3] = {0.5, 0.5, 2.0};
	raytracer.translate(sphere, Vector3D(1, -0.75, -2));	

	raytracer.translate(cylinder, Vector3D(-1, 1, -2));	
	raytracer.scale(cylinder, Point3D(0, 0, 0), factor4);
	raytracer.rotate(cylinder, 'y', 180);

	// raytracer.rotate(sphere, 'x', -45); 
	// raytracer.rotate(sphere, 'z', 45); 
	// raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane, Vector3D(0, 0, -7));	
	// raytracer.rotate(plane, 'z', 90); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane2, Vector3D(-4, 0, -4));	
	raytracer.rotate(plane2, 'y', 90); 
	// raytracer.rotate(plane2, 'z', 90);
	raytracer.scale(plane2, Point3D(0, 0, 0), factor2);
	
	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	printf("Rendering image 1...");
	raytracer.render(width, height, eye, view, up, fov, "images/view1.bmp");
	
	// Render it from a different point of view.
	// Point3D eye2(5, 0, 0);
	// Vector3D view2(-5, 0, -7);
	Point3D eye2(4, -1, 0);
	Vector3D view2(-8, 0, -5);
	printf("Rendering image 2...");
	raytracer.render(width, height, eye2, view2, up, fov, "images/view2.bmp");

}

void scene_part_b_cylinder_cone(int width, int height){

	Raytracer raytracer;

	// Camera parameters.
	Point3D eye(0, 30, 1);
	Vector3D view(0, -17, -20);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
	Colour(0.628281, 0.555802, 0.366065), 
	51.2, 1.0, 0.0, 0.0);
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
	Colour(0.316228, 0.316228, 0.316228), 
	12.8, 1.0, 0.0, 0.0);
	Material randomCol( Colour(.6, .1, 0), Colour(0.12, 0.89, 0.9), 
	Colour(0.1, 0.9, 0.5), 
	51.2, 1.0, 0.0, 0.0);
	Material chrome( Colour(.25, .25, .25), Colour(0.4	,0.4,	0.4), 
	Colour(0.774597,	0.774597,	0.774597), 
	76.8, 1.0, 0.0, 0.0);

	Material glass( Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0), 
		Colour(0.7, 0.7, 0.7),
		0.25, 0.0 , 1.5, 1.0); 
	// gold.reflective=false;
	// jade.reflective=false;
	// chrome.reflective=false;
	// randomCol.reflective=true;

	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(5, 2, -5), 
	Colour(0.9, 0.9, 0.9) ) );
	// raytracer.get_light_list_node()->light->set_render(ren);
	// raytracer.get_light_list_node()->light->sample=1;
	// Add a unit square into the scene with material mat.
	SceneDagNode* cylinder = raytracer.addObject( new Cylinder(), &glass );
	SceneDagNode* cone = raytracer.addObject( new Cone(), &glass );
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &randomCol);
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &jade );
	SceneDagNode* plane2 = raytracer.addObject( new UnitSquare(), &jade );
	SceneDagNode* plane3 = raytracer.addObject( new UnitSquare(), &jade );
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 3.0, 3.0, 8.0 };
	double factor2[3] = { 20.0, 20.0, 15.0 };
	double factor3[3] = { 20.0, 20.0, 7.0 };

	raytracer.translate(sphere, Vector3D(5, 5, -7));

	raytracer.translate(cylinder, Vector3D(5, 10, -20));	
	//raytracer.rotate(cylinder, 'x', -60); 
	raytracer.rotate(cylinder, 'x', 90); 
	// raytracer.rotate(cylinder, 'y', 90);
	raytracer.scale(cylinder, Point3D(3, -1, 0), factor1);

	raytracer.translate(cone, Vector3D(7, 8, -16));	
	raytracer.rotate(cone, 'x', 25); 
	// raytracer.rotate(cone, 'y', 90);
	raytracer.rotate(cone, 'x', 90); 
	raytracer.scale(cone, Point3D(3, -1, 0), factor1);

	raytracer.translate(plane, Vector3D(7,10, -22));	
	//	raytracer.rotate(plane, 'x', 40); 
	//raytracer.rotate(plane, 'z', 90); 
	raytracer.rotate(plane, 'y', -45);
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
	raytracer.translate(plane2, Vector3D(-7, 10, -22));	
	//raytracer.rotate(plane2, 'x', -90); 
	//raytracer.rotate(plane2, 'z', -20); 
	raytracer.rotate(plane2, 'y', 45);
	raytracer.scale(plane2, Point3D(0, 0, 0), factor2);
	raytracer.translate(plane3, Vector3D(0, 0, -15));	
	raytracer.rotate(plane3, 'x', -90); 
	//raytracer.rotate(plane2, 'z', 45); 
	raytracer.rotate(plane3, 'z', -45);
	raytracer.scale(plane3, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	printf("Rendering image 1...");
	raytracer.render(width, height, eye, view, up, fov, "images/cc_view1.bmp");
	// Render it from a different point of view.
	Point3D eye2(4, 2, 5);
	Vector3D view2(-2, 0, -6);
	printf("Rendering image 2...");
	raytracer.render(width, height, eye2, view2, up, fov, "images/cc_view2.bmp");
}

void scene_single_textured_sphere(int width, int height){

	Raytracer raytracer;

	// Read bmp texture
	char texture_filename[] = "textures/checker.bmp";
	unsigned long int tex_w;
	long int tex_h;
	unsigned char *rarray, *garray, *barray;
	bool read_error = true;
	printf("Reading texture (%s)... ", texture_filename);
	read_error = bmp_read (texture_filename, &tex_w, &tex_h, &rarray, &garray, &barray);
	printf("DONE.\n");
	printf("Texture Width: %d, Height: %d\n", tex_w, tex_h);

	if (read_error) {
		printf("Texture read error, exiting...\n");
		return;
	}

	// Camera parameters.
	Point3D eye(0, 30, 1);
	Vector3D view(0, -17, -20);
	Vector3D up(0, 1, 0);
	double fov = 70;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
	Colour(0.628281, 0.555802, 0.366065), 
	51.2, 1.0, 0.0, 0.0);
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
	Colour(0.316228, 0.316228, 0.316228), 
	12.8, 1.0, 0.0, 0.0);
	Material randomCol( Colour(.6, .1, 0), Colour(0.12, 0.89, 0.9), 
	Colour(0.1, 0.9, 0.5), 
	51.2, 1.0, 0.0, 0.0);
	Material chrome( Colour(.25, .25, .25), Colour(0.4	,0.4,	0.4), 
	Colour(0.774597,	0.774597,	0.774597), 
	76.8, 1.0, 0.0, 0.0);
	Material texture( tex_w, tex_h, rarray, garray, barray );
	// gold.reflective=false;
	// jade.reflective=false;
	// chrome.reflective=false;
	// randomCol.reflective=true;

	// Defines a point light source.
	raytracer.addLightSource( new AreaLight(Point3D(0, 60, -15), 
	Vector3D(10,0,0), Vector3D(0,0,10), Colour(0.9, 0.9, 0.9) ) );
	// raytracer.get_light_list_node()->light->set_render(ren);
	// raytracer.get_light_list_node()->light->sample=1;
	// Add a unit square into the scene with material mat.
	// SceneDagNode* cylinder = raytracer.addObject( new Cylinder(), &chrome );
	// SceneDagNode* cone = raytracer.addObject( new Cone(), &gold );
	SceneDagNode* sphere = raytracer.addObject( new TexturedUnitSphere(), &texture);
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &jade );
	SceneDagNode* plane2 = raytracer.addObject( new UnitSquare(), &jade );
	SceneDagNode* plane3 = raytracer.addObject( new UnitSquare(), &jade );
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 3.0, 3.0, 8.0 };
	double factor2[3] = { 20.0, 20.0, 15.0 };
	double factor3[3] = { 7,7,7 };

	raytracer.translate(sphere, Vector3D(0, 10, -15));
	raytracer.scale(sphere, Point3D(0,0,0), factor3);

	// raytracer.translate(cylinder, Vector3D(5, 10, -15));	
	// //raytracer.rotate(cylinder, 'x', -60); 
	// raytracer.rotate(cylinder, 'x', 90); 
	// // raytracer.rotate(cylinder, 'y', 90);
	// raytracer.scale(cylinder, Point3D(3, -1, 0), factor1);

	// raytracer.translate(cone, Vector3D(7, 8, -16));	
	// raytracer.rotate(cone, 'x', 25); 
	// // raytracer.rotate(cone, 'y', 90);
	// raytracer.rotate(cone, 'x', 90); 
	// raytracer.scale(cone, Point3D(3, -1, 0), factor1);

	raytracer.translate(plane, Vector3D(7,10, -22));	
	//	raytracer.rotate(plane, 'x', 40); 
	//raytracer.rotate(plane, 'z', 90); 
	raytracer.rotate(plane, 'y', -45);
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
	raytracer.translate(plane2, Vector3D(-7, 10, -22));	
	//raytracer.rotate(plane2, 'x', -90); 
	//raytracer.rotate(plane2, 'z', -20); 
	raytracer.rotate(plane2, 'y', 45);
	raytracer.scale(plane2, Point3D(0, 0, 0), factor2);
	raytracer.translate(plane3, Vector3D(0, 0, -15));	
	raytracer.rotate(plane3, 'x', -90); 
	//raytracer.rotate(plane2, 'z', 45); 
	raytracer.rotate(plane3, 'z', -45);
	raytracer.scale(plane3, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	printf("Rendering image 1...");
	raytracer.render(width, height, eye, view, up, fov, "images/cc_view1.bmp");
	// Render it from a different point of view.
	// Point3D eye2(4, 2, 5);
	// Vector3D view2(-2, 0, -6);
	// printf("Rendering image 2...");
	// raytracer.render(width, height, eye2, view2, up, fov, "images/cc_view2.bmp");
}

void refraction_reflection_scene(int width, int height) { 
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it frm two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;

	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material blue( Colour(0, 0, 0), Colour(100.0/255, 149.0/255, 237.0/255), 
		Colour(0.31628, 0.31628, 0.31628), 
		12.8, 0.9, 1.3, 0.0);
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2, 1.0 , 0.0, 0.0);
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.0 , 0.0, 0.0);

	Material silver( Colour(0.2, 0.2, 0.2), Colour(0.50754, 0.50754, 0.50754),
			Colour(0.508273, 0.508273, 0.508273), 
			51.2, 1.0 , 0.0, 0.0);

	Material glass( Colour(0.05, 0.05, 0.05), Colour(0.05, 0.05, 0.05), 
			Colour(0.5, 0.5, 0.5),
			12.2, 0.0 , 1.0, 1.0); 


	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0,0,5), 
				Colour(0.9, 0.9, 0.9) ) );

	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere1 = raytracer.addObject( new UnitSphere(), &jade );
	SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* plane = raytracer.addObject(new UnitSquare(), &blue);
	SceneDagNode* sphere = raytracer.addObject( new Cone(), &glass );

	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 1.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	double factor3[3] = { 1.0, 1.0, 1.0 };

	raytracer.translate(sphere, Vector3D(0, 0, -2.0));
	raytracer.rotate(sphere, 'y', 180); 

	raytracer.translate(sphere1, Vector3D(-1, 0, -4));	
	raytracer.scale(sphere1, Point3D(0, 0, 0), factor3);

	raytracer.translate(sphere2, Vector3D(1, 0, -4));	
	raytracer.scale(sphere2, Point3D(0, 0, 0), factor3);


	raytracer.translate(plane, Vector3D(0, 0, -7));
	raytracer.rotate(plane, 'z', 45); 	
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
	
	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	printf("Rendering image 1...");
	raytracer.render(width, height, eye, view, up, fov, "images/rr_view1.bmp");
	
	// Render it from a different point of view.
	Point3D eye2(4, 2, 1);
	Vector3D view2(-4, -2, -6);
	printf("Rendering image 2...");
	raytracer.render(width, height, eye2, view2, up, fov, "images/rr_view2.bmp");

}
