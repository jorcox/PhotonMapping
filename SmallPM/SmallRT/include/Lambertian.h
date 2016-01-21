/*********************************************************************************
Copyright (C) 2014 Adrian Jarabo (ajarabo@unizar.es)
Copyright (C) 2014 Diego Gutierrez (diegog@unizar.es)
All rights reserved.

This is an educational Ray Tracer developed for the course 'Informatica Grafica'
(Computer Graphics) tought at Universidad de Zaragoza (Spain). As such, it does not 
intend to be fast or general, but just to provide an educational tool for undergraduate
students. 

This software is provided as is, and any express or implied warranties are disclaimed.
In no event shall copyright holders be liable for any damage.
**********************************************************************************/

#ifndef _PHONG_H_
#define _PHONG_H_

#include "BSDF.h"
#include "Vector3.h"
#include "Intersection.h"

class Lambertian: public BSDF
{
	//Diffuse coefficient
	Vector3 Kd;

	Vector3 brdf(const Vector3 &wi, const Vector3 &wo, const Vector3 &normal)const;
public:
	Lambertian(World *_w):Kd(.67),BSDF(_w){}
	Lambertian(World *_w, const Vector3 &_Kd):Kd(_Kd),BSDF(_w){}
	Vector3 shade(const Intersection &it)const;
	void get_outgoing_sample_ray(const Intersection &it, Ray &r, Real &pdf)const;
	Vector3 get_albedo(const Intersection &it) const;
	bool is_delta()const;
};

#endif