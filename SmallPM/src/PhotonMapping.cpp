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
#include "PhotonMapping.h"
#include "World.h"
#include "Intersection.h"
#include "Specular.h"
#include "Ray.h"
#include "BSDF.h"

const int MAX_PHOTONS = 100;

//*********************************************************************
// Compute the photons by tracing the Ray 'r' from the light source
// through the scene, and by storing the intersections with matter
// in the lists 'xx_photons', storing the diffuse (global) and caustic
// photons respectively. For efficiency, both are computed at the same
// time, since computing them separately would result into a lost of
// several samples marked as caustic or diffuse.
// Same goes with the boolean 'direct', that specifies if direct 
// photons (from light to surface) are being stored or not. 
// The initial traced photon has energy defined by the tristimulus
// 'p', that accounts for the emitted power of the light source.
// The function will return true when there are more photons (caustic
// or diffuse) to be shot, and false otherwise.
//---------------------------------------------------------------------
bool PhotonMapping::trace_ray(const Ray& r, const Vector3 &p,
	std::list<Photon> &global_photons, std::list<Photon> &caustic_photons, bool direct)
{

	//Check if max number of shots done...
	if (++m_nb_current_shots > m_max_nb_shots)
	{
		return false;
	}

	// Compute irradiance photon's energy
	Vector3 energy(p);

	Ray photon_ray(r);
	photon_ray.shift();

	bool is_caustic_particle = false;

	//Iterate the path
	while (1)
	{
		// Throw ray and update current_it
		Intersection it;
		world->first_intersection(photon_ray, it);

		if (!it.did_hit())
			break;

		//Check if has hit a delta material...
		if (it.intersected()->material()->is_delta())
		{
			// If delta material, then is caustic...
			// Don't store the photon!
			is_caustic_particle = true;
		}
		else if (photon_ray.get_level() > 0 || direct)
		{
			//If non-delta material, store the photon!
			if (is_caustic_particle)
			{
				//If caustic particle, store in caustics
				if (caustic_photons.size() < m_nb_caustic_photons)
					caustic_photons.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy));
			}
			else
			{
				//If non-caustic particle, store in global
				if (global_photons.size() < m_nb_global_photons)
					global_photons.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy));
			}
			is_caustic_particle = false;
		}

		Real pdf;

		Vector3 surf_albedo = it.intersected()->material()->get_albedo(it);
		Real avg_surf_albedo = surf_albedo.avg();

		Real epsilon2 = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);
		while (epsilon2 < 0.)
			epsilon2 = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);

		if (epsilon2 > avg_surf_albedo || photon_ray.get_level() > 20)
			break;

		// Random walk's next step
		// Get sampled direction plus pdf, and update attenuation
		it.intersected()->material()->get_outgoing_sample_ray(it, photon_ray, pdf);

		// Shade...
		energy = energy*surf_albedo;
		if (!it.intersected()->material()->is_delta())
			energy *= dot_abs(it.get_normal(), photon_ray.get_direction()) / 3.14159;

		energy = energy / (pdf*avg_surf_albedo);
	}

	if (caustic_photons.size() == m_nb_caustic_photons &&
		global_photons.size() == m_nb_global_photons)
	{
		m_max_nb_shots = m_nb_current_shots - 1;
		return false;
	}

	return true;
}

//*********************************************************************
// TODO: Implement the preprocess step of photon mapping,
// where the photons are traced through the scene. To do it,
// you need to follow these steps for each shoot:
//  1 - Sample a world's light source in the scene to create
//		the initial direct photon from the light source.
//	2 - Trace the photon through the scene storing the inter-
//		sections between the photons and matter. You can use
//		the function 'trace_ray' for this purpose.
//	3 - Finally, once all the photons have been shot, you'll
//		need to build the photon maps, that will be used later
//		for rendering. 
//		NOTE: Careful with function
//---------------------------------------------------------------------
void PhotonMapping::preprocess()
{

	std::vector<LightSource*> ls = world->light_source_list;
	/*
	Se generan puntos aleatorios dentro de un cubo. Solo se admiten para lanzar fotones desde este punto aquellos
	que estén  dentro de una esfera, la cual esté dentro de ese cubo y su diametro sea igual al lado del cubo.
	Esto se realiza debido a que sino, se concentrarían en las esquinas.
	*/
	Vector3 centro(0, 0, 0);
	int radio = 1;
	Real fmin = -(radio);
	std::list<Photon> *global_photons = new std::list<Photon>();
	std::list<Photon> *caustic_photons = new std::list<Photon>();
	Ray *r;
	/*
	   For each light
	 */
	for (int i = 0; i < ls.size(); i++)
	{
		m_nb_current_shots = 0;
		LightSource *luz = ls[i];
		Vector3 flujo((luz->get_intensities().operator/(m_max_nb_shots/3.5)));
		do{
			

			Real x;
			Real y;
			Real z;
			Real dist = radio + 1;
			// Iterate until it finds a point inside the sphere
			while (dist > radio){
				// Generates points between (-radius,radius)
				x = fmin + ((double)rand() / RAND_MAX) * radio*2;
				y = fmin + ((double)rand() / RAND_MAX) * radio*2;
				z = fmin + ((double)rand() / RAND_MAX) * radio*2;
				dist = sqrt(x*x + y*y + z*z);
			}
			Vector3 random = Vector3(x, y, z);
			r = new Ray(luz->get_position(), random.normalize(), 0);

		// Call to trace_ray until it's not neccesary to shot more photons
		} while (trace_ray(*r, flujo, *global_photons, *caustic_photons, false));

		unsigned int n_phot_g = m_max_nb_shots;

		
		/*
		  Bucle de guardado de fotones globales
		*/
		int n_global = global_photons->size();
		for (int j = 0; j < n_global; j++){
			Photon photon = global_photons->front();
			global_photons->pop_front();
			m_global_map.store(std::vector<Real>(photon.position.data,
				photon.position.data + 3), photon);
		}
				
		/*
		  Bucle de guardado de fotones causticos
		*/
		int n_caustics = caustic_photons->size();
		for (int j = 0; j < n_caustics; j++){
			Photon photon = caustic_photons->front();
			caustic_photons->pop_front();
			m_caustics_map.store(std::vector<Real>(photon.position.data,
				photon.position.data + 3), photon);
		}
	}

	m_global_map.balance();
	if (!m_caustics_map.isEmptyPreBalanced()){
		m_caustics_map.balance();
	}


}


Real PhotonMapping::distancia(Vector3* a, Vector3 b){
	return sqrt(pow(a->getComponent(1) - b.getComponent(1), 2) +
		pow(a->getComponent(2) - b.getComponent(2), 2) +
		pow(a->getComponent(3) - b.getComponent(3), 2)
		);
}

//*********************************************************************
// TODO: Implement the function that computes the rendering equation 
// using radiance estimation with photon mapping, using the photon
// maps computed as a proprocess. Note that you will need to handle
// both direct and global illumination, together with recursive the 
// recursive evaluation of delta materials. For an optimal implemen-
// tation you should be able to do it iteratively.
// In principle, the class is prepared to perform radiance estimation
// using k-nearest neighbors ('m_nb_photons') to define the bandwidth
// of the kernel.
//---------------------------------------------------------------------
Vector3 PhotonMapping::shade(Intersection &it0)const
{
	//Vector3 L(0);
	Intersection it(it0);
	int NUM_PHOTONS = this->m_nb_photons;
	bool repeat = false;

	Vector3 kd = it.intersected()->material()->get_albedo(it);

	Vector3 direct;

	/*********************************************************************************************************/
	/*                                           DIRECT/LOCAL                                                */
	/*********************************************************************************************************/
	int j = 0;
	
	do{

		repeat = false;
		/* Current kd */
		kd = it.intersected()->material()->get_albedo(it);

		Vector3 normal = normalize(it.get_normal());

		std::vector<LightSource*> worldLights = world->light_source_list;
		Vector3 difuseIntensity;
		for (size_t i = 0; i < worldLights.size(); i++)
		{
			LightSource* ls = worldLights[i]; // Current light
			if (world->light(i).is_visible(it.get_position()) || it.intersected()->material()->is_delta()){
				if (it.intersected()->material()->is_delta()){
					Ray newRay;
					Real s;
					it.intersected()->material()->get_outgoing_sample_ray(it, newRay, s);
					newRay.shift();
					world->first_intersection(newRay, it);
					repeat = true;
					j++;
				}
				else {

					/* Difuse */
					Vector3 light = normalize(ls->get_position() - it.get_position());
					double difAngle = dot(normal, light);

					/* If angle < 0 ==> Not illuminated ==> angle = 0 */
					if (difAngle < 0) {
						difAngle = 0;
					}
					difuseIntensity +=  ls->get_incoming_light(it.get_position()).operator*(difAngle);
				}
			}
		}
		direct = kd * difuseIntensity;
	} while (repeat && j<10);

	/*********************************************************************************************************/
	/*                                              GLOBAL                                                   */	
	/*********************************************************************************************************/

	/* Photon mapping algorithm for Global Illumination */
	std::vector<const KDTree<Photon, 3>::Node*> global_photons;
	std::vector<const KDTree<Photon, 3>::Node*> caustic_photons;

	Real maxDistanceGlobal;
	Real maxDistanceCaustic;

	/* Find nearest global photons */
	m_global_map.find(std::vector<Real>(it.get_position().data, it.get_position().data + 3),
		NUM_PHOTONS, global_photons, maxDistanceGlobal);

	/* Find nearest caustic photons */
	m_caustics_map.find(std::vector<Real>(it.get_position().data, it.get_position().data + 3),
		NUM_PHOTONS, caustic_photons, maxDistanceCaustic);



	

	/* Cone area */

	double k = 1.3;
	Real coneAreaGlobal = 1 / ((1 - (2 / (3 * k)))*3.141592*maxDistanceGlobal*maxDistanceGlobal);
	Real coneAreaCaustic = 1 / ((1 - (2 / (3 * k)))*3.141592*maxDistanceCaustic*maxDistanceCaustic);
	/*
	Real sphereAreaGlobal = 1 / ((1 - (2 / (3 * k)))*3.141592*maxDistanceGlobal*maxDistanceGlobal);
	Real sphereAreaCaustic = 1 / ((1 - (2 / (3 * k)))*3.141592*maxDistanceCaustic*maxDistanceCaustic);
	*/
	Vector3 accumulatedGlobal = Vector3(0, 0, 0);
	Vector3 accumulatedCaustic = Vector3(0, 0, 0);

	/*
	*  Radiance estimation for global photons
	*/
	for (int i = 0; i < global_photons.size(); i++)	{
		// Area
		// Flux

		const KDTree<Photon, 3>::Node* nodePhoton = global_photons.at(i);
		Photon photon = nodePhoton->data();
		Vector3 flux = photon.flux;

		// Fr(Wp,Wr) -> BRDF
		//Vector3 vphot = photon.direction;
		Vector3 aux = photon.flux.operator*(kd); 


		// Cone
		Vector3 restaPosition = photon.position.operator-(it.get_position());
		float distance = restaPosition.length2();
		aux = aux.operator*(1 - (distance) / (k * maxDistanceGlobal));

		accumulatedGlobal.operator += (aux);
	}

	accumulatedGlobal.operator *= (coneAreaGlobal);

	

	/*
	*  Radiance estimation for caustic photons
	*/
	for (int i = 0; i < caustic_photons.size(); i++) {
		// Area
		// Flux

		const KDTree<Photon, 3>::Node* nodePhoton = caustic_photons.at(i);
		Photon photon = nodePhoton->data();
		Vector3 flux = photon.flux;

		// Fr(Wp,Wr) -> BRDF
		Vector3 vphot = photon.direction;
		Vector3 aux = photon.flux.operator*(kd);


		// Cone
		Vector3 restaPosition = photon.position.operator-(it.get_position());
		float distance = restaPosition.length2();
		accumulatedCaustic.operator += (aux.operator*(1 - (distance) / (k*maxDistanceCaustic)));
	}

	accumulatedCaustic.operator *= (coneAreaCaustic);


	/*
	 * Final summation
	 */
	Vector3 in = accumulatedCaustic.operator+(direct).operator+(accumulatedGlobal);
	//Vector3 in = accumulatedCaustic.operator+(accumulatedGlobal);

	return in;

	/*
	// Photon mapping algorithm for Caustics
	std::vector<const KDTree<Photon, 3>::Node*> caustic_photons;
	max_distance = 1;
	do {
	m_caustics_map.balance();

	max_distance++;
	} while (caustic_photons.size() < MAX_PHOTONS);
	*/

	/*
	//**********************************************************************
	// The following piece of code is included here for two reasons: first
	// it works as a 'hello world' code to check that everthing compiles
	// just fine, and second, to illustrate some of the functions that you
	// will need when doing the work. Goes without saying: remove the
	// pieces of code that you won't be using.
	//
	unsigned int debug_mode = 1;

	switch (debug_mode)
	{
	case 1:
	// ----------------------------------------------------------------
	// Display Albedo Only
	L = it.intersected()->material()->get_albedo(it);
	break;
	case 2:
	// ----------------------------------------------------------------
	// Display Normal Buffer
	L = it.get_normal();
	break;
	case 3:
	// ----------------------------------------------------------------
	// Display whether the material is specular (or refractive)
	L = Vector3(it.intersected()->material()->is_delta());
	break;

	case 4:
	// ----------------------------------------------------------------
	// Display incoming illumination from light(0)
	L = world->light(0).get_incoming_light(it.get_position());
	break;

	case 5:
	// ----------------------------------------------------------------
	// Display incoming direction from light(0)
	L = world->light(0).get_incoming_direction(it.get_position());
	break;

	case 6:
	// ----------------------------------------------------------------
	// Check Visibility from light(0)
	if (world->light(0).is_visible(it.get_position()))
	L = Vector3(1.);
	break;
	}
	// End of exampled code
	//**********************************************************************

	return L;
	*/
}