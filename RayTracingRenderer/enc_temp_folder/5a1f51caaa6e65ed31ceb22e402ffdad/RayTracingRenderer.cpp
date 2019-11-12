// RayTracingRenderer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "constant.h"
#include "utility.h"

#include <iostream>
#include <scene.h>
#include <xmlload.h>
#include <viewport.h>
#include <objects.h>
#include <omp.h>
#include <cyBVH.h>
#include <time.h>
#include <math.h>

Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
Plane thePlane;
ItemFileList<Object> objList;
LightList lights;
MaterialList materials;
std::vector<NodeMtl> nodeMtlList;
TexturedColor background;
TexturedColor environment;
TextureList textureList;

Vec3f pixelx;
Vec3f pixely;

using namespace Utility;

int main()
{
	//LoadScene(".\\xmlfiles\\playground.xml");
	//LoadScene(".\\xmlfiles\\catscene.xml");
	//LoadScene(".\\xmlfiles\\potscene.xml");
	//LoadScene(".\\xmlfiles\\assignment9.xml");
	LoadScene(".\\xmlfiles\\assignment11.xml");
	//LoadScene(".\\xmlfiles\\assignment11.xml");
	ShowViewport();
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

Color RayTraversing(Node * node, Ray ray, HitInfo & hit, int bounce)
{
	Node * currentnode;
	int numberofchild = node->GetNumChild();
	HitInfo tmphit = HitInfo();

	for (int i = 0; i < numberofchild; i++)
	{
		currentnode = node->GetChild(i);
		Ray currentray = currentnode->ToNodeCoords(ray);

		if (currentnode->GetNodeObj() != nullptr)
		{
			if (currentnode->GetNodeObj()->IntersectRay(currentray, tmphit, 0))
			{
				tmphit.node = currentnode;
				currentnode->FromNodeCoords(tmphit);
				hit = tmphit;
			}
		}

		if (currentnode != nullptr)
		{
			Node childnode;
			RayTraversing(currentnode, currentray, hit, bounce);
			if (hit.node != nullptr && hit.node != tmphit.node)
			{
				currentnode->FromNodeCoords(hit);
				tmphit = hit;
			}
		}
	}

	//Shading
	if (node == &rootNode)
	{
		hit = tmphit;
		if (hit.node != nullptr)
		{
			return materials.Find(hit.node->GetMaterial()->GetName())->Shade(ray, hit, lights, bounce);
		}

		if(bounce != TIMEOFREFRECTION)
			return environment.SampleEnvironment(ray.dir);
	}
	return Color(0, 0, 0);
}

Color BlurEffect(Ray ray, Node * traversingnode, Node * node, float & zbuffer)
{
	Ray cameraraies[RAYPERPIXELFORBLUREFFECT];
	HitInfo hits[RAYPERPIXELFORBLUREFFECT];

	Color averagepixelcolor = Color(0, 0, 0);
	Color pixelcolors[RAYPERPIXELFORBLUREFFECT];
	Vec3f screenpoints[RAYPERPIXELFORBLUREFFECT];

	Color answercolor = Color(0, 0, 0);

	for (int i = 0; i < RAYPERPIXELFORBLUREFFECT; i++)
	{
		screenpoints[i] = ray.p + camera.focaldist * ray.dir;

		Vec3f offset = Vec3f(0, 0, 0);
		CircleUniformSampling(offset, camera.dir, camera.up, camera.dof);
		cameraraies[i].p = ray.p + offset;

		cameraraies[i].dir = screenpoints[i] - cameraraies[i].p;
		cameraraies[i].Normalize();

		pixelcolors[i] = RayTraversing(traversingnode, cameraraies[i], hits[i], TIMEOFREFRECTION);
		answercolor += pixelcolors[i];
	}

	return answercolor / RAYPERPIXELFORBLUREFFECT;
}

Color AdaptiveSampling(Ray ray, float length, uint8_t & samplecount)
{
	samplecount++;
	Ray cameraraies[RAYPERSAMPLING];
	HitInfo hits[RAYPERSAMPLING];

	Color averagepixelcolor = Color(0, 0, 0);
	Color pixelcolors[RAYPERSAMPLING];
	Vec3f screenpoints[RAYPERSAMPLING];

	for (int i = 0; i < RAYPERSAMPLING; i++)
	{
		if (i % RAYPERSAMPLING == 0)
		{
			screenpoints[i] = ray.p + camera.focaldist * ray.dir;
		}
		else if (i % RAYPERSAMPLING == 1)
		{
			screenpoints[i] = ray.p + camera.focaldist * ray.dir - length * pixelx;
		}
		else if (i % RAYPERSAMPLING == 2)
		{
			screenpoints[i] = ray.p + camera.focaldist * ray.dir - length * pixely;
		}
		else if (i % RAYPERSAMPLING == 3)
		{
			screenpoints[i] = ray.p + camera.focaldist * ray.dir - length * pixelx - length * pixely;
		}

		Vec3f offset = Vec3f(0, 0, 0);
		SquareUniformSampling(offset, pixelx, pixely, length/2);
		screenpoints[i] += offset;

		cameraraies[i].p = ray.p;
		cameraraies[i].dir = screenpoints[i] - cameraraies[i].p;
		cameraraies[i].Normalize();

		pixelcolors[i] = RayTraversing(&rootNode, cameraraies[i], hits[i], TIMEOFREFRECTION);
		averagepixelcolor += pixelcolors[i];
	}

	averagepixelcolor /= RAYPERSAMPLING;
	Color variance = Color(0, 0, 0);
	
	for (int i = 0; i < RAYPERSAMPLING; i++)
	{
		variance += ((pixelcolors[i] - averagepixelcolor) * (pixelcolors[i] - averagepixelcolor));
	}
	variance /= RAYPERSAMPLING;

	Color returncolor = Color(0, 0, 0);

	for (int i = 0; i < RAYPERSAMPLING; i++)
	{
		if (variance.r >= SAMPLEVARIENCE || variance.g >= SAMPLEVARIENCE || variance.b >= SAMPLEVARIENCE)
		{
			if (samplecount < MAXSAMPLECOUNT)
			{
				returncolor += AdaptiveSampling(cameraraies[i], length / 2.0f, samplecount);
			}
			else
			{
				returncolor += pixelcolors[i];
			}
		}
		else
		{
			returncolor += pixelcolors[i];
		}
	}

	return returncolor / RAYPERSAMPLING;
}

void BeginRender() {

	time_t time0;   // create timers.
	time_t time1;

	time(&time0);   // get current time.


	Node * startnode = &rootNode;

	float * zbuffers = renderImage.GetZBuffer();
	uint8_t * samplecount = renderImage.GetSampleCount();
	Color24* pixels = renderImage.GetPixels();
	Ray * cameraray = new Ray[renderImage.GetHeight() * renderImage.GetWidth()];

	float l = camera.focaldist;
	float h = 2 * l * tanf((camera.fov / 2) * static_cast<float>(M_PI) / 180);
	float w = camera.imgWidth * (h / camera.imgHeight);

	int H = camera.imgHeight;
	int W = camera.imgWidth;

	Vec3f x = camera.dir.Cross(camera.up);
	Vec3f y = camera.up;

	Vec3f f = camera.pos + l * camera.dir + (h / 2) * y - (w / 2) * x;

	pixelx = (w / W)*x;
	pixely = (h / H)*y;

#pragma omp parallel for
	for (int i = 0; i < renderImage.GetHeight(); i++) {
		for (int j = 0; j < renderImage.GetWidth(); j++) {

			//hit.duvw[0] = (w / W) * x;
			//hit.duvw[1] = (h / H) * y;

			cameraray[i * renderImage.GetWidth() + j].dir = f + (j + HALF) * (w / W)*x - (i + HALF) * (h / H)*y - camera.pos;
			cameraray[i * renderImage.GetWidth() + j].p = camera.pos;
			zbuffers[i * renderImage.GetWidth() + j] = BIGFLOAT;
			samplecount[i * renderImage.GetWidth() + j] = 0;
			cameraray[i * renderImage.GetWidth() + j].Normalize();

			Color resultColor;

		#ifdef ANTIALIASING
			resultColor = AdaptiveSampling(cameraray[i * renderImage.GetWidth() + j], HALF, samplecount[i * renderImage.GetWidth() + j]);
		#elif defined BLUREFFECT
			resultColor = BlurEffect(cameraray[i * renderImage.GetWidth() + j], startnode, node, zbuffers[i * renderImage.GetWidth() + j]);
		
		#elif defined NOANTIALIASING
			HitInfo hit = HitInfo();
			//if (i == 271 && j == 300)
			//{
			//	pixels[i * renderImage.GetWidth() + j] = (Color24)OutRayTraversing(startnode, cameraray[i * renderImage.GetWidth() + j], zbuffers[i * renderImage.GetWidth() + j], hit);
			//}
			resultColor = RayTraversing(startnode, cameraray[i * renderImage.GetWidth() + j], hit, TIMEOFREFRECTION);
			zbuffers[i * renderImage.GetWidth() + j] = hit.z;
		#endif 

			resultColor.r = pow(resultColor.r, 1 / 2.2f);
			resultColor.g = pow(resultColor.g, 1 / 2.2f);
			resultColor.b = pow(resultColor.b, 1 / 2.2f);

			pixels[i * renderImage.GetWidth() + j] = (Color24)resultColor;

			if (pixels[i * renderImage.GetWidth() + j] == (Color24)Color(0, 0, 0))
			{
				Vec3f v((float)j / renderImage.GetWidth(), (float)i / renderImage.GetHeight(), 0.0f);
				pixels[i * renderImage.GetWidth() + j] = (Color24)background.Sample(v);
				//pixels[i * renderImage.GetWidth() + j] = Color24(0, 255, 255);
			}
		}
	}

	time(&time1);
	double seconds = (double)(time1 - time0);

	printf("Done. Time was %f \n", seconds);
	renderImage.SaveImage("saveimage.png");
	delete cameraray;

	return;
}

void StopRender() {
}

/* Unfinished function, do not used this*/
//Color InnerRayTraversing(Node * node, Ray ray, HitInfo & hit, int bounce)
//{
//	Node * currentnode = node;
//	Node * previousnode;
//	Node * parentnode = node;
//
//	std::vector<Node *> parentsnodes;
//
//	Ray currentray = ray;
//	Ray previousray;
//	Ray parentray = ray;
//
//	int numberofchildren = currentnode->GetNumChild();
//	int i = 0;
//	int j = 0;
//
//	while (i < numberofchildren)
//	{
//		previousnode = currentnode;
//		previousray = currentray;
//
//		j = 0;
//		numberofchildren = previousnode->GetNumChild();
//
//		while (j < numberofchildren)
//		{
//			currentnode = previousnode;
//			currentray = previousray;
//
//			currentnode = currentnode->GetChild(j);
//			currentray = currentnode->ToNodeCoords(currentray);
//
//			if (currentnode->GetNodeObj() != nullptr)
//			{
//				if (currentnode->GetNodeObj()->IntersectRay(currentray, hit, 0))
//				{
//					hit.currentnode = currentnode;
//					currentnode->FromNodeCoords(hit);
//				}
//			}
//			j++;
//		}
//
//		if (i == parentnode->GetNumChild())
//		{
//			parentsnodes.push_back(parentnode);
//
//			i = 0;
//			parentnode = parentnode->GetChild(i);
//			parentray = parentnode->ToNodeCoords(parentray);
//
//			numberofchildren = parentnode->GetNumChild();
//			if (numberofchildren == 0)
//			{
//				break;
//			}
//		}
//		
//		currentnode = parentnode->GetChild(i);
//		currentray = currentnode->ToNodeCoords(parentray);
//		i++;
//	}
//
//	//Shading
//	if (currentnode->GetNodeObj() != nullptr)
//	{
//		if (hit.currentnode != nullptr)
//		{
//			if (materials.Find(currentnode->GetMaterial()->GetName()) != nullptr)
//			{
//				return materials.Find(hit.currentnode->GetMaterial()->GetName())->Shade(ray, hit, lights, bounce);
//			}
//		}
//		//return background.Sample(originalray.dir);
//	}
//	return environment.SampleEnvironment(ray.dir);
//}

float FresnelReflections(const HitInfo & hInfo, float refractionIndex, float cos1)
{
	float R0;
	if (hInfo.front)
	{
		R0 =      ((1 - refractionIndex) * (1 - refractionIndex))
			  / //----------------------------------------------
			      ((1 + refractionIndex) * (1 + refractionIndex));
	}
	else
	{
		R0 =    ((refractionIndex - 1) * (refractionIndex - 1))
			 / //-----------------------------------------------
			    ((1 + refractionIndex) * (1 + refractionIndex));
	}

	return R0 + (1 - R0) * pow(1 - cos1, 5);
}

Color Reflection(Ray const & ray, const HitInfo & hInfo, int bounce, const float glossiness)
{
	// bounce 1 is reflect 1 time
	if (bounce <= 0)
	{
		return Color(0, 0, 0);
	}
	else
	{
		Vec3f N = hInfo.N;

		//if(glossiness !=0)
		ConeUniformSampling(N, glossiness);

		Vec3f P = N;
		P = SHADOWBIAS * P; 
		P += hInfo.p;

		Vec3f V = -1 * ray.dir;
		Vec3f R = 2 * (N.Dot(V)) * N - V;
		R.Normalize();

		// S is a starting point from the surface point
		Ray S;
		S.dir = R;
		S.p = P;

		// Start traversing currentnode 
		Node node;
		Node * startnode = &rootNode;
		HitInfo hitinfo;

		Color returnColor = RayTraversing(startnode, S, hitinfo, bounce - 1);
		return returnColor;
	}
}

Color TotalInternalReflection(Ray ray, int bounce)
{
	// Start traversing currentnode 
	Node node;
	Node * startnode = &rootNode;
	HitInfo hitinfo;
	Color returnColor = Color(0, 0, 0);

	returnColor = RayTraversing(startnode, ray, hitinfo, bounce-1);

	return returnColor;
}

Color Refraction(Ray const & ray, const HitInfo & hInfo, int bounce, float refractionIndex, Color refraction, const float glossiness)
{
	float R;
	// Add bounce + 1 because in most case refraction need to go through front and back sides of the object
	if (bounce <= -1)
	{
		return Color(0, 0, 0);
	}
	else
	{
		Vec3f V = -1 * ray.dir;
		V.Normalize();

		Vec3f N = hInfo.N;
		ConeUniformSampling(N, glossiness);

		Vec3f P;

		float cos1, cos2, sin1, sin2;

		if (V.Dot(N) >= 0)
		{
			P = -1 * SHADOWBIAS * N;  P += hInfo.p;
			cos1 = V.Dot(N);
		}
		else
		{
			P = SHADOWBIAS * N;  P += hInfo.p;
			cos1 = V.Dot(-N);
		}

		sin1 = sqrt(1 - (cos1 * cos1));

		

		// Calculate fresnel reflection
		R = FresnelReflections(hInfo, refractionIndex, cos1);

		Vec3f T_h, T_v, T;

		if (V.Dot(N) >= 0)
		{
			sin2 = (1 / refractionIndex) * sin1;
		}
		else
		{
			sin2 = refractionIndex * sin1;
		}

		// S is a starting point from the surface point
		Ray S;

		// Total internal reflection
		if (sin2 > 1)
		{
			T = ray.dir - 2 * (ray.dir.Dot(N)) * N;
			P -= 2 * SHADOWBIAS * N;
			S.dir = T;
			S.p = P;
			S.Normalize();
			return TotalInternalReflection(S, bounce);
		}
		else // Normal procedure
		{
			cos2 = sqrt(1 - (sin2 * sin2));

			if (V.Dot(N) >= 0)
			{
				// Horizontal dirction Vector
				T_h = -cos2 * N;
				// Vertical Direction Vector
				T_v = (V - (V.Dot(N))* N);
			}
			else
			{
				// Horizontal direction Vector
				T_h = -cos2 * -N;
				// Vertical Direction Vector
				T_v = (V - (V.Dot(-N))* -N);
			}

			T_v.Normalize();
			T_v = -sin2 * T_v;
			// Combined horizontal and vertical
			T = T_h + T_v;

			S.dir = T;
			S.p = P;

			// Start traversing currentnode 
			Node node;
			Node * startnode = &rootNode;
			HitInfo hitinfo;

			Color returnColor = Color(0, 0, 0);

			// Refraction part
			returnColor = (1 - R) * refraction * RayTraversing(startnode, S, hitinfo, bounce - 1);
			// Reflection part
			returnColor += R * refraction * Reflection(ray, hInfo, bounce, glossiness);
			return returnColor;
		}
	}
}

Color CalculateAbsorption(Color in, Color absorption, float distance)
{
	float rout = exp(-1 * absorption.r * distance) * in.r;
	float gout = exp(-1 * absorption.g * distance) * in.g;
	float bout = exp(-1 * absorption.b * distance) * in.b;
	return Color(rout, gout, bout);
}

Color MtlBlinn::Shade(Ray const & ray, const HitInfo & hInfo, const LightList & lights, int bounce) const
{
	Vec3f N = hInfo.N;
	Color color = Color();

	Color specularpart;
	Color diffusepart;

	for (auto light = lights.begin(); light != lights.end(); ++light)
	{
		if ((*light)->IsAmbient())
		{
			color += (*light)->Illuminate(hInfo.p, N) * this->diffuse.Sample(hInfo.uvw, hInfo.duvw);
		}
		else if (strcmp((*light)->GetName(), "directLight") == 0)
		{
			Vec3f L = -1 * (*light)->Direction(hInfo.p);
			L.Normalize();
			Vec3f V = ray.p - hInfo.p;
			V.Normalize();
			Vec3f H = (V + L) / (V + L).Length();
			H.Normalize();

			// Incoming light
			Color IR = Color(0, 0, 0);
			if (N.Dot(-1 * (*light)->Direction(hInfo.p)) >= 0)
			{
				IR = (*light)->Illuminate(hInfo.p, N) * N.Dot(L);
			}

			float oneofcos = 1 / hInfo.N.Dot(L);
			Color specularpart = oneofcos * pow(H.Dot(hInfo.N), this->glossiness) * this->specular.Sample(hInfo.uvw, hInfo.duvw);
			Color diffusepart = this->diffuse.Sample(hInfo.uvw, hInfo.duvw);
			color += (diffusepart + specularpart) * IR;
		}
		else if (strcmp((*light)->GetName(), "pointLight") == 0)
		{
			Vec3f L = -1 * (*light)->Direction(hInfo.p);
			L.Normalize();
			Vec3f V = ray.p - hInfo.p;
			V.Normalize();
			Vec3f H = (V + L) / (V + L).Length();
			H.Normalize();

			// Incoming light
			Color IR = Color(0, 0, 0);
			if (N.Dot(-1 * (*light)->Direction(hInfo.p)) >= 0)
			{
				IR = (*light)->Illuminate(hInfo.p, N) * N.Dot(L);
			}

			float oneofcos = 1 / hInfo.N.Dot(L);
			specularpart = oneofcos * pow(H.Dot(hInfo.N), this->glossiness) * this->specular.Sample(hInfo.uvw, hInfo.duvw);
			diffusepart = this->diffuse.Sample(hInfo.uvw, hInfo.duvw);
			color += (diffusepart + specularpart) * IR;
		}
	}

#ifdef MONTECARLOGI
	// For global illumination
	if (bounce > 0)
	{
		Color returnColor = Color(0, 0, 0);
		Vec3f N_dash = N;

		for (int i = 0; i < MONTECARLOGI; i++)
		{
			CosineWeightedHemisphereUniformSampling(N_dash);

			// Start traversing currentnode 
			Node node;
			Node * startnode = &rootNode;
			HitInfo hitinfo;

			Ray ray_gi;
			ray_gi.p = hInfo.p;
			ray_gi.p += SHADOWBIAS * N;
			ray_gi.dir = N_dash;

			returnColor += this->diffuse.Sample(hInfo.uvw, hInfo.duvw) * RayTraversing(startnode, ray_gi, hitinfo, bounce - 2);
		}
		returnColor /= MONTECARLOGI;
		color += returnColor;
	}
#endif

	// Caluculate only relfection part for reflection
	if (this->reflection.Sample(hInfo.uvw) != Color(0, 0, 0))
	{
		Color reflectioncolor = Color(0, 0, 0);
		for (int i = 0; i < RAYPERPIXELFORGLOSSINESS; i++)
		{
			reflectioncolor += this->reflection.Sample(hInfo.uvw) * Reflection(ray, hInfo, bounce, reflectionGlossiness);
		}
		reflectioncolor /= RAYPERPIXELFORGLOSSINESS;

		color += reflectioncolor;
	}

	// Caluculate refraction part
	if (this->refraction.Sample(hInfo.uvw) != Color(0, 0, 0))
	{
		// When it is a back side hit, it means that absorption gonna happen during inside the material the light go through
		if (!hInfo.front)
		{
			color += CalculateAbsorption(color, absorption, hInfo.z);
		}

		Color refractioncolor = Color(0, 0, 0);
		for (int i = 0; i < RAYPERPIXELFORGLOSSINESS; i++)
		{
			refractioncolor += Refraction(ray, hInfo, bounce, ior, refraction.Sample(hInfo.uvw), refractionGlossiness);
		}
		refractioncolor /= RAYPERPIXELFORGLOSSINESS;

		color += refractioncolor;
	}

	return color;
}

bool DetectShadow(Node * traversingnode, Node * node, Ray ray, float t_max)
{
	int numberofchild = traversingnode->GetNumChild();
	Ray currentray = ray;
	for (int i = 0; i < numberofchild; i++)
	{
		node = traversingnode->GetChild(i);
		if (node->GetNodeObj() != nullptr)
		{
			Ray changedray = node->ToNodeCoords(currentray);
			HitInfo fake; fake.z = t_max;
			if (node->GetNodeObj()->IntersectRay(changedray, fake, 1))
				return true;
		}

		if (node != nullptr)
		{
			Node childnode;
			if (DetectShadow(node, &childnode, node->ToNodeCoords(currentray), t_max))
			{
				return true;
			}
		}
	}
	return false;
}

bool Sphere::IntersectRay(Ray const & ray, HitInfo & hInfo, int hitSide) const
{
	float a = ray.dir.Dot(ray.dir);
	float b = 2 * ray.dir.Dot(ray.p);
	float c = ray.p.Dot(ray.p) - 1;

	if (b*b - 4 * a*c >= 0)
	{

		float answer1 = (-1 * b + sqrt(b*b - 4 * a*c)) / (2 * a);
		float answer2 = (-1 * b - sqrt(b*b - 4 * a*c)) / (2 * a);

		float large;
		float small;

		float answer;

		if (answer1 >= answer2)
		{
			large = answer1; small = answer2;
		}
		else
		{
			large = answer2; small = answer1;
		}

		if (small < 0)
		{
			if (large > 0)
			{
				// this is for shadowprocess
				if (hitSide == 1)
				{
					answer = large;
					if ((answer > SHADOWBIAS) & (answer <= hInfo.z))
					{
						return true;
					}
				}
				else
				{
					if(hInfo.z > large ? true : false)
					{
						hInfo.z = large;
						hInfo.front = false;
						hInfo.p = ray.p + large * ray.dir;
						hInfo.N = hInfo.p;
						float u = (1 / 2 * static_cast<float>(M_PI)) * atan2f(hInfo.p.y, hInfo.p.x) + .5f;
						float v = (1 / static_cast<float>(M_PI)) * asinf(hInfo.p.z) + 0.5f;
						hInfo.uvw = Vec3f(u, v, 0.0f);
						return true;
					}
				}
			}
		}
		else
		{
			// this is for shadowprocess
			if (hitSide == 1)
			{
				answer = small;
				if ((answer > SHADOWBIAS) & (answer <= hInfo.z))
				{
					return true;
				}
			}
			else
			{
				if (hInfo.z > small ? true : false)
				{
					hInfo.z = small;
					hInfo.front = true;
					hInfo.p = ray.p + small * ray.dir;
					hInfo.N = hInfo.p;
					float u = (1 / (2 * static_cast<float>(M_PI))) * atan2f(hInfo.p.y, hInfo.p.x) + .5f;
					float v = (1 / static_cast<float>(M_PI)) * asinf(hInfo.p.z) + 0.5f;
					hInfo.uvw = Vec3f(u, v, 0.0f);
					return true;
				}
			}
		}

	}
	return false;
}

bool Plane::IntersectRay(Ray const & ray, HitInfo & hInfo, int hitSide) const
{
	float t = -1 * (ray.p.Dot(Vec3f(0, 0, 1))) / ray.dir.Dot(Vec3f(0, 0, 1));
	Vec3f point = ray.p + t * ray.dir;
	if (t < 0)
		return false;
	if (point.x >= -1 && point.x <= 1)
	{
		if (point.y >= -1 && point.y <= 1)
		{
			// this is for shadowprocess
			if (hitSide == 1)
			{
				if ((t > SHADOWBIAS) & (t <= hInfo.z))
				{
					return true;
				}
				else
				{
					return false;
				}
			}

			if (hInfo.z > t ? true : false)
			{
				hInfo.front = true;
				hInfo.p = point;
				hInfo.z = t;
				hInfo.N = Vec3f(0, 0, 1);
				hInfo.uvw = Vec3f(0.5f * point.x + 0.5f, 0.5f * point.y + 0.5f, point.z);

				//float t2 = -1 * (ray.p + hInfo.duvw[0]).Dot(Vec3f(0, 0, 1)) / ray.dir.Dot(Vec3f(0, 0, 1));
				//Vec3f point2 = ray.p + hInfo.duvw[0] + t2 * (ray.dir);
				//hInfo.duvw[0] = 1.0f * (point2 - point);

				//float t3 = -1 * (ray.p + hInfo.duvw[1]).Dot(Vec3f(0, 0, 1)) / ray.dir.Dot(Vec3f(0, 0, 1));
				//Vec3f point3 = ray.p + hInfo.duvw[1] + t3 * (ray.dir);
				//hInfo.duvw[1] = 1.0f * (point3 - point);

				return true;
			}
		}
	}
	return false;
}

bool CheckBoxCollision(const float* vertices, Ray const & ray, float & answer)
{
	// The t value for each case
	float t_xmin, t_ymin, t_zmin;
	float t_xmax, t_ymax, t_zmax;
	float t_tmpmin, t_tmpmax;

	//if (ray.dir.x == 0)
	//{
	//	t_ymin = (vertices[1] - ray.p.y) / ray.dir.y;
	//	t_ymax = (vertices[4] - ray.p.y) / ray.dir.y;

	//	if (t_ymin > t_ymax)
	//	{
	//		SwapFloat(t_ymin, t_ymax);
	//	}

	//	t_zmin = (vertices[2] - ray.p.z) / ray.dir.z;
	//	t_zmax = (vertices[5] - ray.p.z) / ray.dir.z;

	//	if (t_zmin > t_zmax)
	//	{
	//		SwapFloat(t_zmin, t_zmax);
	//	}

	//	if (t_ymin > t_zmax || t_ymin > t_zmax)
	//	{
	//		return false;
	//	}
	//	return true;
	//} 
	//else if (ray.dir.y == 0)
	//{
	//	t_xmin = (vertices[0] - ray.p.x) / ray.dir.x;
	//	t_xmax = (vertices[3] - ray.p.x) / ray.dir.x;

	//	if (t_xmin > t_xmax)
	//	{
	//		SwapFloat(t_xmin, t_xmax);
	//	}

	//	t_zmin = (vertices[2] - ray.p.z) / ray.dir.z;
	//	t_zmax = (vertices[5] - ray.p.z) / ray.dir.z;

	//	if (t_zmin > t_zmax)
	//	{
	//		SwapFloat(t_zmin, t_zmax);
	//	}

	//	if (t_xmin > t_zmax || t_zmin > t_xmax)
	//	{
	//		return false;
	//	}
	//	return true;
	//} 
	//else if (ray.dir.z == 0)
	//{
	//	t_xmin = (vertices[0] - ray.p.x) / ray.dir.x;
	//	t_xmax = (vertices[3] - ray.p.x) / ray.dir.x;

	//	if (t_xmin > t_xmax)
	//	{
	//		SwapFloat(t_xmin, t_xmax);
	//	}

	//	t_ymin = (vertices[1] - ray.p.y) / ray.dir.y;
	//	t_ymax = (vertices[4] - ray.p.y) / ray.dir.y;

	//	if (t_ymin > t_ymax)
	//	{
	//		SwapFloat(t_ymin, t_ymax);
	//	}

	//	if (t_xmin > t_ymax || t_ymin > t_xmax)
	//	{
	//		return false;
	//	}
	//	return true;
	//}
	//else
	{
		t_xmin = (vertices[0] - ray.p.x) / ray.dir.x;
		t_xmax = (vertices[3] - ray.p.x) / ray.dir.x;

		if (t_xmin > t_xmax)
		{
			SwapFloat(t_xmin, t_xmax);
		}

		t_ymin = (vertices[1] - ray.p.y) / ray.dir.y;
		t_ymax = (vertices[4] - ray.p.y) / ray.dir.y;

		if (t_ymin > t_ymax)
		{
			SwapFloat(t_ymin, t_ymax);
		}

		if (t_xmin > t_ymax || t_ymin > t_xmax)
		{
			return false;
		}

		t_tmpmin = (t_xmin > t_ymin) ? t_ymin : t_xmin;
		t_tmpmax = (t_xmax > t_ymax) ? t_xmax : t_ymax;

		t_zmin = (vertices[2] - ray.p.z) / ray.dir.z;
		t_zmax = (vertices[5] - ray.p.z) / ray.dir.z;

		if (t_zmin > t_zmax)
		{
			SwapFloat(t_zmin, t_zmax);
		}

		if (t_tmpmin > t_zmax || t_zmin > t_tmpmax)
		{
			return false;
		}
		return true;
	}
}

bool TriObj::IntersectRay(Ray const & ray, HitInfo & hInfo, int hitSide) const
{
	const float * vertices = bvh.GetNodeBounds(bvh.GetRootNodeID());
	float answer = 0.0f;
	if (CheckBoxCollision(vertices, ray, answer))
	{
		return TraceBVHNode(ray, hInfo, hitSide, bvh.GetRootNodeID());
	}
	else
	{
		return false;
	}
}

bool TriObj::IntersectTriangle(Ray const & ray, HitInfo & hInfo, int hitSide, unsigned int faceID) const
{
	int i0 = F(faceID).v[0];
	int i1 = F(faceID).v[1];
	int i2 = F(faceID).v[2];

	Vec3f n = (v[i2] - v[i0]).Cross(v[i1] - v[i0]);
	float h = -1 * n.Dot(v[i0]);
	float t = -1 * (ray.p.Dot(n) + h) / ray.dir.Dot(n);

	// This is for shadow process
	if (hitSide == 1)
	{
		if (t - SHADOWBIAS < 0 || t > hInfo.z)
			return false;
	}

	if (t < 0)
		return false;

	Vec3f point = ray.p + t * ray.dir;

	Vec2f v0, v1, v2, x;
	if (std::abs(n.x) >= std::abs(n.y) && std::abs(n.x) >= std::abs(n.z))
	{
		v0 = Vec2f(v[i0].y, v[i0].z); v1 = Vec2f(v[i1].y, v[i1].z); v2 = Vec2f(v[i2].y, v[i2].z);
		x = Vec2f(point.y, point.z);
	}
	else if (std::abs(n.y) >= std::abs(n.x) && std::abs(n.y) >= std::abs(n.z))
	{
		v0 = Vec2f(v[i0].x, v[i0].z); v1 = Vec2f(v[i1].x, v[i1].z); v2 = Vec2f(v[i2].x, v[i2].z);
		x = Vec2f(point.x, point.z);
	}
	else if (std::abs(n.z) >= std::abs(n.x) && std::abs(n.z) >= abs(n.y))
	{
		v0 = Vec2f(v[i0].x, v[i0].y); v1 = Vec2f(v[i1].x, v[i1].y); v2 = Vec2f(v[i2].x, v[i2].y);
		x = Vec2f(point.x, point.y);
	}

	float a0, a1, a2;
	a0 = (v1 - x).Cross(v2 - x);
	a1 = (v2 - x).Cross(v0 - x);
	a2 = (v0 - x).Cross(v1 - x);

	if ((a0 >= 0 && a1 >= 0 && a2 >= 0) || (a0 < 0 && a1 < 0 && a2 < 0))
	{
		// This is for shadow process
		if (hitSide == 1)
		{
			return true;
		}

		if (hInfo.z > t ? true : false)
		{
			float a = (v1 - v0).Cross(v2 - v0);
			float beta0, beta1, beta2;
			beta0 = std::abs(a0 / a);
			beta1 = std::abs(a1 / a);
			beta2 = std::abs(a2 / a);

			Vec3f bc = Vec3f(beta0, beta1, beta2);
			Vec3f normal = GetNormal(faceID, bc);
			
			if (a0 >= 0 && a1 >= 0 && a2 >= 0)
			{
				hInfo.front = true;
			}
			else if (a0 < 0 && a1 < 0 && a2 < 0)
			{
				hInfo.front = false;
			}
			hInfo.N = normal;
			hInfo.p = point;
			hInfo.z = t;
			// This is for texturing
			{
				hInfo.uvw = GetTexCoord(faceID, bc);
			}
			return true;
		}
	}
	return false;
}

bool TriObj::TraceBVHNode(Ray const & ray, HitInfo & hInfo, int hitSide, unsigned int nodeID) const
{
	if (bvh.IsLeafNode(nodeID))
	{
		unsigned int nodecount = bvh.GetNodeElementCount(nodeID);
		const unsigned int* nodeelementslist = bvh.GetNodeElements(nodeID);

		bool hit = false;
		for (int i = 0; i < (signed)nodecount; i++)
		{
			if (IntersectTriangle(ray, hInfo, hitSide, nodeelementslist[i]))
			{
				hit = true;
			}
		}

		return hit;
	}
	else
	{
		unsigned int child1id = bvh.GetFirstChildNode(nodeID);
		unsigned int child2id = bvh.GetSecondChildNode(nodeID);
		const float * vertices1 = bvh.GetNodeBounds(child1id);
		const float * vertices2 = bvh.GetNodeBounds(child2id);

		float t_child1, t_child2;

		if (!CheckBoxCollision(vertices1, ray, t_child1) & !(CheckBoxCollision(vertices2, ray, t_child2)))
		{
			return false;
		}

		if (t_child1 < t_child2)
		{
			if (t_child1 <= hInfo.z)
			{
				bool hitcheck = false;
				if (TraceBVHNode(ray, hInfo, hitSide, child1id))
				{
					hitcheck = true;

					if (t_child2 <= hInfo.z)
					{
						if (TraceBVHNode(ray, hInfo, hitSide, child2id))
						{
							hitcheck = true;
						}
					}
				}
				else if (t_child2 <= hInfo.z)
				{
					if (TraceBVHNode(ray, hInfo, hitSide, child2id))
					{
						hitcheck = true;
					}
				}
				return hitcheck;
			}
		}
		else if (t_child1 > t_child2)
		{
			if (t_child2 <= hInfo.z)
			{
				bool hitcheck = false;
				if (TraceBVHNode(ray, hInfo, hitSide, child2id))
				{
					hitcheck = true;

					if (t_child1 <= hInfo.z)
					{
						if (TraceBVHNode(ray, hInfo, hitSide, child1id))
						{
							hitcheck = true;
						}
					}
				}
				else if (t_child1 <= hInfo.z)
				{
					if (TraceBVHNode(ray, hInfo, hitSide, child1id))
					{
						hitcheck = true;
					}
				}
				return hitcheck;
			}
		}
		else if (t_child1 == t_child2)
		{
			bool hitcheck = false;
			if (TraceBVHNode(ray, hInfo, hitSide, child1id))
			{
				hitcheck = true;
			}
			if (TraceBVHNode(ray, hInfo, hitSide, child2id))
			{
				hitcheck = true;
			}
			return hitcheck;
		}
	}
	return false;
}