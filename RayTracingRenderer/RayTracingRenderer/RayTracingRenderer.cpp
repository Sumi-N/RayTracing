// RayTracingRenderer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "constant.h"

#include <iostream>
#include <scene.h>
#include <xmlload.h>
#include <viewport.h>
#include <objects.h>
#include <omp.h>
#include <cyBVH.h>
#include <time.h>
#include <math.h>

#include "utility.h"
#include "reflectionandrefraction.h"

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
using namespace ReflectionAndRefraction;

int main()
{
	//LoadScene(".\\xmlfiles\\playground5.xml");
	//LoadScene(".\\xmlfiles\\catscene.xml");
	//LoadScene(".\\xmlfiles\\potscene.xml");
	//LoadScene(".\\xmlfiles\\playground6.xml");
	//LoadScene(".\\xmlfiles\\assignment11_2.xml");
	//LoadScene(".\\xmlfiles\\bosonscene.xml");
	LoadScene(".\\xmlfiles\\assignment11.xml");
	//LoadScene(".\\xmlfiles\\assignment4.xml");
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

void RayTraversing(Node * node, Ray ray, HitInfo & hit)
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
			RayTraversing(currentnode, currentray, hit);
			if (hit.node != nullptr && hit.node != tmphit.node)
			{
				currentnode->FromNodeCoords(hit);
				tmphit = hit;
			}
		}
	}

	if (node == &rootNode)
	{
		hit = tmphit;
	}
}

Color InitialRayTraverse(Ray ray)
{
	Node * node = &rootNode;
	HitInfo hit = HitInfo();

	RayTraversing(node, ray, hit);

	if (hit.node != nullptr)
	{
		Color returnColor = materials.Find(hit.node->GetMaterial()->GetName())->Shade(ray, hit, lights, 0);
		return returnColor;
	}
	else
	{
		return Color(0, 0, 0);
	}
}

Color PathTracingTraverse(Ray ray)
{
	Vec3f screenpoint = ray.p + camera.focaldist * ray.dir;
	screenpoint -= HALF * (pixelx + pixely);
	Color returnColor = Color(0, 0, 0);

	for (int k = 0; k < RAYPERPT; k++)
	{
		Vec3f offset = Vec3f(0, 0, 0);
		SquareUniformSampling(offset, pixelx, pixely, 2 * HALF);
		Vec3f fixedpoint = screenpoint + offset;

		ray.dir = fixedpoint - ray.p;
		ray.Normalize();

		returnColor += InitialRayTraverse(ray);
	}
	returnColor /= RAYPERPT;
	return returnColor;
}

Color GlobalIlluminationTraverse(Ray ray, int bounce)
{
	Node * node = &rootNode;
	HitInfo hit = HitInfo();

	RayTraversing(node, ray, hit);
	
	if (hit.node != nullptr)
	{
		return materials.Find(hit.node->GetMaterial()->GetName())->Shade(ray, hit, lights, bounce);
	}
	else
	{
		return environment.SampleEnvironment(ray.dir);
	}
}

Color BlurEffect(Ray ray)
{
	Ray cameraraies[RAYPERPIXELFORBLUREFFECT];
	HitInfo hits[RAYPERPIXELFORBLUREFFECT];

	Color averagepixelcolor = Color(0, 0, 0);
	Color pixelcolors[RAYPERPIXELFORBLUREFFECT];
	Vec3f screenpoints[RAYPERPIXELFORBLUREFFECT];

	Color returnColor = Color(0, 0, 0);

	for (int i = 0; i < RAYPERPIXELFORBLUREFFECT; i++)
	{
		screenpoints[i] = ray.p + camera.focaldist * ray.dir;

		Vec3f offset = Vec3f(0, 0, 0);
		CircleUniformSampling(offset, camera.dir, camera.up, camera.dof);
		cameraraies[i].p = ray.p + offset;

		cameraraies[i].dir = screenpoints[i] - cameraraies[i].p;
		cameraraies[i].Normalize();

		pixelcolors[i] = InitialRayTraverse(cameraraies[i]);
		returnColor += pixelcolors[i];
	}

	return returnColor / RAYPERPIXELFORBLUREFFECT;
}
void BeginRender() {

	time_t time0;   // create timers.
	time_t time1;

	time(&time0);   // get current time.

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
	for (int j = 0; j < renderImage.GetWidth(); j++)
	{
		for (int i = 0; i < renderImage.GetHeight(); i++) 
		{

			cameraray[i * renderImage.GetWidth() + j].dir = f + (j + HALF) * (w / W)*x - (i + HALF) * (h / H)*y - camera.pos;
			cameraray[i * renderImage.GetWidth() + j].p = camera.pos;
			zbuffers[i * renderImage.GetWidth() + j] = BIGFLOAT;
			samplecount[i * renderImage.GetWidth() + j] = 0;
			cameraray[i * renderImage.GetWidth() + j].Normalize();

			Color resultColor = Color(0, 0, 0);

		#ifdef ENABLEPT
			resultColor = PathTracingTraverse(cameraray[i * renderImage.GetWidth() + j]);
			//resultColor = PathTracingTraverse(cameraray[i * renderImage.GetWidth() + j]);
		#elif defined BLUREFFECT
			resultColor = BlurEffect(cameraray[i * renderImage.GetWidth() + j]);		
		#elif
			resultColor = InitialRayTraverse(cameraray[i * renderImage.GetWidth() + j]);
			//zbuffers[i * renderImage.GetWidth() + j] = hit.z;
		#endif 

			if (isnan(resultColor.r) || isnan(resultColor.g) || isnan(resultColor.b))
			{	
				printf("There is some NaN in the calculation process\n");
			}

			// gamma correction
		#ifdef ENABLEGAMMA
				resultColor.r = powf(resultColor.r, 1 / 2.2f);
				resultColor.g = powf(resultColor.g, 1 / 2.2f);
				resultColor.b = powf(resultColor.b, 1 / 2.2f);
		#endif

			pixels[i * renderImage.GetWidth() + j] = (Color24)resultColor;

			if (pixels[i * renderImage.GetWidth() + j] == (Color24)Color(0, 0, 0))
			{
				Vec3f v((float)j / renderImage.GetWidth(), (float)i / renderImage.GetHeight(), 0.0f);
				pixels[i * renderImage.GetWidth() + j] = (Color24)background.Sample(v);
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

Color MtlBlinn::Shade(Ray const & ray, const HitInfo & hInfo, const LightList & lights, int bounce) const
{
	Vec3f N = hInfo.N;
	Color color = Color(0, 0, 0);
	Color specularpart = Color(0, 0, 0);
	Color diffusepart = Color(0, 0, 0);

	if (ray.p == hInfo.p)
	{
		return color;
	}

	for (auto light = lights.begin(); light != lights.end(); ++light)
	{
		if (strcmp((*light)->GetName(), "directLight") == 0)
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
			Color diffusepart = this->diffuse.Sample(hInfo.uvw, hInfo.duvw);
			if (hInfo.N.Dot(L) != 0)
			{
				specularpart = ((this->glossiness + 2) / 2) * oneofcos * powf(H.Dot(hInfo.N), this->glossiness) * this->specular.Sample(hInfo.uvw, hInfo.duvw);
			}
			else
			{
				specularpart = Color(0, 0, 0);
			}
			ClampColorValue(specularpart);
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
			diffusepart = (1 / static_cast<float>(M_PI)) * this->diffuse.Sample(hInfo.uvw, hInfo.duvw);
			if (hInfo.N.Dot(L) != 0)
			{
				specularpart = ((this->glossiness + 2) / 2) * oneofcos * powf(H.Dot(hInfo.N), this->glossiness) * this->specular.Sample(hInfo.uvw, hInfo.duvw);
			}
			else
			{
				specularpart = Color(0, 0, 0);
			}
			ClampColorValue(specularpart);
			color += (diffusepart + specularpart) * IR;
		}
	}
#ifdef ENABLEPT

	float theta = GetUniformRamdomFloat();
	float phy = GetUniformRamdomFloat();

	float rand = GetUniformRamdomFloat();
	float rand2 = GetUniformRamdomFloat();
	float rand3 = GetUniformRamdomFloat();

#ifdef ENABLEGIMIS

	float pdf_diffuse = sinf(static_cast<float>(M_PI) / 2 * rand);
	float pdf_specular = ((this->glossiness + 2) / 2) * cosf(pow(static_cast<float>(M_PI) / 2 * rand2, 1 / (this->glossiness + 1)));
	float borderline = (pdf_diffuse) / (pdf_diffuse + pdf_specular);

	if (bounce < GIBOUNCE)
	{
		int bouncetime = bounce + 1;
		Color returnColor = Color(0, 0, 0);

		// Diffuse part for GI
		if (rand3 < borderline)
		{
			if (this->diffuse.Sample(hInfo.uvw) != Color(0, 0, 0))
			{
				Vec3f N_dash = CosineWeightedHemisphereUniformSampling(N, theta, phy);

				Ray ray_gi;
				ray_gi.p = hInfo.p;
				ray_gi.p += SHADOWBIAS * N;
				ray_gi.dir = N_dash;

				returnColor = this->diffuse.Sample(hInfo.uvw, hInfo.duvw) * GlobalIlluminationTraverse(ray_gi, bouncetime);		
			}
		}
		else
		{
			// Specular part for GI
			if (this->diffuse.Sample(hInfo.uvw) != Color(0, 0, 0))
			{
				Vec3f V = -1 * ray.dir;
				Vec3f R = 2 * (N.Dot(V)) * N - V;
				R.Normalize();

				Vec3f D_dash = SpecularWeightedHemisphereSampling(R, this->glossiness, theta, phy);

				Ray ray_gi;
				ray_gi.p = hInfo.p;
				ray_gi.p += SHADOWBIAS * N;
				ray_gi.dir = D_dash;

				returnColor = this->specular.Sample(hInfo.uvw, hInfo.duvw) * GlobalIlluminationTraverse(ray_gi, bouncetime);
			}
		}
		color += returnColor;
	}

#else
	// Diffuse part for GI
	if (bounce < GIBOUNCE)
	{
		int bouncetime = bounce + 1;
		Color returnColor = Color(0, 0, 0);

		Vec3f N_dash = CosineWeightedHemisphereUniformSampling(N, theta, phy);

		Ray ray_gi;
		ray_gi.p = hInfo.p;
		ray_gi.p += SHADOWBIAS * N;
		ray_gi.dir = N_dash;

		returnColor += this->diffuse.Sample(hInfo.uvw, hInfo.duvw) * GlobalIlluminationTraverse(ray_gi, bouncetime);


		Vec3f V = -1 * ray.dir;
		Vec3f R = 2 * (N.Dot(V)) * N - V;
		R.Normalize();

		Vec3f D_dash = SpecularWeightedHemisphereSampling(R, this->glossiness, theta, phy);

		Ray ray_gi2;
		ray_gi2.p = hInfo.p;
		ray_gi2.p += SHADOWBIAS * N;
		ray_gi2.dir = D_dash;

		returnColor += this->specular.Sample(hInfo.uvw, hInfo.duvw) * GlobalIlluminationTraverse(ray_gi2, bouncetime);
		color += returnColor;
	}
#endif // End Enable GIMIS
#endif // End Enable PT

	// Calculate refraction part
	if (this->refraction.Sample(hInfo.uvw) != Color(0, 0, 0))
	{
		if (bounce < GIBOUNCE)
		{
			int bouncetime = bounce + 1;

			//// When it is a back side hit, it means that absorption gonna happen during inside the material the light go through
			//if (!hInfo.front)
			//{
			//	color += Color(exp(-1 * absorption.r * hInfo.z) * color.r, exp(-1 * absorption.g * hInfo.z) * color.g, exp(-1 * absorption.b * hInfo.z) * color.b);
			//}

			color += Refraction(ray, hInfo, bouncetime, refraction.Sample(hInfo.uvw), ior, refractionGlossiness);
		}
	}
	if (isnan(color.r) || isnan(color.g) || isnan(color.b))
	{		
		return Color(0, 0, 0);
	}
	else
	{
		return color;
	}
	
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
				// this is for shadow process
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
			// this is for shadow process
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
			// this is for shadow process
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

	if (ray.dir.x == 0)
	{
		t_ymin = (vertices[1] - ray.p.y) / ray.dir.y;
		t_ymax = (vertices[4] - ray.p.y) / ray.dir.y;

		if (t_ymin > t_ymax)
		{
			SwapFloat(t_ymin, t_ymax);
		}

		t_zmin = (vertices[2] - ray.p.z) / ray.dir.z;
		t_zmax = (vertices[5] - ray.p.z) / ray.dir.z;

		if (t_zmin > t_zmax)
		{
			SwapFloat(t_zmin, t_zmax);
		}

		if (t_ymin > t_zmax || t_ymin > t_zmax)
		{
			return false;
		}
		return true;
	}
	else if (ray.dir.y == 0)
	{
		t_xmin = (vertices[0] - ray.p.x) / ray.dir.x;
		t_xmax = (vertices[3] - ray.p.x) / ray.dir.x;

		if (t_xmin > t_xmax)
		{
			SwapFloat(t_xmin, t_xmax);
		}

		t_zmin = (vertices[2] - ray.p.z) / ray.dir.z;
		t_zmax = (vertices[5] - ray.p.z) / ray.dir.z;

		if (t_zmin > t_zmax)
		{
			SwapFloat(t_zmin, t_zmax);
		}

		if (t_xmin > t_zmax || t_zmin > t_xmax)
		{
			return false;
		}
		return true;
	}
	else if (ray.dir.z == 0)
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
		return true;
	}
	else
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