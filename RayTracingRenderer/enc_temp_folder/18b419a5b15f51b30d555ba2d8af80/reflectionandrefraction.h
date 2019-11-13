#pragma once
#include <constant.h>
#include <scene.h>
#include <utility.h>

extern Node rootNode;
extern MaterialList materials;
extern TexturedColor environment;
extern LightList lights;

void RayTraversing(Node * node, Ray ray, HitInfo & hit, int bounce);

namespace ReflectionAndRefraction 
{
	using namespace Utility;

	Color ReflectionRefractionTraverse(Ray ray, int bounce)
	{
		Node * node = &rootNode;
		HitInfo hit = HitInfo();

		RayTraversing(node, ray, hit, bounce);

		if (hit.node != nullptr)
		{
			return materials.Find(hit.node->GetMaterial()->GetName())->Shade(ray, hit, lights, bounce);
		}
		else
		{
			return environment.SampleEnvironment(ray.dir);
		}
	}

	Color CalculateAbsorption(Color in, Color absorption, float distance)
	{
		float rout = exp(-1 * absorption.r * distance) * in.r;
		float gout = exp(-1 * absorption.g * distance) * in.g;
		float bout = exp(-1 * absorption.b * distance) * in.b;
		return Color(rout, gout, bout);
	}

	float FresnelReflections(const HitInfo & hInfo, float refractionIndex, float cos1)
	{
		float R0;
		if (hInfo.front)
		{
			R0 = ((1 - refractionIndex) * (1 - refractionIndex))
				/ //----------------------------------------------
				((1 + refractionIndex) * (1 + refractionIndex));
		}
		else
		{
			R0 = ((refractionIndex - 1) * (refractionIndex - 1))
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

			Color returnColor = ReflectionRefractionTraverse(S, bounce - 1);
			return returnColor;
		}
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



			// Calculate Fresnel reflection
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
				return ReflectionRefractionTraverse(S, bounce);
			}
			else // Normal procedure
			{
				cos2 = sqrt(1 - (sin2 * sin2));

				if (V.Dot(N) >= 0)
				{
					// Horizontal direction Vector
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

				// Refraction part
				Color returnColor = (1 - R) * refraction * ReflectionRefractionTraverse(S, bounce - 1);
				// Reflection part
				returnColor += R * refraction * Reflection(ray, hInfo, bounce, glossiness);
				return returnColor;
			}
		}
	}

}