#pragma once
#include <constant.h>
#include <scene.h>
#include <utility.h>

extern Node rootNode;
extern MaterialList materials;
extern TexturedColor environment;
extern LightList lights;

void RayTraversing(Node * node, Ray ray, HitInfo & hit);

namespace ReflectionAndRefraction 
{
	using namespace Utility;

	Color ReflectionRefractionTraverse(Ray ray, int bounce)
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
		if (bounce > REFLECTIONBOUNCE)
		{
			return Color(0, 0, 0);
		}
		else
		{
			bounce++;
			Vec3f N = hInfo.N;

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

			Color returnColor = ReflectionRefractionTraverse(S, bounce);
			return returnColor;
		}
	}

	Color Refraction(Ray const & ray, const HitInfo & hit, int bounce, float refractionIndex, Color refraction, const float glossiness)
	{
		if (bounce > REFRACTIONBOUNCE)
		{
			return Color(0, 0, 0);
		}
		else
		{
			bounce++;
			Vec3f V = -1 * ray.dir;

			Vec3f N = hit.N;
			ConeUniformSampling(N, glossiness);

			Vec3f P;

			float cos1, cos2, sin1, sin2;

			if (V.Dot(N) >= 0)
			{
				P = -1 * SHADOWBIAS * N;  P += hit.p;
				cos1 = V.Dot(N);
			}
			else
			{
				P = SHADOWBIAS * N;  P += hit.p;
				cos1 = V.Dot(-N);
			}

			sin1 = sqrt(1 - (cos1 * cos1));

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

				// Calculate Fresnel reflection
				float R = FresnelReflections(hit, refractionIndex, cos1);

				// Refraction part
				Color returnColor = (1 - R) * refraction * ReflectionRefractionTraverse(S, bounce);
				// Reflection part
				returnColor += R * refraction * Reflection(ray, hit, bounce, glossiness);
				return returnColor;
			}
		}
	}
}