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

		return R0 + (1 - R0) * powf(1 - cos1, 5);
	}

	Color Reflection(Ray const & ray, const HitInfo & hInfo, int bounce, const float alpha)
	{
		Vec3f N = hInfo.N;

		float theta = GetUniformRamdomFloat();
		float phy = GetUniformRamdomFloat();
		Vec3f D_dash = SpecularWeightedHemisphereSampling(N, alpha, theta, phy);

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

	Color Refraction(const Ray & ray, const HitInfo & hit, int bounce, Color refractioncolor, const float & refractionindex, const float & alpha)
	{
		Vec3f V = -1 * ray.dir;

		Vec3f N = hit.N;

		float theta = GetUniformRamdomFloat();
		float phy = GetUniformRamdomFloat();
		SpecularWeightedHemisphereSampling(N, alpha, theta, phy);

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
			sin2 = (1 / refractionindex) * sin1;
		}
		else
		{
			sin2 = refractionindex * sin1;
		}

		// S is a starting point from the surface point
		Ray S;

		// Total internal reflection
		if (sin2 > 1)
		{
			return refractioncolor * Reflection(ray, hit, bounce, alpha);
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
			float R = FresnelReflections(hit, refractionindex, cos1);

			float rand = GetUniformRamdomFloat();
			Color returncolor = Color(0, 0, 0);

			if (rand <= R)
			{
				// Reflection part
				returncolor = refractioncolor * Reflection(ray, hit, bounce, alpha);

			}
			else
			{
				// Refraction part
				returncolor = refractioncolor * ReflectionRefractionTraverse(S, bounce);
			}
			return returncolor;
		}
	}
} //<- End of namespace