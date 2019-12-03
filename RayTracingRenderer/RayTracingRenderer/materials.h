
//-------------------------------------------------------------------------------
///
/// \file       materials.h 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    13.0
/// \date       August 21, 2019
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifndef _MATERIALS_H_INCLUDED_
#define _MATERIALS_H_INCLUDED_

#include "scene.h"
#include "constant.h"

namespace Utility 
{
	float GetUniformRamdomFloat();
	Vec3f CosineWeightedHemisphereUniformSampling(const Vec3f &, float, float);
	Vec3f SpecularWeightedHemisphereSampling(const Vec3f &, const float &, float, float);
}

//-------------------------------------------------------------------------------

class MtlBlinn : public Material
{
public:
	MtlBlinn() : diffuse(0.5f, 0.5f, 0.5f), specular(0.7f, 0.7f, 0.7f), glossiness(20.0f), emission(0, 0, 0),
		reflection(0, 0, 0), refraction(0, 0, 0), absorption(0, 0, 0), ior(1),
		reflectionGlossiness(0), refractionGlossiness(0)
	{
	}
	virtual Color Shade(Ray const &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const;

	void SetDiffuse(Color dif)
	{
		diffuse.SetColor(dif);
	}
	void SetSpecular(Color spec)
	{
		specular.SetColor(spec);
	}
	void SetGlossiness(float gloss)
	{
		glossiness = gloss;
	}
	void SetEmission(Color e)
	{
		emission.SetColor(e);
	}

	void SetReflection(Color reflect)
	{
		reflection.SetColor(reflect);
	}
	void SetRefraction(Color refract)
	{
		refraction.SetColor(refract);
	}
	void SetAbsorption(Color absorp)
	{
		absorption = absorp;
	}
	void SetRefractionIndex(float _ior)
	{
		ior = _ior;
	}

	void SetDiffuseTexture(TextureMap *map)
	{
		diffuse.SetTexture(map);
	}
	void SetSpecularTexture(TextureMap *map)
	{
		specular.SetTexture(map);
	}
	void SetEmissionTexture(TextureMap *map)
	{
		emission.SetTexture(map);
	}
	void SetReflectionTexture(TextureMap *map)
	{
		reflection.SetTexture(map);
	}
	void SetRefractionTexture(TextureMap *map)
	{
		refraction.SetTexture(map);
	}
	void SetReflectionGlossiness(float gloss)
	{
		reflectionGlossiness = gloss;
	}
	void SetRefractionGlossiness(float gloss)
	{
		refractionGlossiness = gloss;
	}

	virtual void SetViewportMaterial(int subMtlID = 0) const; // used for OpenGL display

	// Photon Extensions
	virtual bool IsPhotonSurface(int subMtlID = 0) const
	{
		return diffuse.GetColor().Gray() > 0;
	} // if this method returns true, the photon will be stored
	virtual bool RandomPhotonBounce(Ray &r, Color &c, const HitInfo &hInfo) const  // if this method returns true, a new photon with the given direction and color will be traced
	{
		float rand = Utility::GetUniformRamdomFloat();

		float xi1 = Utility::GetUniformRamdomFloat();
		float xi2 = Utility::GetUniformRamdomFloat();

		if (0 < rand && rand <= diffuse.GetColor().Max())
		{
			c *= 1 / diffuse.GetColor().Max() * diffuse.GetColor();

			Vec3f N_dash = Utility::CosineWeightedHemisphereUniformSampling(hInfo.N, xi1, xi2);

			Ray ray_gi;
			ray_gi.p = hInfo.p;
			ray_gi.p += SHADOWBIAS * hInfo.N;
			ray_gi.dir = N_dash;
			r = ray_gi;
			return true;
		}
		else if (diffuse.GetColor().Max() < rand && rand <= specular.GetColor().Max() + diffuse.GetColor().Max())
		{
			float rand2 = Utility::GetUniformRamdomFloat();
			if (0 < rand2 && rand2 <= reflection.GetColor().Max())
			{
				c *= 1 / reflection.GetColor().Max() * reflection.GetColor();

				Vec3f V = -1 * r.dir;
				Vec3f R = 2 * (hInfo.N.Dot(V)) * hInfo.N - V;
				R.Normalize();

				//Vec3f D_dash = SpecularWeightedHemisphereSampling(R, glossiness, xi1, xi2);

				Ray ray_gi;
				ray_gi.p = hInfo.p;
				ray_gi.p += SHADOWBIAS * hInfo.N;
				ray_gi.dir = R;
				r = ray_gi;
				return true;
			}
			else if (reflection.GetColor().Max() < rand2 && rand2 <= reflection.GetColor().Max() + refraction.GetColor().Max())
			{
				c *= 1 / refraction.GetColor().Max() * refraction.GetColor();

				Vec3f V = -1 * r.dir;

				Vec3f N = hInfo.N;

				//float xi1 = GetUniformRamdomFloat();
				//float xi2 = GetUniformRamdomFloat();
				//SpecularWeightedHemisphereSampling(N, alpha, xi1, xi2);

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

				Vec3f T_h, T_v, T;

				if (V.Dot(N) >= 0)
				{
					sin2 = (1 / ior) * sin1;
				}
				else
				{
					sin2 = ior * sin1;
				}

				// S is a starting point from the surface point
				Ray S;

				// Total internal reflection
				if (sin2 > 1)
				{
					return false;
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
					//float R = FresnelReflections(hit, refractionindex, cos1);

					//float rand = GetUniformRamdomFloat();
					//Color returncolor = Color(0, 0, 0);
				}

				Ray ray_gi;
				ray_gi = S;
				return true;
			}
			else
			{
				c *= 1 / (1 - reflection.GetColor().Max() - refraction.GetColor().Max()) * specular.GetColor();

				Vec3f V = -1 * r.dir;
				Vec3f R = 2 * (hInfo.N.Dot(V)) * hInfo.N - V;
				R.Normalize();

				Vec3f D_dash = Utility::SpecularWeightedHemisphereSampling(R, glossiness, xi1, xi2);

				Ray ray_gi;
				ray_gi.p = hInfo.p;
				ray_gi.p += SHADOWBIAS * hInfo.N;
				ray_gi.dir = D_dash;
				r = ray_gi;
			}
			return true;
		}
		else
		{
			return false;
		}
	}

private:
	TexturedColor diffuse, specular, reflection, refraction, emission;
	float glossiness;
	Color absorption;
	float ior;  // index of refraction
	float reflectionGlossiness, refractionGlossiness;
};

//-------------------------------------------------------------------------------

class MultiMtl : public Material
{
public:
	virtual ~MultiMtl()
	{
		for (unsigned int i = 0; i < mtls.size(); i++) delete mtls[i];
	}

	virtual Color Shade(Ray const &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const
	{
		return hInfo.mtlID < (int)mtls.size() ? mtls[hInfo.mtlID]->Shade(ray, hInfo, lights, bounceCount) : Color(1, 1, 1);
	}

	virtual void SetViewportMaterial(int subMtlID = 0) const
	{
		if (subMtlID < (int)mtls.size()) mtls[subMtlID]->SetViewportMaterial();
	}

	void AppendMaterial(Material *m)
	{
		mtls.push_back(m);
	}

	// Photon Extensions
	virtual bool IsPhotonSurface(int subMtlID = 0) const
	{
		return mtls[subMtlID]->IsPhotonSurface();
	}
	virtual bool RandomPhotonBounce(Ray &r, Color &c, const HitInfo &hInfo) const
	{
		return hInfo.mtlID < (int)mtls.size() ? mtls[hInfo.mtlID]->RandomPhotonBounce(r, c, hInfo) : false;
	}

private:
	std::vector<Material*> mtls;
};

//-------------------------------------------------------------------------------

#endif