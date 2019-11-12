#pragma once
#define _USE_MATH_DEFINES

#include <math.h>
#include <cmath>
#include <cyCodeBase/cyVector.h>

namespace Utility
{
	using cy::Vec3f;

	class utility
	{
	public:
		static void ConeUniformSampling(Vec3f & dir, const float radius);
	};

	inline void utility::ConeUniformSampling(Vec3f & dir, const float radius)
	{
		Vec3f u = Vec3f(-1 * dir.y, dir.x, 0);
		Vec3f v = dir.Cross(u);
		u.Normalize();
		v.Normalize();

		float randomrand = 2 * static_cast<float>(M_PI) * (static_cast<float>(rand()) / (RAND_MAX));
		Vec3f randomradiusdir = sinf(randomrand) * u + cosf(randomrand) * v;

		float randomradiuslength = (static_cast<float>(rand()) / (RAND_MAX));
		randomradiuslength = radius * sqrt(randomradiuslength);

		dir += randomradiuslength * randomradiusdir;
		dir.Normalize();
	}
}

