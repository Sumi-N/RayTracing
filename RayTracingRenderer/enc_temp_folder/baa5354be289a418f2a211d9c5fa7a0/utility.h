#pragma once
#define _USE_MATH_DEFINES

#include <math.h>
#include <cmath>
#include <cyCodeBase/cyVector.h>

namespace Utility
{
	using namespace cy;

	inline void SwapFloat(float & a, float & b)
	{
		float tmp = a;
		a = b;
		b = tmp;
	}

	inline void SquareUniformSampling(Vec3f & point, Vec3f x_axis, Vec3f y_axis, const float & length)
	{
		float x_length = length * ((float)rand() / (RAND_MAX));
		float y_length = length * ((float)rand() / (RAND_MAX));
		point = x_length * x_axis + y_length * y_axis;
	}

	inline void CircleUniformSampling(Vec3f & point, const Vec3f forward_dir, const Vec3f upward_dir, const float & radius)
	{
		float randomradius = (static_cast<float>(rand()) / (RAND_MAX));
		float randomrand = 2 * static_cast<float>(M_PI) * (static_cast<float>(rand()) / (RAND_MAX));
		randomradius = radius * sqrt(randomradius);

		Vec3f right_dir = forward_dir.Cross(upward_dir);
		point = randomradius * (cos(randomrand) * upward_dir + sin(randomrand) * right_dir);
	}

	inline void ConeUniformSampling(Vec3f & dir, const float & radius)
	{
		Vec3f randomvec = Vec3f((static_cast<float>(rand()) / (RAND_MAX)), (static_cast<float>(rand()) / (RAND_MAX)), (static_cast<float>(rand()) / (RAND_MAX)));
		Vec3f u = dir.Cross(randomvec);
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

	inline void HemisphereUniformSampling(Vec3f & dir)
	{
		Vec3f randomvec = Vec3f((static_cast<float>(rand()) / (RAND_MAX)), (static_cast<float>(rand()) / (RAND_MAX)), (static_cast<float>(rand()) / (RAND_MAX)));
		Vec3f x_dir = dir.Cross(randomvec);
		Vec3f y_dir = dir.Cross(x_dir);
		x_dir.Normalize();
		y_dir.Normalize();

		float theta = (static_cast<float>(rand()) / (RAND_MAX));
		float phy = (static_cast<float>(rand()) / (RAND_MAX));

		float x_length = cosf(2 * static_cast<float>(M_PI) * phy) * sqrtf(1 - (1- theta) * (1 - theta));
		float y_length = sinf(2 * static_cast<float>(M_PI) * phy) * sqrtf(1 - (1 - theta) * (1 - theta));
		float z_length = 1 - theta;

		dir = x_length * x_dir + y_length * y_dir + z_length * dir;
	}

	inline void CosineWeightedHemisphereUniformSampling(Vec3f & dir)
	{
		Vec3f randomvec = Vec3f((static_cast<float>(rand()) / (RAND_MAX)), (static_cast<float>(rand()) / (RAND_MAX)), (static_cast<float>(rand()) / (RAND_MAX)));
		Vec3f x_dir = dir.Cross(randomvec);
		Vec3f y_dir = dir.Cross(x_dir);
		x_dir.Normalize();
		y_dir.Normalize();

		float theta = (static_cast<float>(rand()) / (RAND_MAX));
		float phy = (static_cast<float>(rand()) / (RAND_MAX));

		float x_length = cosf(2 * static_cast<float>(M_PI) * phy) * sqrtf(theta);
		float y_length = sinf(2 * static_cast<float>(M_PI) * phy) * sqrtf(theta);
		float z_length = sqrtf(1 - theta);

		dir = x_length * x_dir + y_length * y_dir + z_length * dir;
	}

	inline void SpecularWeightedHemisphereSampling(Vec3f & dir, float glossiness)
	{
		Vec3f randomvec = Vec3f((static_cast<float>(rand()) / (RAND_MAX)), (static_cast<float>(rand()) / (RAND_MAX)), (static_cast<float>(rand()) / (RAND_MAX)));
		Vec3f x_dir = dir.Cross(randomvec);
		Vec3f y_dir = dir.Cross(x_dir);
		x_dir.Normalize();
		y_dir.Normalize();

		float theta = (static_cast<float>(rand()) / (RAND_MAX));
		float phy = (static_cast<float>(rand()) / (RAND_MAX));

		float x_length = cosf(phy) * pow(theta, 1 / (glossiness + 1));
		float y_length = sinf(phy) * pow(theta, 1 / (glossiness + 1));
		float z_length = sqrtf(1 - pow(theta, 1 / (glossiness + 1)) * pow(theta, 1 / (glossiness + 1)));

		dir = x_length * x_dir + y_length * y_dir + z_length * dir;
	}
}

