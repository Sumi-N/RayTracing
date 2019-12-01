#pragma once
#define _USE_MATH_DEFINES

#include <math.h>
#include <stdio.h> 
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

	inline float GetUniformRamdomFloat()
	{
		return (static_cast<float>(rand()) / (RAND_MAX));
	}

	inline Vec3f HemisphereUniformSampling(const Vec3f & z_dir, float theta, float phy)
	{
		Vec3f x_dir;
		if (z_dir.y != 0 && z_dir.x != 0)
		{
			x_dir = Vec3f(-1 * z_dir.y, z_dir.x, 0);
		}
		else
		{
			x_dir = Vec3f(0, -1 * z_dir.z, z_dir.y);
		}
		Vec3f y_dir = z_dir.Cross(x_dir);
		x_dir.Normalize();
		y_dir.Normalize();

		float x_length = cosf(2 * static_cast<float>(M_PI) * phy) * sqrtf(1 - (1 - theta) * (1 - theta));
		float y_length = sinf(2 * static_cast<float>(M_PI) * phy) * sqrtf(1 - (1 - theta) * (1 - theta));
		float z_length = 1 - theta;

		return x_length * x_dir + y_length * y_dir + z_length * z_dir;
	}

	inline Vec3f CosineWeightedHemisphereUniformSampling(const Vec3f & z_dir, float xi1, float xi2)
	{
		Vec3f x_dir;
		if (z_dir.y != 0 && z_dir.x != 0)
		{
			x_dir = Vec3f(-1 * z_dir.y, z_dir.x, 0);
		}
		else
		{
			x_dir = Vec3f(0, -1 * z_dir.z, z_dir.y);
		}
		Vec3f y_dir = z_dir.Cross(x_dir);
		x_dir.Normalize();
		y_dir.Normalize();

		float x_length = cosf(2 * static_cast<float>(M_PI) * xi2) * sqrtf(xi1);
		float y_length = sinf(2 * static_cast<float>(M_PI) * xi2) * sqrtf(xi1);
		float z_length = sqrtf(1 - xi1);

		return x_length * x_dir + y_length * y_dir + z_length * z_dir;
	}

	inline Vec3f SpecularWeightedHemisphereSampling(const Vec3f & z_dir, const float & alpha, float xi1, float xi2)
	{
		Vec3f x_dir;
		if (z_dir.y != 0 && z_dir.x != 0)
		{
			x_dir = Vec3f(-1 * z_dir.y, z_dir.x, 0);
		}
		else
		{
			x_dir = Vec3f(0, -1 * z_dir.z, z_dir.y);
		}
		Vec3f y_dir = z_dir.Cross(x_dir);
		x_dir.Normalize();
		y_dir.Normalize();

		float x_length = cosf(2 * static_cast<float>(M_PI) * xi2) * sqrtf(1 - powf(xi1, 1 / (1 + alpha) * powf(xi1, 1 / (1 + alpha))));
		float y_length = sinf(2 * static_cast<float>(M_PI) * xi2) * sqrtf(1 - powf(xi1, 1 / (1 + alpha) * powf(xi1, 1 / (1 + alpha))));
		float z_length = powf(xi1, 1 / (1 + alpha));

		return x_length * x_dir + y_length * y_dir + z_length * z_dir;
	}

	inline void ClampColorValue(cy::Color & color)
	{
		if (color.r > 1)
		{
			color.r = 1.0f;
		}
		if (color.g > 1)
		{
			color.g = 1.0f;
		}
		if (color.b > 1)
		{
			color.b = 1.0f;
		}
		return;
	}

	inline float AverageRGBValue(const cy::Color & color) 
	{
		return (color.r + color.g + color.b) / 3;
	}

	inline float TotalRGBValue(const cy::Color & color) 
	{
		return (color.r + color.g + color.b);
	}
}

