#pragma once

#include "constant.h"
#include "scene.h"
#include "utility.h"

extern Camera camera;
extern Node * rootNode;

namespace CameraMode 
{
	using namespace cy;
	using namespace Utility;

	Color BlurMode(Ray ray)
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

			pixelcolors[i] = RayTraversing(&rootNode, cameraraies[i], hits[i], BOUNCINGTIME);
			returnColor += pixelcolors[i];
		}

		return returnColor / RAYPERPIXELFORBLUREFFECT;
	}

	Color AdaptiveSampling()
	{

	}
}