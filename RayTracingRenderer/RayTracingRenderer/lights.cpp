#include "lights.h"
#include "constant.h"
#include "utility.h"

using namespace Utility;
extern Node rootNode;
bool DetectShadow(Node *, Node *, Ray, float);

float GenLight::Shadow(Ray ray, float t_max)
{
	Node node;
	Node * startnode = &rootNode;
	if (DetectShadow(startnode, &node, ray, t_max))
	{
		return 0.0f;
	}
	else
	{
		return 1.0f;
	}
}

Color PointLight::Illuminate(Vec3f const & p, Vec3f const & N) const
{
	Vec3f direction;
	float ratio = 0.0f;
	for (int i = 0; i < RAYPERPIXEL; i++)
	{
		direction = position - p;
		float length = direction.Length();
		ConeUniformSampling(direction, size / 2);
		direction.Normalize();
		if (Shadow(Ray(p, direction)) == 0.0f)
		{
			ratio = 0.0f;
			for (int j = 0; j < RAYPERPIXELFORSHADOW; j++)
			{
				direction = position - p;
				ConeUniformSampling(direction, size / 2);
				direction.Normalize();
				ratio += Shadow(Ray(p, direction), length);
			}
			ratio /= RAYPERPIXELFORSHADOW;
			break;
		}
		else
		{
			ratio = 1.0f;
		}
	}

	//for (int j = 0; j < RAYPERPIXELFORSHADOW; j++)
	//{
	//	direction = position - p;
	//	float length = direction.Length();
	//	ConeUniformSampling(direction, size / 2);
	//	direction.Normalize();
	//	ratio += Shadow(Ray(p, direction), length);
	//}
	//ratio /= RAYPERPIXELFORSHADOW;

	return ratio * intensity;
}