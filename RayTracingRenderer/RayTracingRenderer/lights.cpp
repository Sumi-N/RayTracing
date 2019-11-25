#include "lights.h"
#include "constant.h"
#include "utility.h"

using namespace Utility;
extern Node rootNode;

bool ShadowTraversing(Node * traversingnode, Node * node, Ray ray, float t_max)
{
	int numberofchild = traversingnode->GetNumChild();
	Ray currentray = ray;
	for (int i = 0; i < numberofchild; i++)
	{
		node = traversingnode->GetChild(i);
		if (node->GetNodeObj() != nullptr)
		{
			Ray changedray = node->ToNodeCoords(currentray);
			HitInfo fake; fake.z = t_max;
			if (node->GetNodeObj()->IntersectRay(changedray, fake, 1))
				return true;
		}

		if (node != nullptr)
		{
			Node childnode;
			if (ShadowTraversing(node, &childnode, node->ToNodeCoords(currentray), t_max))
			{
				return true;
			}
		}
	}
	return false;
}

float GenLight::Shadow(Ray ray, float t_max)
{
	Node node;
	Node * startnode = &rootNode;
	if (ShadowTraversing(startnode, &node, ray, t_max))
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

	direction = position - p;
	float length = direction.Length();
	ConeUniformSampling(direction, size / 2);
	direction.Normalize();
	if (Shadow(Ray(p, direction), length) == 0.0f)
	{
		ratio = 0.0f;
	}
	else
	{
		ratio = 1.0f;
	}

	return ratio * intensity;
}