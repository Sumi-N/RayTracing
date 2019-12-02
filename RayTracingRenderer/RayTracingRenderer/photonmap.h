#pragma once

#include <cyPhotonMap.h>

extern Node rootNode;
extern cyPhotonMap photonMap;
extern LightList lights;

inline void PhotonTraversing(Node * node, Ray ray, HitInfo & hit)
{
	Node * currentnode;
	int numberofchild = node->GetNumChild();
	HitInfo tmphit = HitInfo();

	for (int i = 0; i < numberofchild; i++)
	{
		currentnode = node->GetChild(i);
		Ray currentray = currentnode->ToNodeCoords(ray);

		if (currentnode->GetNodeObj() != nullptr)
		{
			if (currentnode->GetNodeObj()->IntersectRay(currentray, tmphit, 0))
			{
				tmphit.node = currentnode;
				currentnode->FromNodeCoords(tmphit);
				hit = tmphit;
			}
		}

		if (currentnode != nullptr)
		{
			Node childnode;
			RayTraversing(currentnode, currentray, hit);
			if (hit.node != nullptr && hit.node != tmphit.node)
			{
				currentnode->FromNodeCoords(hit);
				tmphit = hit;
			}
		}
	}

	if (node == &rootNode)
	{
		hit = tmphit;
	}
}

inline void InitialPhotonRay(Ray ray, PointLight & light)
{
	Node * node = &rootNode;
	HitInfo hit = HitInfo();

	PhotonTraversing(node, ray, hit);

	if (hit.node != nullptr)
	{
		photonMap.AddPhoton(hit.p, ray.dir, light.GetPhotonIntensity());
	}
}

inline void SetupPhotonMap()
{
	photonMap.Clear();
	photonMap.Resize(NUMOFPHOTONS);

	for (auto light = lights.begin(); light != lights.end(); ++light)
	{
		if ((*light)->IsPhotonSource())
		{
			PointLight * point = (PointLight *)*light;

#pragma omp parallel for
			for (int i = 0; i < NUMOFPHOTONS; i++)
			{
				float x_dir = Utility::GetUniformRamdomFloat() - 0.5f;
				float y_dir = Utility::GetUniformRamdomFloat() - 0.5f;
				float z_dir = Utility::GetUniformRamdomFloat() - 0.5f;
				Vec3f ray_dir = Vec3f(x_dir, y_dir, z_dir);
				ray_dir.Normalize();
				Ray ray;
				ray.dir = ray_dir;
				ray.p = point->GetLightPosition();

				InitialPhotonRay(ray, *point);
			}
		}
	}

	FILE *fp = fopen("..//PhotonMapViz//photonmap.dat", "wb");
	fwrite(photonMap.GetPhotons(), sizeof(cyPhotonMap::Photon), photonMap.NumPhotons(), fp);
	fclose(fp);
	return;
}