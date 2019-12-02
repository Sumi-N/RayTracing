#pragma once

#include <cyPhotonMap.h>

extern cyPhotonMap photonMap;
extern LightList lights;


inline void SetupPhotonMap()
{
	photonMap.Clear();

	for (auto light = lights.begin(); light != lights.end(); ++light)
	{
		if ((*light)->IsPhotonSource())
		{
			PointLight * point = (PointLight *)*light;
			photonMap.AddPhoton(point->GetLightPosition(), Vec3f(0,0,0), (*light)->GetPhotonIntensity());
		}
	}

	FILE *fp = fopen("..//PhotonMapViz//photonmap.dat", "wb");
	fwrite(photonMap.GetPhotons(), sizeof(cyPhotonMap::Photon), photonMap.NumPhotons(), fp);
	fclose(fp);
	return;
}