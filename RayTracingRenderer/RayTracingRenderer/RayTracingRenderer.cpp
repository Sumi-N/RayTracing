// RayTracingRenderer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <iostream>
#include <scene.h>
#include <xmlload.h>
#include <viewport.h>
#include <objects.h>
#include <omp.h>
#include <cyBVH.h>
#include <time.h>

Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
Plane thePlane;
ItemFileList<Object> objList;
LightList lights;
MaterialList materials;
std::vector<NodeMtl> nodeMtlList;
TexturedColor background;
TexturedColor environment;
TextureList textureList;

#define TIMEOFREFRECTION 5
#define RAYPERPIXEL 4
#define SHADOWBIAS 0.0005f
#define MAXSAMPLECOUNT 32
#define SAMPLEVARIENCE 0.005f


int main()
{
	//LoadScene(".\\xmlfiles\\assignment5.xml");
	LoadScene(".\\xmlfiles\\assignment6.xml");
	ShowViewport();
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

Color RayTraversing(Node * traversingnode, Node * node, Ray ray, float & zbuffer, Ray originalray, HitInfo & hit) {

	int numberofchild = traversingnode->GetNumChild();
	HitInfo hitinfo = HitInfo();

	//hitinfo.duvw[0] = node->VectorTransformTo(hit.duvw[0]);
	//hitinfo.duvw[1] = node->VectorTransformTo(hit.duvw[1]);

	for (int i = 0; i < numberofchild; i++) {
		node = traversingnode->GetChild(i);
		Ray changedray = node->ToNodeCoords(ray);

		if (node->GetNodeObj() != nullptr) {
			if(node->GetNodeObj()->IntersectRay(changedray, hitinfo, 0))
			{
				hitinfo.node = node;
				node->FromNodeCoords(hitinfo);
			}
			hit = hitinfo;
		}

		if (node != nullptr) {
			Node * childnode = new Node();
			RayTraversing(node, childnode, changedray, zbuffer, originalray, hit);
			if (hit.node != nullptr && hit.node != hitinfo.node)
			{
				node->FromNodeCoords(hit);
				hitinfo = hit;
			}
			delete childnode;
		}
	}

	//hitinfo.duvw[0] = node->VectorTransformFrom(hitinfo.duvw[0]);
	//hitinfo.duvw[1] = node->VectorTransformFrom(hitinfo.duvw[1]);

	if (node->GetNodeObj() != nullptr)
	{
		if (hitinfo.node != nullptr)
		{
			// Get Z buffer
			zbuffer = hitinfo.z;
			//Shading
			if (materials.Find(node->GetMaterial()->GetName()) != nullptr)
			{
				//pixel = (Color24)materials.Find(hitinfo.node->GetMaterial()->GetName())->Shade(originalray, hitinfo, lights, TIMEOFREFRECTION);
				return materials.Find(hitinfo.node->GetMaterial()->GetName())->Shade(originalray, hitinfo, lights, TIMEOFREFRECTION);
			}
		}
		else
		{
			return environment.SampleEnvironment(originalray.dir);
		}
	}
}

Color AdaptiveSampling(Ray ray, float radiusrate, Vec3f xaxis, Vec3f yaxis, Node * traversingnode, Node * node, float & zbuffer, Ray originalray, uint8_t & samplecount)
{
	samplecount++;
	Ray cameraraies[RAYPERPIXEL];
	HitInfo hits[RAYPERPIXEL];
	
	// Assuming RAYPERPIXEL is 4
	cameraraies[0].dir = ray.dir + radiusrate * xaxis + radiusrate * yaxis;
	cameraraies[1].dir = ray.dir + radiusrate * xaxis - radiusrate * yaxis;
	cameraraies[2].dir = ray.dir - radiusrate * xaxis + radiusrate * yaxis;
	cameraraies[3].dir = ray.dir - radiusrate * xaxis - radiusrate * yaxis;

	Color averagepixelcolor = Color(0, 0, 0);
	Color pixelcolors[RAYPERPIXEL];

	for (int i = 0; i < RAYPERPIXEL; i++)
	{
		cameraraies[i].p = ray.p;
		cameraraies[i].Normalize();

		pixelcolors[i] = RayTraversing(traversingnode, node, cameraraies[i], zbuffer, cameraraies[i], hits[i]);
		averagepixelcolor += RayTraversing(traversingnode, node, cameraraies[i], zbuffer, cameraraies[i], hits[i]);
	}

	averagepixelcolor /= RAYPERPIXEL;
	Color variance = Color(0, 0, 0);
	
	for (int i = 0; i < RAYPERPIXEL; i++)
	{
		variance += ((pixelcolors[i] - averagepixelcolor) * (pixelcolors[i] - averagepixelcolor));
	}
	variance /= RAYPERPIXEL;

	Color answercolor = Color(0, 0, 0);

	for (int i = 0; i < RAYPERPIXEL; i++)
	{
		if (abs(variance.r) >= SAMPLEVARIENCE || abs(variance.g) >= SAMPLEVARIENCE || abs(variance.b) >= SAMPLEVARIENCE)
		{
			if (samplecount < MAXSAMPLECOUNT)
			{
				answercolor += AdaptiveSampling(cameraraies[i], radiusrate / 2.0f, xaxis, yaxis, traversingnode, node, zbuffer, originalray, samplecount);
			}
			else
			{
				answercolor += pixelcolors[i];
			}
		}
		else
		{
			answercolor += pixelcolors[i];
		}
	}

	return answercolor / RAYPERPIXEL;
}

void BeginRender() {

	time_t time0;   // create timers.
	time_t time1;

	time(&time0);   // get current time.


	Node * node = new Node();
	Node * startnode = &rootNode;

	float * zbuffers = renderImage.GetZBuffer();
	uint8_t * samplecount = renderImage.GetSampleCount();
	Color24* pixels = renderImage.GetPixels();
	Ray * cameraray = new Ray[renderImage.GetHeight() * renderImage.GetWidth()];

	float l = 1.0f;
	float h = 2 * l * tanf((camera.fov / 2) * 3.14f / 180);
	float w = camera.imgWidth * (h / camera.imgHeight);

	int H = camera.imgHeight;
	int W = camera.imgWidth;

	Vec3f x = camera.dir.Cross(camera.up);
	Vec3f y = camera.up;

	Vec3f f = camera.pos + l * camera.dir + (h / 2) * y - (w / 2) * x;

	for (int i = 0; i < renderImage.GetHeight(); i++) 
	{
		for (int j = 0; j < renderImage.GetWidth(); j++) 
		{
			cameraray[i * renderImage.GetWidth() + j].dir = f + (j + 0.5f) * (w / W)*x - (i + 0.5f) * (h / H)*y - camera.pos;
			cameraray[i * renderImage.GetWidth() + j].p = camera.pos;
			zbuffers[i * renderImage.GetWidth() + j] = BIGFLOAT;
			samplecount[i * renderImage.GetWidth() + j] = 0;
			cameraray[i * renderImage.GetWidth() + j].Normalize();
		}
	}
#pragma omp parallel for
	for (int i = 0; i < renderImage.GetHeight(); i++) {
		for (int j = 0; j < renderImage.GetWidth(); j++) {

			//hit.duvw[0] = (w / W) * x;
			//hit.duvw[1] = (h / H) * y;

			pixels[i * renderImage.GetWidth() + j] = (Color24)AdaptiveSampling(cameraray[i * renderImage.GetWidth() + j], 0.25f, (w / W)*x, (h / H)*y, startnode, node, zbuffers[i * renderImage.GetWidth() + j], cameraray[i * renderImage.GetWidth() + j], samplecount[i * renderImage.GetWidth() + j]);


			//HitInfo hit = HitInfo();
			//if (i == 125 && j == 22)
			//{
			//	RayTraversing(startnode, node, cameraray[i * renderImage.GetWidth() + j], pixels[i * renderImage.GetWidth() + j], zbuffers[i * renderImage.GetWidth() + j], cameraray[i * renderImage.GetWidth() + j], hit);
			//}

			//if (i == 432 && j == 6)
			//{
			//	RayTraversing(startnode, node, cameraray[i * renderImage.GetWidth() + j], pixels[i * renderImage.GetWidth() + j], zbuffers[i * renderImage.GetWidth() + j], cameraray[i * renderImage.GetWidth() + j], hit);
			//}
			//pixels[i * renderImage.GetWidth() + j] = (Color24)RayTraversing(startnode, node, cameraray[i * renderImage.GetWidth() + j],  zbuffers[i * renderImage.GetWidth() + j], cameraray[i * renderImage.GetWidth() + j], hit);
		}
	}

	time(&time1);
	double seconds = (double)(time1 - time0);

	printf("Done. Time was %f \n", seconds);
	renderImage.SaveImage("saveimage.png");
	delete node;
	delete cameraray;

	return;
}

void StopRender() {
}



Color FindReflectionAndRefraction(Node * traversingnode, Node * node, Ray originalray, Ray ray, HitInfo & hit, int bounce)
{
	int numberofchild = traversingnode->GetNumChild();
	HitInfo hitinfo = HitInfo();
	for (int i = 0; i < numberofchild; i++)
	{
		node = traversingnode->GetChild(i);
		Ray changedray = node->ToNodeCoords(ray);

		if (node->GetNodeObj() != nullptr)
			if (node->GetNodeObj()->IntersectRay(changedray, hitinfo, 0))
			{
				hitinfo.node = node;
				node->FromNodeCoords(hitinfo);
			}
		hit = hitinfo;

		if (node != nullptr)
		{
			Node * childnode = new Node();
			FindReflectionAndRefraction(node, childnode, originalray, changedray, hit, bounce);
			if (hit.node != nullptr && hit.node != hitinfo.node)
			{
				node->FromNodeCoords(hit);
				hitinfo = hit;
			}
			delete childnode;
		}
	}

	//Shading
	if (node->GetNodeObj() != nullptr)
	{
		if (hitinfo.node != nullptr)
		{
			if (materials.Find(node->GetMaterial()->GetName()) != nullptr)
			{
				return materials.Find(hitinfo.node->GetMaterial()->GetName())->Shade(originalray, hitinfo, lights, bounce);
			}
		}
		return environment.SampleEnvironment(originalray.dir);
	}
	else
	{
		return environment.SampleEnvironment(originalray.dir);
	}
}

float FresnelReflections(const HitInfo & hInfo, float refractionIndex, float cos1)
{
	float R0;
	if (hInfo.front)
	{
		R0 = (1 - refractionIndex) / (1 + refractionIndex) * (1 - refractionIndex) / (1 + refractionIndex);
	}
	else
	{
		R0 = (refractionIndex - 1) / (1 + refractionIndex) * (refractionIndex - 1) / (1 + refractionIndex);
	}

	return R0 + (1 - R0) * pow(1 - cos1, 5);
}

Color Reflection(Ray const & ray, const HitInfo & hInfo, int bounce)
{
	// bounce 1 is reflect 1 time
	if (bounce <= 0)
	{
		return Color(0, 0, 0);
	}
	else
	{
		Vec3f P = hInfo.N;
		P.Normalize(); P = SHADOWBIAS * P; P += hInfo.p;

		Vec3f V = -1 * ray.dir;
		Vec3f R = 2 * (hInfo.N.Dot(V)) * hInfo.N - V;
		R.Normalize();

		// S is a starting point from the surface point
		Ray S;
		S.dir = R;
		S.p = P;

		// Start traversing node 
		Node * node = new Node();
		Node * startnode = &rootNode;
		HitInfo hitinfo;

		Color returnColor = FindReflectionAndRefraction(startnode, node, S, S, hitinfo, bounce - 1);
		delete node;
		return returnColor;
	}
}

Color Refraction(Ray const & ray, const HitInfo & hInfo, int bounce, float refractionIndex, Color refraction)
{
	float R;
	// bounce 1 is refract 1 time
	if (bounce <= -1)
	{
		return Color(0, 0, 0);
	}
	else
	{
		Vec3f P = hInfo.N;
		P.Normalize();

		Vec3f V = -1 * ray.dir;
		V.Normalize();

		Vec3f N = hInfo.N;
		N.Normalize();

		float cos1, cos2, sin1, sin2;

		if (hInfo.front)
		{
			P = -1 * SHADOWBIAS * P;  P += hInfo.p;
			cos1 = V.Dot(N); sin1 = sqrt(1 - (cos1 * cos1));
		}
		else
		{
			P = SHADOWBIAS * P;  P += hInfo.p;
			cos1 = V.Dot(-N); sin1 = sqrt(1 - (cos1 * cos1));
		}

		// Calculate fresnel reflection
		R = FresnelReflections(hInfo, refractionIndex, cos1);

		Vec3f T_h, T_v, T;

		if (hInfo.front)
		{
			sin2 = (1 / refractionIndex) * sin1;
		}
		else
		{
			sin2 = refractionIndex * sin1;
		}
		cos2 = sqrt(1 - (sin2 * sin2));

		// Total internal reflection
		if (sin2 > 1)
		{
			return Color(0, 0, 0);
		}

		if (hInfo.front)
		{
			// Horizontal dirction Vector
			T_h = -cos2 * N;
			// Vertical Direction Vector
			T_v = (V - (V.Dot(N))* N);
		}
		else
		{
			// Horizontal dirction Vector
			T_h = -cos2 * -N;
			// Vertical Direction Vector
			T_v = (V - (V.Dot(-N))* -N);
		}
		T_v.Normalize();
		T_v = -sin2 * T_v;
		// Combined horizontal and vertical
		T = T_h + T_v;
		float s = T.Length();

		// S is a starting point from the surface point
		Ray S;
		S.dir = T;
		S.p = P;

		// Start traversing node 
		Node * node = new Node();
		Node * startnode = &rootNode;
		HitInfo hitinfo;

		// Refraction part
		Color returnColor = (1 - R) * refraction * FindReflectionAndRefraction(startnode, node, S, S, hitinfo, bounce - 1);
		// Reflection part
		returnColor += R * refraction * Reflection(ray, hInfo, bounce);
		delete node;
		return returnColor;
	}
}

Color CalculateAbsorption(Color in, Color absorption, float distance)
{
	float rout = exp(-1 * absorption.r * distance) * in.r;
	float gout = exp(-1 * absorption.g * distance) * in.g;
	float bout = exp(-1 * absorption.b * distance) * in.b;
	return Color(rout, gout, bout);
}

Color MtlBlinn::Shade(Ray const & ray, const HitInfo & hInfo, const LightList & lights, int bounce) const
{
	Vec3f N = hInfo.N;
	Color color = Color();

	for (auto light = lights.begin(); light != lights.end(); ++light)
	{
		if ((*light)->IsAmbient())
		{
			color += (*light)->Illuminate(hInfo.p, N) * this->diffuse.Sample(hInfo.uvw, hInfo.duvw);
		}
		else if (strcmp((*light)->GetName(), "directLight") == 0)
		{
			Vec3f L = -1 * (*light)->Direction(hInfo.p);
			L.Normalize();
			Vec3f V = ray.p - hInfo.p;
			V.Normalize();
			Vec3f H = (V + L) / (V + L).Length();
			H.Normalize();

			// Incoming light
			Color IR;
			IR.r = 0; IR.g = 0; IR.b = 0;
			if (N.Dot(-1 * (*light)->Direction(hInfo.p)) >= 0)
			{
				IR = (*light)->Illuminate(hInfo.p, N) * N.Dot(L);
			}

			float oneofcos = 1 / hInfo.N.Dot(L);
			Color specularpart = oneofcos * pow(H.Dot(hInfo.N), this->glossiness) * this->specular.Sample(hInfo.uvw, hInfo.duvw);
			Color diffusepart = this->diffuse.Sample(hInfo.uvw, hInfo.duvw);
			color += (diffusepart + specularpart) * IR;
		}
		else if (strcmp((*light)->GetName(), "pointLight") == 0)
		{
			Vec3f L = -1 * (*light)->Direction(hInfo.p);
			L.Normalize();
			Vec3f V = ray.p - hInfo.p;
			V.Normalize();
			Vec3f H = (V + L) / (V + L).Length();
			H.Normalize();

			// Incoming light
			Color IR;
			IR.r = 0; IR.g = 0; IR.b = 0;
			if (N.Dot(-1 * (*light)->Direction(hInfo.p)) >= 0)
			{
				IR = (*light)->Illuminate(hInfo.p, N) * N.Dot(L);
			}

			float oneofcos = 1 / hInfo.N.Dot(L);
			Color specularpart = oneofcos * pow(H.Dot(hInfo.N), this->glossiness) * this->specular.Sample(hInfo.uvw, hInfo.duvw);
			Color diffusepart = this->diffuse.Sample(hInfo.uvw, hInfo.duvw);
			color += (diffusepart + specularpart) * IR;
		}
	}

	// Caluculate only relfection part for reflection
	if (this->reflection.Sample(hInfo.uvw) != Color(0, 0, 0))
	{
		color += this->reflection.Sample(hInfo.uvw) * Reflection(ray, hInfo, bounce);
	}

	// Caluculate refraction part
	if (this->refraction.Sample(hInfo.uvw) != Color(0, 0, 0))
	{
		// When it is a back side hit, it means that absorption gonna happen during inside the material the light go through
		if (!hInfo.front)
		{
			color += CalculateAbsorption(color, absorption, hInfo.z);
		}
		color += Refraction(ray, hInfo, bounce, ior, refraction.Sample(hInfo.uvw));
	}

	return color;
}

bool DetectShadow(Node * traversingnode, Node * node, Ray ray, float t_max)
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
			Node * childnode = new Node();
			DetectShadow(node, childnode, node->ToNodeCoords(currentray), t_max);
			delete childnode;
		}
	}
	return false;
}

float GenLight::Shadow(Ray ray, float t_max)
{
	Node * node = new Node();
	Node * startnode = &rootNode;
	if (DetectShadow(startnode, node, ray, t_max))
	{
		delete node;
		return 0.0f;
	}
	else
	{
		delete node;
		return 1.0f;
	}
}

bool CheckZbuffer(float & zbuffer, float answer)
{
	if (answer < zbuffer)
	{
		return true;
	}
	return false;
}

bool Sphere::IntersectRay(Ray const & ray, HitInfo & hInfo, int hitSide) const
{
	float a = ray.dir.Dot(ray.dir);
	float b = 2 * ray.dir.Dot(ray.p);
	float c = ray.p.Dot(ray.p) - 1;

	if (b*b - 4 * a*c >= 0)
	{

		float answer1 = (-1 * b + sqrt(b*b - 4 * a*c)) / (2 * a);
		float answer2 = (-1 * b - sqrt(b*b - 4 * a*c)) / (2 * a);

		float large;
		float small;

		float answer;

		if (answer1 >= answer2)
		{
			large = answer1;
			small = answer2;
		}
		else
		{
			large = answer2;
			small = answer1;
		}

		if (small < 0)
		{
			if (large > 0)
			{

				// this is for shadowprocess
				if (hitSide == 1)
				{
					answer = large;
					if ((answer > SHADOWBIAS) & (answer <= hInfo.z))
					{
						return true;
					}
				}
				else
				{
					if (CheckZbuffer(hInfo.z, large))
					{
						hInfo.z = large;
						hInfo.front = false;
						hInfo.p = ray.p + large * ray.dir;
						hInfo.N = hInfo.p;
						float u = (1 / 2 * 3.14f) * atan2f(hInfo.p.y, hInfo.p.x) + .5f;
						float v = (1 / 3.14f) * asinf(hInfo.p.z) + 0.5f;
						hInfo.uvw = Vec3f(u, v, 0.0f);
						return true;
					}
				}
			}
		}
		else
		{
			// this is for shadowprocess
			if (hitSide == 1)
			{
				answer = small;
				if ((answer > SHADOWBIAS) & (answer <= hInfo.z))
				{
					return true;
				}
			}
			else
			{
				if (CheckZbuffer(hInfo.z, small))
				{
					hInfo.z = small;
					hInfo.front = true;
					hInfo.p = ray.p + small * ray.dir;
					hInfo.N = hInfo.p;
					float u = (1 / (2 * 3.14f)) * atan2f(hInfo.p.y, hInfo.p.x) + .5f;
					float v = (1 / 3.14f) * asinf(hInfo.p.z) + 0.5f;
					hInfo.uvw = Vec3f(u, v, 0.0f);
					return true;
				}
			}
		}

	}
	return false;
}

bool Plane::IntersectRay(Ray const & ray, HitInfo & hInfo, int hitSide) const
{
	float t = -1 * (ray.p.Dot(Vec3f(0, 0, 1))) / ray.dir.Dot(Vec3f(0, 0, 1));
	Vec3f point = ray.p + t * ray.dir;
	if (t < 0)
		return false;
	if (point.x >= -1 && point.x <= 1)
	{
		if (point.y >= -1 && point.y <= 1)
		{
			// this is for shadowprocess
			if (hitSide == 1)
			{
				if ((t > SHADOWBIAS) & (t <= hInfo.z))
				{
					return true;
				}
				else
				{
					return false;
				}
			}

			if (CheckZbuffer(hInfo.z, t))
			{
				hInfo.front = true;
				hInfo.p = point;
				hInfo.z = t;
				hInfo.N = Vec3f(0, 0, 1);
				hInfo.uvw = Vec3f(0.5f * point.x + 0.5f, 0.5f * point.y + 0.5f, point.z);

				//float t2 = -1 * (ray.p + hInfo.duvw[0]).Dot(Vec3f(0, 0, 1)) / ray.dir.Dot(Vec3f(0, 0, 1));
				//Vec3f point2 = ray.p + hInfo.duvw[0] + t2 * (ray.dir);
				//hInfo.duvw[0] = 1.0f * (point2 - point);

				//float t3 = -1 * (ray.p + hInfo.duvw[1]).Dot(Vec3f(0, 0, 1)) / ray.dir.Dot(Vec3f(0, 0, 1));
				//Vec3f point3 = ray.p + hInfo.duvw[1] + t3 * (ray.dir);
				//hInfo.duvw[1] = 1.0f * (point3 - point);

				return true;
			}
		}
	}
	return false;
}

bool CheckBoxCollision(const float* vertices, Ray const & ray, float & answer)
{
	// The t value for each case
	float xmin, ymin, zmin;
	float xmax, ymax, zmax;
	float maxofmin, minofmax;
	bool detected = false;
	bool secondcheck = false;

	float x, y, z;

	xmin = (vertices[0] - ray.p.x) / ray.dir.x;
	y = xmin * ray.dir.y + ray.p.y;
	z = xmin * ray.dir.z + ray.p.z;
	if (y >= vertices[1] && y <= vertices[4] && z >= vertices[2] && z <= vertices[5])
	{
		if (detected)
		{
			maxofmin = xmin;
			secondcheck = true;
		}
		else
		{
			minofmax = xmin;
			detected = true;
		}
	}

	ymin = (vertices[1] - ray.p.y) / ray.dir.y;
	x = ymin * ray.dir.x + ray.p.x;
	z = ymin * ray.dir.z + ray.p.z;
	if (x >= vertices[0] && x <= vertices[3] && z >= vertices[2] && z <= vertices[5])
	{
		if (detected)
		{
			maxofmin = ymin;
			secondcheck = true;
		}
		else
		{
			minofmax = ymin;
			detected = true;
		}
	}

	zmin = (vertices[2] - ray.p.z) / ray.dir.z;
	x = zmin * ray.dir.x + ray.p.x;
	y = zmin * ray.dir.y + ray.p.y;
	if (x >= vertices[0] && x <= vertices[3] && y >= vertices[1] && y <= vertices[4])
	{
		if (detected)
		{
			maxofmin = zmin;
			secondcheck = true;
		}
		else
		{
			minofmax = zmin;
			detected = true;
		}
	}

	xmax = (vertices[3] - ray.p.x) / ray.dir.x;
	y = xmax * ray.dir.y + ray.p.y;
	z = xmax * ray.dir.z + ray.p.z;
	if (y >= vertices[1] && y <= vertices[4] && z >= vertices[2] && z <= vertices[5])
	{
		if (detected)
		{
			maxofmin = xmax;
			secondcheck = true;
		}
		else
		{
			minofmax = xmax;
			detected = true;
		}
	}

	ymax = (vertices[4] - ray.p.y) / ray.dir.y;
	x = ymax * ray.dir.x + ray.p.x;
	z = ymax * ray.dir.z + ray.p.z;
	if (x >= vertices[0] && x <= vertices[3] && z >= vertices[2] && z <= vertices[5])
	{
		if (detected)
		{
			maxofmin = ymax;
			secondcheck = true;
		}
		else
		{
			minofmax = ymax;
			detected = true;
		}
	}

	zmax = (vertices[5] - ray.p.z) / ray.dir.z;
	x = zmax * ray.dir.x + ray.p.x;
	y = zmax * ray.dir.y + ray.p.y;
	if (x >= vertices[0] && x <= vertices[3] && y >= vertices[1] && y <= vertices[4])
	{
		if (detected)
		{
			maxofmin = zmax;
			secondcheck = true;
		}
		else
		{
			minofmax = zmax;
			detected = true;
		}
	}

	if (!detected | !secondcheck)
	{
		return false;
	}

	if (maxofmin <= minofmax)
	{
		answer = maxofmin;
		return true;
	}
	else
	{
		answer = minofmax;
		return true;
	}
	return false;
}

bool TriObj::IntersectRay(Ray const & ray, HitInfo & hInfo, int hitSide) const
{
	const float * vertices = bvh.GetNodeBounds(bvh.GetRootNodeID());
	float answer = 0.0f;
	if (CheckBoxCollision(vertices, ray, answer))
	{
		return TraceBVHNode(ray, hInfo, hitSide, bvh.GetRootNodeID());
	}
	else
	{
		return false;
	}
}

bool TriObj::IntersectTriangle(Ray const & ray, HitInfo & hInfo, int hitSide, unsigned int faceID) const
{
	int i0 = F(faceID).v[0];
	int i1 = F(faceID).v[1];
	int i2 = F(faceID).v[2];

	Vec3f n = (v[i2] - v[i0]).Cross(v[i1] - v[i0]);
	float h = -1 * n.Dot(v[i0]);
	float t = -1 * (ray.p.Dot(n) + h) / ray.dir.Dot(n);

	// This is for shadow process
	if (hitSide == 1)
	{
		if (t - SHADOWBIAS < 0 || t > hInfo.z)
			return false;
	}

	if (t < 0)
		return false;

	Vec3f point = ray.p + t * ray.dir;

	Vec2f v0, v1, v2, x;
	if (std::abs(n.x) >= std::abs(n.y) && std::abs(n.x) >= std::abs(n.z))
	{
		v0 = Vec2f(v[i0].y, v[i0].z); v1 = Vec2f(v[i1].y, v[i1].z); v2 = Vec2f(v[i2].y, v[i2].z);
		x = Vec2f(point.y, point.z);
	}
	else if (std::abs(n.y) >= std::abs(n.x) && std::abs(n.y) >= std::abs(n.z))
	{
		v0 = Vec2f(v[i0].x, v[i0].z); v1 = Vec2f(v[i1].x, v[i1].z); v2 = Vec2f(v[i2].x, v[i2].z);
		x = Vec2f(point.x, point.z);
	}
	else if (std::abs(n.z) >= std::abs(n.x) && std::abs(n.z) >= abs(n.y))
	{
		v0 = Vec2f(v[i0].x, v[i0].y); v1 = Vec2f(v[i1].x, v[i1].y); v2 = Vec2f(v[i2].x, v[i2].y);
		x = Vec2f(point.x, point.y);
	}

	float a0, a1, a2;
	a0 = (v1 - x).Cross(v2 - x);
	a1 = (v2 - x).Cross(v0 - x);
	a2 = (v0 - x).Cross(v1 - x);

	if ((a0 >= 0 && a1 >= 0 && a2 >= 0) || (a0 < 0 && a1 < 0 && a2 < 0))
	{
		// This is for shadow process
		if (hitSide == 1)
		{
			return true;
		}

		if (CheckZbuffer(hInfo.z, t))
		{
			float a = (v1 - v0).Cross(v2 - v0);
			float beta0, beta1, beta2;
			beta0 = std::abs(a0 / a);
			beta1 = std::abs(a1 / a);
			beta2 = std::abs(a2 / a);

			Vec3f normal = beta0 * vn[i0] + beta1 * vn[i1] + beta2 * vn[i2];

			hInfo.front = true;
			hInfo.N = normal;
			hInfo.p = point;
			hInfo.z = t;
			// This is for texturing
			{
				Vec3f bc = Vec3f(beta0, beta1, beta2);
				hInfo.uvw = GetTexCoord(faceID, bc);
			}
			return true;
		}
	}
	return false;
}

bool TriObj::TraceBVHNode(Ray const & ray, HitInfo & hInfo, int hitSide, unsigned int nodeID) const
{
	if (bvh.IsLeafNode(nodeID))
	{
		unsigned int nodecount = bvh.GetNodeElementCount(nodeID);
		const unsigned int* nodeelementslist = bvh.GetNodeElements(nodeID);

		bool hit = false;
		for (int i = 0; i < (signed)nodecount; i++)
		{
			if (IntersectTriangle(ray, hInfo, hitSide, nodeelementslist[i]))
			{
				hit = true;
			}
		}

		return hit;
	}
	else
	{
		unsigned int child1id = bvh.GetFirstChildNode(nodeID);
		unsigned int child2id = bvh.GetSecondChildNode(nodeID);
		const float * vertices1 = bvh.GetNodeBounds(child1id);
		const float * vertices2 = bvh.GetNodeBounds(child2id);

		float t_child1, t_child2;

		if (!CheckBoxCollision(vertices1, ray, t_child1) & !(CheckBoxCollision(vertices2, ray, t_child2)))
		{
			return false;
		}

		if (t_child1 < t_child2)
		{
			if (t_child1 <= hInfo.z)
			{
				bool hitcheck = false;
				if (TraceBVHNode(ray, hInfo, hitSide, child1id))
				{
					hitcheck = true;

					if (t_child2 <= hInfo.z)
					{
						if (TraceBVHNode(ray, hInfo, hitSide, child2id))
						{
							hitcheck = true;
						}
					}
				}
				else if (t_child2 <= hInfo.z)
				{
					if (TraceBVHNode(ray, hInfo, hitSide, child2id))
					{
						hitcheck = true;
					}
				}
				return hitcheck;
			}
		}
		else if (t_child1 > t_child2)
		{
			if (t_child2 <= hInfo.z)
			{
				bool hitcheck = false;
				if (TraceBVHNode(ray, hInfo, hitSide, child2id))
				{
					hitcheck = true;

					if (t_child1 <= hInfo.z)
					{
						if (TraceBVHNode(ray, hInfo, hitSide, child1id))
						{
							hitcheck = true;
						}
					}
				}
				else if (t_child1 <= hInfo.z)
				{
					if (TraceBVHNode(ray, hInfo, hitSide, child1id))
					{
						hitcheck = true;
					}
				}
				return hitcheck;
			}
		}
		else if (t_child1 == t_child2)
		{
			bool hitcheck = false;
			if (TraceBVHNode(ray, hInfo, hitSide, child1id))
			{
				hitcheck = true;
			}
			if (TraceBVHNode(ray, hInfo, hitSide, child2id))
			{
				hitcheck = true;
			}
			return hitcheck;
		}
	}
	return false;
}