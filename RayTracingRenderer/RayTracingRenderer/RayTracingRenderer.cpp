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

void RayTraversing(Node * traversingnode, Node * node, Ray ray, Color24 & pixel, float & zbuffer, Ray originalray, HitInfo & hit) {

	int numberofchild = traversingnode->GetNumChild();
	HitInfo hitinfo = HitInfo();
	for (int i = 0; i < numberofchild; i++) {
		node = traversingnode->GetChild(i);
		Ray changedray = node->ToNodeCoords(ray);

		//hitinfo.duvw[0] = node->VectorTransformTo(hit.duvw[0]);
		//hitinfo.duvw[1] = node->VectorTransformTo(hit.duvw[1]);

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
			RayTraversing(node, childnode, changedray, pixel, zbuffer, originalray, hit);
			if (hit.node != nullptr && hit.node != hitinfo.node)
			{
				node->FromNodeCoords(hit);
				hitinfo = hit;
			}
			delete childnode;
		}
	}


	if (node->GetNodeObj() != nullptr)
	{
		if (hitinfo.node != nullptr)
		{
			// Get Z buffer
			zbuffer = hitinfo.z;
			//Shading
			if (materials.Find(node->GetMaterial()->GetName()) != nullptr)
			{
				pixel = (Color24)materials.Find(hitinfo.node->GetMaterial()->GetName())->Shade(originalray, hitinfo, lights, TIMEOFREFRECTION);
			}
			else
			{
				pixel = (Color24)environment.SampleEnvironment(originalray.dir);
			}
		}
	}
}

void BeginRender() {

	time_t time0;   // create timers.
	time_t time1;

	time(&time0);   // get current time.


	Node * node = new Node();
	Node * startnode = &rootNode;

	float * zbuffers = renderImage.GetZBuffer();
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

	for (int i = 0; i < renderImage.GetHeight(); i++) {
		for (int j = 0; j < renderImage.GetWidth(); j++) {
			cameraray[i * renderImage.GetWidth() + j].dir = f + (j + 0.5f) * (w / W)*x - (i + 0.5f) * (h / H)*y - camera.pos;
			cameraray[i * renderImage.GetWidth() + j].p = camera.pos;
			zbuffers[i * renderImage.GetWidth() + j] = BIGFLOAT;
			cameraray[i * renderImage.GetWidth() + j].Normalize();
		}
	}
#pragma omp parallel for
	for (int i = 0; i < renderImage.GetHeight(); i++) {
		for (int j = 0; j < renderImage.GetWidth(); j++) {
			HitInfo hit = HitInfo();
			hit.duvw[0] = (w / W) * x;
			hit.duvw[1] = (h / H) * y;
			RayTraversing(startnode, node, cameraray[i * renderImage.GetWidth() + j], pixels[i * renderImage.GetWidth() + j], zbuffers[i * renderImage.GetWidth() + j], cameraray[i * renderImage.GetWidth() + j], hit);
		}
	}

	time(&time1);
	double seconds = time1 - time0;

	printf("Done. Time was %f \n", seconds);
	renderImage.SaveImage("saveimage.png");
	delete node;
	delete cameraray;
	return;
}

void StopRender() {
}