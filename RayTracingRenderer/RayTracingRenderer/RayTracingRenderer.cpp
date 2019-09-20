// RayTracingRenderer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <iostream>
#include <scene.h>
#include <xmlload.h>
#include <viewport.h>
#include <objects.h>

Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
Plane thePlane;
ItemFileList<Object> objList;
LightList lights;
MaterialList materials;
std::vector<NodeMtl> nodeMtlList;


int main()
{
	//LoadScene(".\\objfiles\\plane.obj");
	//LoadScene(".\\objfiles\\teapot.obj");
	//LoadScene(".\\objfiles\\teapot-low.obj");
	LoadScene(".\\xmlfiles\\assignment5.xml");
	//LoadScene(".\\xmlfiles\\SimpleTestScene.xml");
	//printf("%d", rootNode.GetNumChild());
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

bool CheckZbuffer(float & zbuffer, float answer) {
	if (answer < zbuffer) {
		return true;
	}
	return false;
}

void UpdateHitInfo(Ray ray, float & zbuffer, HitInfo & hitinfo, Node * node) {
	float a = ray.dir.Dot(ray.dir);
	float b = 2 * ray.dir.Dot(ray.p);
	float c = ray.p.Dot(ray.p) - 1;

	if (b*b - 4 * a*c >= 0) {

		float answer1 = (-1 * b + sqrt(b*b - 4 * a*c)) / (2 * a);
		float answer2 = (-1 * b - sqrt(b*b - 4 * a*c)) / (2 * a);

		float large;
		float small;

		if (answer1 >= answer2) {
			large = answer1;
			small = answer2;
		}
		else
		{
			large = answer2;
			small = answer1;
		}

		if (small < 0) {
			if (large > 0) {
				if (CheckZbuffer(hitinfo.z, large)) {
					zbuffer = large;
					hitinfo.z = large;
					hitinfo.front = false;
					hitinfo.p = ray.p + large * ray.dir;
					hitinfo.N = hitinfo.p;
					hitinfo.node = node;
					node->FromNodeCoords(hitinfo);
				}
			}
		}
		else
		{
			if (CheckZbuffer(hitinfo.z, small)) {
				zbuffer = small;
				hitinfo.z = small;
				hitinfo.front = true;
				hitinfo.p = ray.p + small * ray.dir;
				hitinfo.N = hitinfo.p;
				hitinfo.node = node;
				node->FromNodeCoords(hitinfo);
			}
		}

	}
}

void ConvertRayCoordination(Node * traversingnode, Node * node, Ray ray, Color24 & pixel, float & zbuffer, Ray originalray, HitInfo & hit) {

	int numberofchild = traversingnode->GetNumChild();
	HitInfo hitinfo = HitInfo();
	for (int i = 0; i < numberofchild; i++) {
		node = traversingnode->GetChild(i);
		Ray changedray = node->ToNodeCoords(ray);

		if (node->GetNodeObj() != nullptr) {
			UpdateHitInfo(changedray, zbuffer, hitinfo, node);
			hit = hitinfo;
		}

		if (node != nullptr) {
			Node * childnode = new Node();
			ConvertRayCoordination(node, childnode, changedray, pixel, zbuffer, originalray, hit);
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
				pixel = (Color24)materials.Find(hitinfo.node->GetMaterial()->GetName())->Shade(originalray, hitinfo, lights,5);
			}
		}
	}
}

void BeginRender() {

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

	//test
	//ConvertRayCordinationTest(startnode, node, cameraray[0], pixels[0], zbuffers[0]);
	for (int i = 0; i < renderImage.GetHeight(); i++) {
		for (int j = 0; j < renderImage.GetWidth(); j++) {
			HitInfo hit = HitInfo();
			ConvertRayCoordination(startnode, node, cameraray[i * renderImage.GetWidth() + j], pixels[i * renderImage.GetWidth() + j], zbuffers[i * renderImage.GetWidth() + j], cameraray[i * renderImage.GetWidth() + j], hit);
		}
	}

	printf("Ready \n");
	renderImage.SaveImage("saveimage.png");
	//renderImage.SaveZImage("savezimage.png");
	delete node;
	delete cameraray;
	return;
}

void StopRender() {
}