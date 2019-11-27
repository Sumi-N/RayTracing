#pragma once
//-------------------------------------------------------------------------------
///
/// \file       xmlload.cpp 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    11.0
/// \date       August 21, 2019
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#include "scene.h"
#include "objects.h"
#include "materials.h"
#include "lights.h"
#include "texture.h"
#include "tinyxml/tinyxml.h"

//-------------------------------------------------------------------------------

extern Node rootNode;
extern Camera camera;
extern RenderImage renderImage;
extern MaterialList materials;
extern LightList lights;
extern ObjFileList objList;
extern TexturedColor background;
extern TexturedColor environment;
extern TextureList textureList;

//-------------------------------------------------------------------------------

#ifdef WIN32
#define COMPARE(a,b) (_stricmp(a,b)==0)
#else
#define COMPARE(a,b) (strcasecmp(a,b)==0)
#endif

//-------------------------------------------------------------------------------

void LoadScene(TiXmlElement *element);
void LoadNode(Node *node, TiXmlElement *element, int level = 0);
void LoadTransform(Transformation *trans, TiXmlElement *element, int level);
void LoadMaterial(TiXmlElement *element);
void LoadLight(TiXmlElement *element);
void ReadVector(TiXmlElement *element, Vec3f &v);
void ReadColor(TiXmlElement *element, Color  &c);
void ReadFloat(TiXmlElement *element, float  &f, char const *name = "value");
TextureMap* ReadTexture(TiXmlElement *element);
Texture* ReadTexture(char const *filename);

//-------------------------------------------------------------------------------

struct NodeMtl
{
	Node *node;
	char const *mtlName;
};

//-------------------------------------------------------------------------------

extern std::vector<NodeMtl> nodeMtlList;

int LoadScene(char const *filename);