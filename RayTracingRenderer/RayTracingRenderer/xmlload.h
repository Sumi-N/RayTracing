#pragma once

#include "scene.h"
#include "objects.h"
#include "tinyxml/tinyxml.h"

extern Node rootNode;
extern Camera camera;
extern RenderImage renderImage;

int LoadScene(char const *filename);
void LoadScene(TiXmlElement *element);
void LoadNode(Node *node, TiXmlElement *element, int level = 0);
void LoadTransform(Transformation *trans, TiXmlElement *element, int level);
void ReadVector(TiXmlElement *element, Vec3f &v);
void ReadFloat(TiXmlElement *element, float  &f, char const *name = "value");