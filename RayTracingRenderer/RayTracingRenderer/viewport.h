//-------------------------------------------------------------------------------
///
/// \file       viewport.cpp 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    9.0
/// \date       August 21, 2019
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "scene.h"
#include "objects.h"
#include "lights.h"
#include "materials.h"
#include "texture.h"
#include <stdlib.h>
#include <time.h>

#ifdef USE_GLUT
# ifdef __APPLE__
#  include <GLUT/glut.h>
# else
#  include <GL/glut.h>
# endif
#else
# include <GL/freeglut.h>
#endif

//-------------------------------------------------------------------------------

void BeginRender(); // Called to start rendering (renderer must run in a separate thread)
void StopRender();  // Called to end rendering (if it is not already finished)

extern Node rootNode;
extern Camera camera;
extern RenderImage renderImage;
extern LightList lights;
extern TexturedColor background;

//-------------------------------------------------------------------------------

enum Mode {
	MODE_READY,         // Ready to render
	MODE_RENDERING,     // Rendering the image
	MODE_RENDER_DONE    // Rendering is finished
};

enum ViewMode
{
	VIEWMODE_OPENGL,
	VIEWMODE_IMAGE,
	VIEWMODE_Z,
	VIEWMODE_SAMPLECOUNT,
};

enum MouseMode {
	MOUSEMODE_NONE,
	MOUSEMODE_DEBUG,
	MOUSEMODE_ROTATE,
};

static Mode     mode = MODE_READY;       // Rendering mode
static ViewMode viewMode = VIEWMODE_OPENGL;  // Display mode
static MouseMode mouseMode = MOUSEMODE_NONE;   // Mouse mode
static int      startTime;                      // Start time of rendering
static int mouseX = 0, mouseY = 0;
static float viewAngle1 = 0, viewAngle2 = 0;
static GLuint viewTexture;

static int      dofDrawCount = 0;
static Color    *dofImage = nullptr;
static Color24 *dofBuffer = nullptr;

#define MAX_DOF_DRAW    32

//-------------------------------------------------------------------------------

void GlutDisplay();
void GlutReshape(int w, int h);
void GlutIdle();
void GlutKeyboard(unsigned char key, int x, int y);
void GlutMouse(int button, int state, int x, int y);
void GlutMotion(int x, int y);


// Custom reference I set
//--------------------------------------------------------------------------------
void ShowViewport();