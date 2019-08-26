//-------------------------------------------------------------------------------
///
/// \file       viewport.cpp 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    1.0
/// \date       August 21, 2019
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#include <viewport.h>

//-------------------------------------------------------------------------------

void ShowViewport()
{
	int argc = 1;
	char argstr[] = "raytrace";
	char *argv = argstr;
	glutInit(&argc, &argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	if (glutGet(GLUT_SCREEN_WIDTH) > 0 && glutGet(GLUT_SCREEN_HEIGHT) > 0) {
		glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - camera.imgWidth) / 2, (glutGet(GLUT_SCREEN_HEIGHT) - camera.imgHeight) / 2);
	}
	else glutInitWindowPosition(50, 50);
	glutInitWindowSize(camera.imgWidth, camera.imgHeight);

	glutCreateWindow("Ray Tracer - CS 6620");
	glutDisplayFunc(GlutDisplay);
	glutReshapeFunc(GlutReshape);
	glutIdleFunc(GlutIdle);
	glutKeyboardFunc(GlutKeyboard);
	glutMouseFunc(GlutMouse);
	glutMotionFunc(GlutMotion);

	glClearColor(0, 0, 0, 0);

	glPointSize(3.0);
	glEnable(GL_CULL_FACE);

#define LIGHTAMBIENT 0.1f
	glEnable(GL_LIGHT0);
	float lightamb[4] = { LIGHTAMBIENT, LIGHTAMBIENT, LIGHTAMBIENT, 1.0f };
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightamb);

#define LIGHTDIF0 1.0f
	float lightdif0[4] = { LIGHTDIF0, LIGHTDIF0, LIGHTDIF0, 1.0f };
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightdif0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightdif0);

	glEnable(GL_NORMALIZE);

	glLineWidth(2);

	glGenTextures(1, &viewTexture);
	glBindTexture(GL_TEXTURE_2D, viewTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	glutMainLoop();
}

//-------------------------------------------------------------------------------

void GlutReshape(int w, int h)
{
	if (w != camera.imgWidth || h != camera.imgHeight) {
		glutReshapeWindow(camera.imgWidth, camera.imgHeight);
	}
	else {
		glViewport(0, 0, w, h);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		float r = (float)w / float(h);
		gluPerspective(camera.fov, r, 0.02, 1000.0);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
}

//-------------------------------------------------------------------------------

void DrawNode(Node *node)
{
	glPushMatrix();

	Matrix3f tm = node->GetTransform();
	Vec3f p = node->GetPosition();
	float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, p.x,p.y,p.z,1 };
	glMultMatrixf(m);

	Object *obj = node->GetNodeObj();
	if (obj) obj->ViewportDisplay();

	for (int i = 0; i < node->GetNumChild(); i++) {
		DrawNode(node->GetChild(i));
	}

	glPopMatrix();
}

//-------------------------------------------------------------------------------

void DrawScene()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);

	glPushMatrix();
	Vec3f p = camera.pos;
	Vec3f t = camera.pos + camera.dir;
	Vec3f u = camera.up;
	gluLookAt(p.x, p.y, p.z, t.x, t.y, t.z, u.x, u.y, u.z);

	glRotatef(viewAngle1, 1, 0, 0);
	glRotatef(viewAngle2, 0, 0, 1);

	DrawNode(&rootNode);

	glPopMatrix();

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
}

//-------------------------------------------------------------------------------

void DrawProgressBar(float done)
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glBegin(GL_LINES);
	glColor3f(1, 1, 1);
	glVertex2f(-1, -1);
	glVertex2f(done * 2 - 1, -1);
	glColor3f(0, 0, 0);
	glVertex2f(done * 2 - 1, -1);
	glVertex2f(1, -1);
	glEnd();

	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

//-------------------------------------------------------------------------------

void DrawImage(void *data, GLenum type, GLenum format)
{
	glBindTexture(GL_TEXTURE_2D, viewTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, renderImage.GetWidth(), renderImage.GetHeight(), 0, format, type, data);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	glEnable(GL_TEXTURE_2D);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glColor3f(1, 1, 1);
	glBegin(GL_QUADS);
	glTexCoord2f(0, 1);
	glVertex2f(-1, -1);
	glTexCoord2f(1, 1);
	glVertex2f(1, -1);
	glTexCoord2f(1, 0);
	glVertex2f(1, 1);
	glTexCoord2f(0, 0);
	glVertex2f(-1, 1);
	glEnd();

	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

	glDisable(GL_TEXTURE_2D);
}

//-------------------------------------------------------------------------------

void DrawRenderProgressBar()
{
	int rp = renderImage.GetNumRenderedPixels();
	int np = renderImage.GetWidth() * renderImage.GetHeight();
	if (rp >= np) return;
	float done = (float)rp / (float)np;
	DrawProgressBar(done);
}

//-------------------------------------------------------------------------------

void GlutDisplay()
{
	switch (viewMode) {
	case VIEWMODE_OPENGL:
		DrawScene();
		break;
	case VIEWMODE_IMAGE:
		//glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
		//glDrawPixels( renderImage.GetWidth(), renderImage.GetHeight(), GL_RGB, GL_UNSIGNED_BYTE, renderImage.GetPixels() );
		DrawImage(renderImage.GetPixels(), GL_UNSIGNED_BYTE, GL_RGB);
		DrawRenderProgressBar();
		break;
	case VIEWMODE_Z:
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
		if (!renderImage.GetZBufferImage()) renderImage.ComputeZBufferImage();
		//glDrawPixels( renderImage.GetWidth(), renderImage.GetHeight(), GL_LUMINANCE, GL_UNSIGNED_BYTE, renderImage.GetZBufferImage() );
		DrawImage(renderImage.GetZBufferImage(), GL_UNSIGNED_BYTE, GL_LUMINANCE);
		break;
	}

	glutSwapBuffers();
}

//-------------------------------------------------------------------------------

void GlutIdle()
{
	static int lastRenderedPixels = 0;
	if (mode == MODE_RENDERING) {
		int nrp = renderImage.GetNumRenderedPixels();
		if (lastRenderedPixels != nrp) {
			lastRenderedPixels = nrp;
			if (renderImage.IsRenderDone()) {
				mode = MODE_RENDER_DONE;
				int endTime = (int)time(nullptr);
				int t = endTime - startTime;
				int h = t / 3600;
				int m = (t % 3600) / 60;
				int s = t % 60;
				printf("\nRender time is %d:%02d:%02d.\n", h, m, s);
			}
			glutPostRedisplay();
		}
	}
}

//-------------------------------------------------------------------------------

void GlutKeyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 27:    // ESC
		exit(0);
		break;
	case ' ':
		switch (mode) {
		case MODE_READY:
			mode = MODE_RENDERING;
			viewMode = VIEWMODE_IMAGE;
			DrawScene();
			glReadPixels(0, 0, renderImage.GetWidth(), renderImage.GetHeight(), GL_RGB, GL_UNSIGNED_BYTE, renderImage.GetPixels());
			{
				Color24 *c = renderImage.GetPixels();
				for (int y0 = 0, y1 = renderImage.GetHeight() - 1; y0 < y1; y0++, y1--) {
					int i0 = y0 * renderImage.GetWidth();
					int i1 = y1 * renderImage.GetWidth();
					for (int x = 0; x < renderImage.GetWidth(); x++, i0++, i1++) {
						Color24 t = c[i0]; c[i0] = c[i1]; c[i1] = t;
					}
				}
			}
			startTime = (int)time(nullptr);
			BeginRender();
			break;
		case MODE_RENDERING:
			mode = MODE_READY;
			StopRender();
			glutPostRedisplay();
			break;
		case MODE_RENDER_DONE:
			mode = MODE_READY;
			viewMode = VIEWMODE_OPENGL;
			glutPostRedisplay();
			break;
		}
		break;
	case '1':
		viewAngle1 = viewAngle2 = 0;
		viewMode = VIEWMODE_OPENGL;
		glutPostRedisplay();
		break;
	case '2':
		viewMode = VIEWMODE_IMAGE;
		glutPostRedisplay();
		break;
	case '3':
		viewMode = VIEWMODE_Z;
		glutPostRedisplay();
		break;
	}
}

//-------------------------------------------------------------------------------

void PrintPixelData(int x, int y)
{
	if (x < renderImage.GetWidth() && y < renderImage.GetHeight()) {
		Color24 *colors = renderImage.GetPixels();
		float *zbuffer = renderImage.GetZBuffer();
		int i = (renderImage.GetHeight() - y - 1) *renderImage.GetWidth() + x;
		printf("Pixel [ %d, %d ] Color24: %d, %d, %d   Z: %f\n", x, y, colors[i].r, colors[i].g, colors[i].b, zbuffer[i]);
	}
	else {
		printf("-- Invalid pixel (%d,%d) --\n", x, y);
	}
}

//-------------------------------------------------------------------------------

void GlutMouse(int button, int state, int x, int y)
{
	if (state == GLUT_UP) {
		mouseMode = MOUSEMODE_NONE;
	}
	else {
		switch (button) {
		case GLUT_LEFT_BUTTON:
			mouseMode = MOUSEMODE_DEBUG;
			PrintPixelData(x, y);
			break;
		case GLUT_RIGHT_BUTTON:
			mouseMode = MOUSEMODE_ROTATE;
			mouseX = x;
			mouseY = y;
			break;
		}
	}
}

//-------------------------------------------------------------------------------

void GlutMotion(int x, int y)
{
	switch (mouseMode) {
	case MOUSEMODE_DEBUG:
		PrintPixelData(x, y);
		break;
	case GLUT_RIGHT_BUTTON:
		viewAngle1 -= 0.2f * (mouseY - y);
		viewAngle2 -= 0.2f * (mouseX - x);
		mouseX = x;
		mouseY = y;
		glutPostRedisplay();
		break;
	}
}

//-------------------------------------------------------------------------------
// Viewport Methods for various classes
//-------------------------------------------------------------------------------
void Sphere::ViewportDisplay() const
{
	static GLUquadric *q = nullptr;
	if (q == nullptr) {
		q = gluNewQuadric();
	}
	gluSphere(q, 1, 50, 50);
}

bool Sphere::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide) const {
	return true;
}
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

void TraverseNode(Node * traversingnode, Node ** nodes, int & cn ) {
	int numberofchild = traversingnode->GetNumChild();
	for (int i = cn; i < numberofchild + cn; i++) {
		nodes[i] = traversingnode->GetChild(i);
		cn++;
		if (nodes[i] != nullptr) {
			TraverseNode(nodes[i],nodes,cn);
		}
	}
	//cn += numberofchild;
}

void BeginRender() {

	Color24* pixels = renderImage.GetPixels();

	Vec3f * cameraray = new Vec3f[renderImage.GetHeight() * renderImage.GetWidth()];

	float l = 1.0f;
	float h = 2 * l * tanf((camera.fov/2) * 3.14 /180);
	float w = camera.imgWidth * ( h / camera.imgHeight);

	int H = camera.imgHeight;
	int W = camera.imgWidth;

	Vec3f x = camera.up.Cross(camera.dir);
	Vec3f y = camera.up;

	Vec3f f = camera.pos + l * camera.dir + (h / 2) * y - (w / 2) * x;

	for (int i = 0; i < renderImage.GetHeight(); i++) {
		for (int j = 0; j < renderImage.GetWidth(); j++) {
			Vec3f mockcicle = Vec3f(0, 50, -25);
			cameraray[i * renderImage.GetWidth() + j] = f + (j + 0.5f) * (w / W)*x - (i + 0.5f) * (h / H)*y - camera.pos;

			Vec3f d = cameraray[i * renderImage.GetWidth() + j];
			Vec3f e = camera.pos - mockcicle;

			float a = d.Dot(d);
			float b = 2 * d.Dot(e);
			float c = e.Dot(e) - 625;

			if (b*b - 4*a*c >= 0) {
				pixels[i * renderImage.GetWidth() + j].r = 0;
				pixels[i * renderImage.GetWidth() + j].b = 0;
				pixels[i * renderImage.GetWidth() + j].g = 255;
			}
		}
	}


	for (int i = 0; i < renderImage.GetHeight(); i++) {
		for (int j = 0; j < renderImage.GetWidth(); j++) {
			Vec3f mockcicle = Vec3f(0, 50, 5.1);
			cameraray[i * renderImage.GetWidth() + j] = f + (j + 0.5f) * (w / W)*x - (i + 0.5f) * (h / H)*y - camera.pos;

			Vec3f d = cameraray[i * renderImage.GetWidth() + j];
			Vec3f e = camera.pos - mockcicle;

			float a = d.Dot(d);
			float b = 2 * d.Dot(e);
			float c = e.Dot(e) - 25;

			if (b*b - 4 * a*c >= 0) {
				pixels[i * renderImage.GetWidth() + j].r = 0;
				pixels[i * renderImage.GetWidth() + j].b = 255;
				pixels[i * renderImage.GetWidth() + j].g = 0;
			}
		}
	}

	for (int i = 0; i < renderImage.GetHeight(); i++) {
		for (int j = 0; j < renderImage.GetWidth(); j++) {
			Vec3f mockcicle = Vec3f(0, 50, 11.1);
			cameraray[i * renderImage.GetWidth() + j] = f + (j + 0.5f) * (w / W)*x - (i + 0.5f) * (h / H)*y - camera.pos;

			Vec3f d = cameraray[i * renderImage.GetWidth() + j];
			Vec3f e = camera.pos - mockcicle;

			float a = d.Dot(d);
			float b = 2 * d.Dot(e);
			float c = e.Dot(e) - 1; //6.8 is the magic number

			if (b*b - 4 * a*c >= 0) {
				pixels[i * renderImage.GetWidth() + j].r = 255;
				pixels[i * renderImage.GetWidth() + j].b = 0;
				pixels[i * renderImage.GetWidth() + j].g = 0;
			}
		}
	}

	renderImage.SaveImage("hello.png");

	return;
}

void StopRender() {
	return;
}