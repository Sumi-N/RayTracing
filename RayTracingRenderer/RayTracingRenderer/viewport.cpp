#include "viewport.h"
#include <math.h>

extern MaterialList materials;

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

	float zero[] = { 0,0,0,0 };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, zero);

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

	const Material *mtl = node->GetMaterial();
	if (mtl) mtl->SetViewportMaterial();

	Matrix3f tm = node->GetTransform();
	Vec3f p = node->GetPosition();
	float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, p.x,p.y,p.z,1 };
	glMultMatrixf(m);

	Object *obj = node->GetNodeObj();
	if (obj) obj->ViewportDisplay(mtl);

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

	if (lights.size() > 0) {
		for (unsigned int i = 0; i < lights.size(); i++) {
			lights[i]->SetViewportLight(i);
		}
	}
	else {
		float white[] = { 1,1,1,1 };
		float black[] = { 0,0,0,0 };
		Vec4f p(camera.pos, 1);
		glEnable(GL_LIGHT0);
		glLightfv(GL_LIGHT0, GL_AMBIENT, black);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
		glLightfv(GL_LIGHT0, GL_SPECULAR, white);
		glLightfv(GL_LIGHT0, GL_POSITION, &p.x);
	}

	DrawNode(&rootNode);

	glPopMatrix();

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
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
		DrawImage(renderImage.GetPixels(), GL_UNSIGNED_BYTE, GL_RGB);
		DrawRenderProgressBar();
		break;
	case VIEWMODE_Z:
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
		if (!renderImage.GetZBufferImage()) renderImage.ComputeZBufferImage();
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

bool Sphere::IntersectRay(Ray const & ray, HitInfo & hInfo, int hitSide) const
{
	return false;
}

//-------------------------------------------------------------------------------
// Viewport Methods for various classes
//-------------------------------------------------------------------------------

Color FindReflection(Node * traversingnode, Node * node, Ray originalray, Ray ray, HitInfo & hit, int bounce)
{
	int numberofchild = traversingnode->GetNumChild();
	HitInfo hitinfo = HitInfo();
	for (int i = 0; i < numberofchild; i++)
	{
		node = traversingnode->GetChild(i);
		Ray changedray = node->ToNodeCoords(ray);

		if (node->GetNodeObj() != nullptr)
		{
			// This zbuffer is a fake, I don't use this buffer to do something
			float fakezbuffer;
			UpdateHitInfo(changedray, fakezbuffer, hitinfo, node);
			hit = hitinfo;
		}

		if (node != nullptr)
		{
			Node * childnode = new Node();
			FindReflection(node, childnode, originalray, changedray, hit, bounce);
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
	}
	else
	{
		return Color(0, 0, 0);
	}
}

Color Reflection(Ray const & ray, const HitInfo & hInfo, int bounce) 
{
	// bounce 1 is reflect 1 time
	if (bounce <= 0)
	{
		return Color(0,0,0);
	}
	else
	{
		//P is a surface point, 0.00003f is a bias
		Vec3f P = hInfo.N;
		P.Normalize(); P += 0.00003f; P += hInfo.p;

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

		Color returnColor = FindReflection(startnode, node, S, S, hitinfo, bounce - 1);
		delete node;
		return returnColor;
	}
}

Color FindRefraction(Node * traversingnode, Node * node, Ray originalray, Ray ray, HitInfo & hit, int bounce)
{
	int numberofchild = traversingnode->GetNumChild();
	HitInfo hitinfo = HitInfo();
	for (int i = 0; i < numberofchild; i++)
	{
		node = traversingnode->GetChild(i);
		Ray changedray = node->ToNodeCoords(ray);

		if (node->GetNodeObj() != nullptr)
		{
			// This zbuffer is a fake, I don't use this buffer to do something
			float fakezbuffer;
			UpdateHitInfo(changedray, fakezbuffer, hitinfo, node);
			hit = hitinfo;
		}

		if (node != nullptr)
		{
			Node * childnode = new Node();
			FindReflection(node, childnode, originalray, changedray, hit, bounce);
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
				Color returnColor = materials.Find(hitinfo.node->GetMaterial()->GetName())->Shade(originalray, hitinfo, lights, bounce);
				return returnColor;
			}
		}
	}
	else
	{
		return Color(0, 0, 0);
	}
}

Color Refraction(Ray const & ray, const HitInfo & hInfo, int bounce, float refractionIndex)
{
	// bounce 1 is refract 1 time
	if (bounce <= 0)
	{
		return Color(0, 0, 0);
	}
	else
	{
		// P is a surface point, 0.00003f is a bias
		Vec3f P = hInfo.N;
		P.Normalize(); P -= 0.00003f;  hInfo.p;

		Vec3f V = -1 * ray.dir;
		V.Normalize();
		Vec3f N = hInfo.N;
		N.Normalize();

		// Calculate angles 
		float cos1 = V.Dot(N); float sin1 = sqrt(1 - (cos1 * cos1));
		float sin2;
		if (hInfo.front)
		{
			sin2 = (1 / refractionIndex) * sin1;
		}
		else
		{
			sin2 = refractionIndex * sin1;
		}
		float cos2 = sqrt(1 - (sin2 * sin2));

		// Total internal reflection
		if (sin2 > 1)
		{
			return Color(0, 0, 0);	
		}

		// Horizontal dirction Vector
		Vec3f T_h = -cos2 * N;
		// Vertical Direction Vector
		Vec3f T_v = -sin2 * (V - (V.Dot(N))*N);
		Vec3f T = T_h + T_v;

		// S is a starting point from the surface point
		Ray S;
		S.dir = T;
		S.p = P;

		// Start traversing node 
		Node * node = new Node();
		Node * startnode = &rootNode;
		HitInfo hitinfo;

		Color returnColor = FindReflection(startnode, node, S, S, hitinfo, bounce - 1);
		delete node;
	}
}

void Sphere::ViewportDisplay(const Material *mtl) const
{
	static GLUquadric *q = nullptr;
	if (q == nullptr) {
		q = gluNewQuadric();
	}
	gluSphere(q, 1, 50, 50);
}
Color MtlBlinn::Shade(Ray const & ray, const HitInfo & hInfo, const LightList & lights, int bounce) const
{
	Vec3f N = hInfo.N;
	//N.Normalize();
	Color color = Color();
	for (auto light = lights.begin(); light != lights.end(); ++light) 
	{
		if ((*light)->IsAmbient()) 
		{
			color += (*light)->Illuminate(hInfo.p, N) * this->diffuse;
		}
		else if(strcmp((*light)->GetName(), "directionalLight") == 0)
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

			float oneofcos = 1/hInfo.N.Dot(L);
			Color specularpart = oneofcos * pow(H.Dot(hInfo.N), this->glossiness) * this->specular;
			Color diffusepart = this->diffuse;
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

			float oneofcos = 1/hInfo.N.Dot(L);
			Color specularpart = oneofcos * pow(H.Dot(hInfo.N), this->glossiness) * this->specular;
			Color diffusepart = this->diffuse;
			color += (diffusepart + specularpart) * IR;
		}
	}

	// Caluculate only relfection part for reflection
	if (this->reflection != Color(0, 0, 0))
	{
		color += this->reflection * Reflection(ray, hInfo, bounce);
	}

	//// Caluculate refraction part
	//if (this->refraction != Color(0, 0, 0))
	//{
	//	if (!hInfo.front)
	//	{

	//	}
	//	color += this->refraction * Refraction(ray, hInfo, bounce, ior);
	//}

	return color;
}
void MtlBlinn::SetViewportMaterial(int subMtlID) const
{
	ColorA c;
	c = ColorA(diffuse);
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, &c.r);
	c = ColorA(specular);
	glMaterialfv(GL_FRONT, GL_SPECULAR, &c.r);
	glMaterialf(GL_FRONT, GL_SHININESS, glossiness*1.5f);
}
void GenLight::SetViewportParam(int lightID, ColorA ambient, ColorA intensity, Vec4f pos) const
{
	glEnable(GL_LIGHT0 + lightID);
	glLightfv(GL_LIGHT0 + lightID, GL_AMBIENT, &ambient.r);
	glLightfv(GL_LIGHT0 + lightID, GL_DIFFUSE, &intensity.r);
	glLightfv(GL_LIGHT0 + lightID, GL_SPECULAR, &intensity.r);
	glLightfv(GL_LIGHT0 + lightID, GL_POSITION, &pos.x);
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

			float a = changedray.dir.Dot(changedray.dir);
			float b = 2 * changedray.dir.Dot(changedray.p);
			float c = changedray.p.Dot(changedray.p) - 1;

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
						answer = large;

						float delta = 0.00003f;
						if (answer > delta && answer <= t_max)
						{
							return true;
						}
					}
				}
				else
				{
					answer = small;

					float delta = 0.00003f;
					if (answer >= delta && answer <= t_max)
					{
						return true;
					}
				}
			}
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
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------