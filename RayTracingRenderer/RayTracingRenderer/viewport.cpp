#include "viewport.h"
#include "objects.h"
#include <math.h>

extern MaterialList materials;
extern TexturedColor environment;

#define SHADOWBIAS 0.0005f

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

	Color bg = background.GetColor();
	glClearColor(bg.r, bg.g, bg.b, 0);

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

	const TextureMap *bgMap = background.GetTexture();
	if (bgMap) {
		glDepthMask(GL_FALSE);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		Color c = background.GetColor();
		glColor3f(c.r, c.g, c.b);
		if (bgMap->SetViewportTexture()) {
			glEnable(GL_TEXTURE_2D);
			glMatrixMode(GL_TEXTURE);
			Matrix3f tm = bgMap->GetInverseTransform();
			Vec3f p = tm * bgMap->GetPosition();
			float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, -p.x,-p.y,-p.z,1 };
			glLoadMatrixf(m);
			glMatrixMode(GL_MODELVIEW);
		}
		else {
			glDisable(GL_TEXTURE_2D);
		}
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
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glDepthMask(GL_TRUE);

		glDisable(GL_TEXTURE_2D);
	}

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
	glDisable(GL_TEXTURE_2D);
}

//-------------------------------------------------------------------------------

void DrawImage(void *data, GLenum type, GLenum format)
{
	glBindTexture(GL_TEXTURE_2D, viewTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, renderImage.GetWidth(), renderImage.GetHeight(), 0, format, type, data);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	glEnable(GL_TEXTURE_2D);

	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();

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
	case 27:	// ESC
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
void Sphere::ViewportDisplay(const Material *mtl) const
{
	static GLUquadric *q = nullptr;
	if (q == nullptr) {
		q = gluNewQuadric();
		gluQuadricTexture(q, true);
	}
	gluSphere(q, 1, 50, 50);
}
void Plane::ViewportDisplay(const Material *mtl) const
{
	const int resolution = 32;
	float xyInc = 2.0f / resolution;
	float uvInc = 1.0f / resolution;
	glPushMatrix();
	glNormal3f(0, 0, 1);
	glBegin(GL_QUADS);
	float y1 = -1, y2 = xyInc - 1, v1 = 0, v2 = uvInc;
	for (int y = 0; y < resolution; y++) {
		float x1 = -1, x2 = xyInc - 1, u1 = 0, u2 = uvInc;
		for (int x = 0; x < resolution; x++) {
			glTexCoord2f(u1, v1);
			glVertex3f(x1, y1, 0);
			glTexCoord2f(u2, v1);
			glVertex3f(x2, y1, 0);
			glTexCoord2f(u2, v2);
			glVertex3f(x2, y2, 0);
			glTexCoord2f(u1, v2);
			glVertex3f(x1, y2, 0);
			x1 = x2; x2 += xyInc; u1 = u2; u2 += uvInc;
		}
		y1 = y2; y2 += xyInc; v1 = v2; v2 += uvInc;
	}
	glEnd();
	glPopMatrix();
}
void TriObj::ViewportDisplay(const Material *mtl) const
{
	unsigned int nextMtlID = 0;
	unsigned int nextMtlSwith = NF();
	if (mtl && NM() > 0) {
		mtl->SetViewportMaterial(0);
		nextMtlSwith = GetMaterialFaceCount(0);
		nextMtlID = 1;
	}
	glBegin(GL_TRIANGLES);
	for (unsigned int i = 0; i < NF(); i++) {
		while (i >= nextMtlSwith) {
			if (nextMtlID >= NM()) nextMtlSwith = NF();
			else {
				glEnd();
				nextMtlSwith += GetMaterialFaceCount(nextMtlID);
				mtl->SetViewportMaterial(nextMtlID);
				nextMtlID++;
				glBegin(GL_TRIANGLES);
			}
		}
		for (int j = 0; j < 3; j++) {
			if (HasTextureVertices()) glTexCoord3fv(&VT(FT(i).v[j]).x);
			if (HasNormals()) glNormal3fv(&VN(FN(i).v[j]).x);
			glVertex3fv(&V(F(i).v[j]).x);
		}
	}
	glEnd();
}
void MtlBlinn::SetViewportMaterial(int subMtlID) const
{
	ColorA c;
	c = ColorA(diffuse.GetColor());
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, &c.r);
	c = ColorA(specular.GetColor());
	glMaterialfv(GL_FRONT, GL_SPECULAR, &c.r);
	glMaterialf(GL_FRONT, GL_SHININESS, glossiness*1.5f);
	const TextureMap *dm = diffuse.GetTexture();
	if (dm && dm->SetViewportTexture()) {
		glEnable(GL_TEXTURE_2D);
		glMatrixMode(GL_TEXTURE);
		Matrix3f tm = dm->GetInverseTransform();
		Vec3f p = tm * dm->GetPosition();
		float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, -p.x,-p.y,-p.z,1 };
		glLoadMatrixf(m);
		glMatrixMode(GL_MODELVIEW);
	}
	else {
		glDisable(GL_TEXTURE_2D);
	}
}
void GenLight::SetViewportParam(int lightID, ColorA ambient, ColorA intensity, Vec4f pos) const
{
	glEnable(GL_LIGHT0 + lightID);
	glLightfv(GL_LIGHT0 + lightID, GL_AMBIENT, &ambient.r);
	glLightfv(GL_LIGHT0 + lightID, GL_DIFFUSE, &intensity.r);
	glLightfv(GL_LIGHT0 + lightID, GL_SPECULAR, &intensity.r);
	glLightfv(GL_LIGHT0 + lightID, GL_POSITION, &pos.x);
}

bool TextureFile::SetViewportTexture() const
{
	if (viewportTextureID == 0) {
		glGenTextures(1, &viewportTextureID);
		glBindTexture(GL_TEXTURE_2D, viewportTextureID);
		gluBuild2DMipmaps(GL_TEXTURE_2D, 3, width, height, GL_RGB, GL_UNSIGNED_BYTE, &data[0].r);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	}
	glBindTexture(GL_TEXTURE_2D, viewportTextureID);
	return true;
}

bool TextureChecker::SetViewportTexture() const
{
	if (viewportTextureID == 0) {
		const int texSize = 256;
		glGenTextures(1, &viewportTextureID);
		glBindTexture(GL_TEXTURE_2D, viewportTextureID);
		Color24 c[2] = { Color24(color1), Color24(color2) };
		Color24 *tex = new Color24[texSize*texSize];
		for (int i = 0; i < texSize*texSize; i++) {
			int ix = (i%texSize) < 128 ? 0 : 1;
			if (i / 256 >= 128) ix = 1 - ix;
			tex[i] = c[ix];
		}
		gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texSize, texSize, GL_RGB, GL_UNSIGNED_BYTE, &tex[0].r);
		delete[] tex;
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	}
	glBindTexture(GL_TEXTURE_2D, viewportTextureID);
	return true;
}
//-------------------------------------------------------------------------------
// Custom class I made
//-------------------------------------------------------------------------------

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
		return Color(0,0,0);
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
			P = - 1 * SHADOWBIAS * P;  P += hInfo.p;
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
		Color returnColor = (1-R) * refraction * FindReflectionAndRefraction(startnode, node, S, S, hitinfo, bounce -1);
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
		else if(strcmp((*light)->GetName(), "directLight") == 0)
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

			float oneofcos = 1/hInfo.N.Dot(L);
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
					if (answer > SHADOWBIAS & answer <= hInfo.z)
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
						float v = (1 / 3.14) * asinf(hInfo.p.z) + 0.5f;
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
				if (answer > SHADOWBIAS & answer <= hInfo.z)
				{
					return true;
				}
			}else
			{
				if (CheckZbuffer(hInfo.z, small))
				{
					hInfo.z = small;
					hInfo.front = true;
					hInfo.p = ray.p + small * ray.dir;
					hInfo.N = hInfo.p;
					float u = (1 / (2 * 3.14f)) * atan2f(hInfo.p.y, hInfo.p.x) + .5f;
					float v = (1 / 3.14) * asinf(hInfo.p.z) + 0.5f;
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
				if (t > SHADOWBIAS & t <= hInfo.z)
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


				//float t2 = -1 * ((ray.p + hInfo.duvw[0]).Dot(Vec3f(0, 0, 1))) / (ray.p + hInfo.duvw[0]).Dot(Vec3f(0, 0, 1));
				//Vec3f point2 = ray.p + t2 * (ray.dir + hInfo.duvw[0]);
				//hInfo.duvw[0] = point2 - point;

				//float t3 = -1 * ((ray.p + hInfo.duvw[1]).Dot(Vec3f(0, 0, 1))) / (ray.p + hInfo.duvw[1]).Dot(Vec3f(0, 0, 1));
				//Vec3f point3 = ray.p + t3 * (ray.dir + hInfo.duvw[1]);
				//hInfo.duvw[1] = point3 - point;

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

	float x, y, z;

	xmin = (vertices[0] - ray.p.x) / ray.dir.x;
	y = xmin * ray.dir.y + ray.p.y;
	z = xmin * ray.dir.z + ray.p.z;
	if (y >= vertices[1] && y <= vertices[4] && z >= vertices[2] && z <= vertices[5])
	{
		if (detected)
		{
			maxofmin = xmin;
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
		}
		else
		{
			minofmax = ymax;
			detected = true;
		}
	}

	if (!detected)
	{
		return false;
	}

	zmax = (vertices[5] - ray.p.z) / ray.dir.z;
	x = zmax * ray.dir.x + ray.p.x;
	y = zmax * ray.dir.y + ray.p.y;
	if (x >= vertices[0] && x <= vertices[3] && y >= vertices[1] && y <= vertices[4])
	{
		if (detected)
		{
			maxofmin = zmax;
		}
		else
		{
			minofmax = zmax;
			detected = true;
		}
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
		for (int i = 0; i < nodecount; i++)
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
		
		float child1answer, child2answer;

		if (!CheckBoxCollision(vertices1, ray, child1answer) & !(CheckBoxCollision(vertices2, ray, child2answer)))
		{
			return false;
		}

		if (child1answer < child2answer)
		{
			if (child1answer <= hInfo.z)
			{
				bool hit2 = false;
				if (TraceBVHNode(ray, hInfo, hitSide, child1id))
				{
					hit2 = true;

					if (child2answer <= hInfo.z)
					{
						if (TraceBVHNode(ray, hInfo, hitSide, child1id))
						{
							hit2 = true;
						}
					}
				}
				else if(child2answer <= hInfo.z)
				{
					if (TraceBVHNode(ray, hInfo, hitSide, child2id))
					{
						hit2 = true;
					}
				}
				return hit2;
			}
		}
		else if (child1answer > child2answer)
		{
			if (child2answer <= hInfo.z)
			{
				bool hit2 = false;
				if (TraceBVHNode(ray, hInfo, hitSide, child2id))
				{
					hit2 = true;

					if (child2answer <= hInfo.z)
					{
						if (TraceBVHNode(ray, hInfo, hitSide, child1id))
						{
							hit2 = true;
						}
					}
				}
				else if(child2answer <= hInfo.z)
				{
					if (TraceBVHNode(ray, hInfo, hitSide, child1id))
					{
						hit2 = true;
					}
				}
				return hit2;
			}
		}
		else if (child1answer == child2answer)
		{
			bool hit2 = false;
			if (TraceBVHNode(ray, hInfo, hitSide, child1id))
			{
				hit2 = true;
			}
			if (TraceBVHNode(ray, hInfo, hitSide, child2id))
			{
				hit2 = true;
			}
			return hit2;
		}
	}
	return false;
}


