// glfw_30.cpp : Defines the entry point for the console application.
//  http://www.glfw.org/docs/latest/quick.html
#include "stdafx.h"
#define _CRT_SECURE_NO_DEPRECATE
#include <time.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <GL/GL.h>
#include <GL/glaux.h>
#pragma comment(lib, "glaux.lib")
#include <GLFW/glfw3.h>
#include <vector>
#include <Windows.h>
#include <iostream>

#define SCREEN_WIDTH 1000
#define SCREEN_HEIGHT 1000
#define PC_MODE 16
#define CONSOLE_MODE 0

#define MAXRAD 4
#define MINRAD 1
#define MAXHEIGHT 7
#define MINHEIGHT 2
#define MAXFREQ_X 40
#define MINFREQ_X 3
#define ANGLE 2
#define TWINKIE_SPEED 20
#define MAIN_LIST 1


//Optimization goes here:
//enable PERFORM_TESTS and enable optimization
//results (FPS after optimization) after 10 tests will be stored in x_file.txt

//#define PERFORM_TESTS

//#define X_NOSMOOTH			//disable built-in antialiasing
//#define X_CULLFACE			//hide surfaces that cant be seen
//#define X_SCREENSIZE			//reduce screen size
//#define X_NORAMLIZE			//disable built-in normalization of normals
//#define X_TEXSIZE			//reduce texture size
//#define X_LIGHTDECREASE		//disable 2-sided light model
//#define X_DISPLIST			//use display list for drawing


typedef struct {
	GLfloat x, y, z;
} point;

std::vector<std::vector<point>> cone_top;
std::vector<std::vector<point>> cone_bot;
std::vector<std::vector<point>> cone_side;

std::vector<std::vector<point>> cyl_top;
std::vector<std::vector<point>> cyl_bot;
std::vector<std::vector<point>> cyl_side;

bool shift, twinkie_flag, tex_flag, light_flag;
int type, stop_i, stop_j, freq_x, freq_y, list;
float t, twinkie_steps = 100, dt = 5 / twinkie_steps;
GLfloat light_pos_x, light_pos_y, light_pos_z;
static float alpha, beta, gamma, h, r;


static void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}

point count_normal(point p1, point p2, point p3)
{
	point a, b, n;
	GLfloat l;

	a.x = p2.x - p1.x;
	a.y = p2.y - p1.y;
	a.z = p2.z - p1.z;

	b.x = p3.x - p1.x;
	b.y = p3.y - p1.y;
	b.z = p3.z - p1.z;

	n.x = (a.y * b.z) - (a.z * b.y);
	n.y = (a.z * b.x) - (a.x * b.z);
	n.z = (a.x * b.y) - (a.y * b.x);

	// Normalize (divide by root of dot product)
	l = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
	n.x /= l;
	n.y /= l;
	n.z /= l;

	return n;
}

GLuint texture[1];

void load_textures()
{
	AUX_RGBImageRec *texture1;

#ifdef X_TEXSIZE
	texture1 = auxDIBImageLoadA("tex_cage_small.bmp");
#endif

#ifndef X_TEXSIZE
	texture1 = auxDIBImageLoadA("tex_cage.bmp");
#endif

	glGenTextures(1, &texture[0]);
	glBindTexture(GL_TEXTURE_2D, texture[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, texture1->sizeX, texture1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, texture1->data);

	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
}

void light_enable()
{

#ifdef X_LIGHTDECRESE
	glShadeModel(GL_FLAT);
#endif

#ifndef X_LIGHTDECREASE
	glShadeModel(GL_SMOOTH);
#endif

#ifndef X_LIGHTDECREASE
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
#endif

	GLfloat material_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);


	GLfloat light0_diffuse[] = { 0.4, 0.7, 0.2 };
	GLfloat light0_position[] = { light_pos_x, light_pos_y, light_pos_z, 1.0 };
	//GLfloat spot_direction[] = { -1.0, -1.0, 0.0 };


	//glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 45.0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	//glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, spot_direction);
	//glLightfv(GL_LIGHT0, GL_SPOT_EXPONENT, );


	glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.5);
	glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05);
	glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.05);

#ifdef X_NORMALIZE
	glEnable(GL_NORMALIZE);
#endif
}

void count_cone()
{
	int i, j, top_rad = r, bot_rad = r * 2;
	float angle, dangle, height, dheight, drad0, drad1, drad2;



	height = h;
	dheight = height / (float)freq_y;
	dangle = 2 * M_PI / (float)freq_x;
	drad0 = top_rad / (float)freq_y;
	drad1 = (bot_rad - top_rad) / (float)freq_y;
	drad2 = bot_rad / (float)freq_y;



	cone_top.resize(freq_y + 1);
	for (i = 0; i < freq_y + 1; i++) cone_top[i].resize(freq_x + 1);

	cone_bot.resize(freq_y + 1);
	for (i = 0; i < freq_y + 1; i++) cone_bot[i].resize(freq_x + 1);

	cone_side.resize(freq_y + 1);
	for (i = 0; i < freq_y + 1; i++) cone_side[i].resize(freq_x + 1);



	for (i = 0; i < freq_y + 1; i += 1) {
		for (j = 0, angle = 0; angle <= 2 * M_PI; j += 1, angle += dangle) {
			cone_top[i][j].x = drad0 * i * cos(angle);
			cone_top[i][j].y = height;
			cone_top[i][j].z = drad0 * i * sin(angle);

			cone_side[i][j].x = (top_rad + drad1 * i) * cos(angle);
			cone_side[i][j].y = height - dheight * i;
			cone_side[i][j].z = (top_rad + drad1 * i) * sin(angle);

			cone_bot[i][j].x = drad2 * i * cos(angle);
			cone_bot[i][j].y = 0;
			cone_bot[i][j].z = drad2 * i * sin(angle);
		}
	}
	stop_i = i;
	stop_j = j;

	printf("cone_renewed: h=%.2f, top_r=%d, fr_x=%d, fr_y=%d s_i=%d s_j=%d\n", height, top_rad, freq_x, freq_y, stop_i, stop_j);
}

void count_cylinder()
{
	int i, j, rad = r;
	float angle, dangle, height, dheight, drad0;


	height = h;
	dheight = height / (float)freq_y;
	dangle = 2 * M_PI / (float)freq_x;
	drad0 = rad / (float)freq_y;



	cyl_top.resize(freq_y + 1);
	for (i = 0; i < freq_y + 1; i++) cyl_top[i].resize(freq_x + 1);

	cyl_bot.resize(freq_y + 1);
	for (i = 0; i < freq_y + 1; i++) cyl_bot[i].resize(freq_x + 1);

	cyl_side.resize(freq_y + 1);
	for (i = 0; i < freq_y + 1; i++) cyl_side[i].resize(freq_x + 1);



	for (i = 0; i < freq_y + 1; i += 1) {
		for (j = 0, angle = 0; angle <= 2 * M_PI; j += 1, angle += dangle) {
			cyl_top[i][j].x = drad0 * i * cos(angle);
			cyl_top[i][j].y = height;
			cyl_top[i][j].z = drad0 * i * sin(angle);

			cyl_side[i][j].x = rad * cos(angle);
			cyl_side[i][j].y = height - dheight * i;
			cyl_side[i][j].z = rad * sin(angle);

			cyl_bot[i][j].x = drad0 * i * cos(angle);
			cyl_bot[i][j].y = 0;
			cyl_bot[i][j].z = drad0 * i * sin(angle);
		}
	}
	stop_i = i;
	stop_j = j;

	printf("cyl_renewed: h=%.2f, r=%d, fr_x=%d, fr_y=%d s_i=%d s_j=%d\n", height, rad, freq_x, freq_y, stop_i, stop_j);
}

void draw_figure(
	std::vector<std::vector<point>>& top,
	std::vector<std::vector<point>>& side,
	std::vector<std::vector<point>>& bot
	)
{
	int i, j;
	point normal;
	GLfloat tex_i, tex_j, dtex_i, dtex_j;

	dtex_i = 1.0 / (float)(freq_y);
	dtex_j = 1.0 / (float)(freq_x);

	glRotatef(alpha, 1, 0, 0);
	glRotatef(beta, 0, 1, 0);
	glRotatef(gamma, 0, 0, 1);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glBindTexture(GL_TEXTURE_2D, texture[0]);

#ifdef X_DISPLIST
	list = glGenLists(1);
	if (list != 0) {
		glNewList(list, GL_COMPILE);
#endif
		for (i = 0, tex_i = 0; i < stop_i - 1; i += 1, tex_i += dtex_i) {
			for (j = 0, tex_j = 0; j < stop_j; j += 1, tex_j += dtex_j) {
				//top
				glBegin(GL_POLYGON);
				normal = count_normal(top[i][j], top[i + 1][(j + 1) % stop_j], top[i + 1][j]);
				glNormal3f(normal.x, normal.y, normal.z);

				glColor3f(1, 1, 1);
				glVertex3f(top[i][j].x, top[i][j].y, top[i][j].z);
				glVertex3f(top[i][(j + 1) % stop_j].x, top[i][(j + 1) % stop_j].y, top[i][(j + 1) % stop_j].z);
				glVertex3f(top[i + 1][(j + 1) % stop_j].x, top[i + 1][(j + 1) % stop_j].y, top[i + 1][(j + 1) % stop_j].z);
				glVertex3f(top[i + 1][j].x, top[i + 1][j].y, top[i + 1][j].z);
				glEnd();


				//sidelines
				glBegin(GL_POLYGON);
				normal = count_normal(side[i][j], side[i + 1][(j + 1) % stop_j], side[i + 1][j]);
				glNormal3f(normal.x, normal.y, normal.z);

				glColor3f(0.6, 0.6, 0.6);
				glTexCoord2f(tex_i, tex_j);
				glVertex3f(side[i][j].x, side[i][j].y, side[i][j].z);
				glTexCoord2f(tex_i, tex_j + dtex_j);
				glVertex3f(side[i][(j + 1) % stop_j].x, side[i][(j + 1) % stop_j].y, side[i][(j + 1) % stop_j].z);
				glTexCoord2f(tex_i + dtex_i, tex_j + dtex_j);
				glVertex3f(side[i + 1][(j + 1) % stop_j].x, side[i + 1][(j + 1) % stop_j].y, side[i + 1][(j + 1) % stop_j].z);
				glTexCoord2f(tex_i + dtex_i, tex_j);
				glVertex3f(side[i + 1][j].x, side[i + 1][j].y, side[i + 1][j].z);
				glEnd();


				//bot
				glBegin(GL_POLYGON);
				normal = count_normal(bot[i][j], bot[i + 1][(j + 1) % stop_j], bot[i + 1][j]);
				glNormal3f(normal.x, normal.y, normal.z);

				glColor3f(0.3, 0.3, 0.3);
				glTexCoord2f(0.0f, 0.0f);
				glVertex3f(bot[i][j].x, bot[i][j].y, bot[i][j].z);
				glTexCoord2f(0.0f, 1.0f);
				glVertex3f(bot[i][(j + 1) % stop_j].x, bot[i][(j + 1) % stop_j].y, bot[i][(j + 1) % stop_j].z);
				glTexCoord2f(1.0f, 1.0f);
				glVertex3f(bot[i + 1][(j + 1) % stop_j].x, bot[i + 1][(j + 1) % stop_j].y, bot[i + 1][(j + 1) % stop_j].z);
				glTexCoord2f(1.0f, 0.0f);
				glVertex3f(bot[i + 1][j].x, bot[i + 1][j].y, bot[i + 1][j].z);
				glEnd();
			}
		}

		glBegin(GL_LINES);
		glColor3f(1, 0, 0);
		glVertex3f(light_pos_x, light_pos_y, light_pos_z);
		glVertex3f(0, 0, 0);
		glEnd();

#ifdef X_DISPLIST
	}
	glEndList();
#endif
}

GLFWwindow* WINDOW;

void tweenking()
{
	int i, j, speed;
	//Bezier cube:
	//(1-t)^3*p0 + 3t(1-t)^2*p1 + 3*t^2*(1-t)*p2 + t^3*p3

	for (t = 0, speed = 0; speed <= 20; t += dt, speed++) {
		for (i = 0; i < stop_i; i += 1) {
			for (j = 0; j < stop_j; j += 1) {
				cone_top[i][j].x = (1 - t)*(1 - t)*(1 - t)*cone_top[i][j].x + 3 * t*(1 - t)*(1 - t)*cone_top[i][j].x + 3 * t*t*(1 - t)*cone_top[i][j].x + t*t*t*cyl_top[i][j].x;
				cone_top[i][j].y = (1 - t)*(1 - t)*(1 - t)*cone_top[i][j].y + 3 * t*(1 - t)*(1 - t)*cone_top[i][j].y + 3 * t*t*(1 - t)*cone_top[i][j].y + t*t*t*cyl_top[i][j].y;
				cone_top[i][j].z = (1 - t)*(1 - t)*(1 - t)*cone_top[i][j].z + 3 * t*(1 - t)*(1 - t)*cone_top[i][j].z + 3 * t*t*(1 - t)*cone_top[i][j].z + t*t*t*cyl_top[i][j].z;

				cone_side[i][j].x = (1 - t)*(1 - t)*(1 - t)*cone_side[i][j].x + 3 * t*(1 - t)*(1 - t)*cone_side[i][j].x + 3 * t*t*(1 - t)*cone_side[i][j].x + t*t*t*cyl_side[i][j].x;
				cone_side[i][j].y = (1 - t)*(1 - t)*(1 - t)*cone_side[i][j].y + 3 * t*(1 - t)*(1 - t)*cone_side[i][j].y + 3 * t*t*(1 - t)*cone_side[i][j].y + t*t*t*cyl_side[i][j].y;
				cone_side[i][j].z = (1 - t)*(1 - t)*(1 - t)*cone_side[i][j].z + 3 * t*(1 - t)*(1 - t)*cone_side[i][j].z + 3 * t*t*(1 - t)*cone_side[i][j].z + t*t*t*cyl_side[i][j].z;

				cone_bot[i][j].x = (1 - t)*(1 - t)*(1 - t)*cone_bot[i][j].x + 3 * t*(1 - t)*(1 - t)*cone_bot[i][j].x + 3 * t*t*(1 - t)*cone_bot[i][j].x + t*t*t*cyl_bot[i][j].x;
				cone_bot[i][j].y = (1 - t)*(1 - t)*(1 - t)*cone_bot[i][j].y + 3 * t*(1 - t)*(1 - t)*cone_bot[i][j].y + 3 * t*t*(1 - t)*cone_bot[i][j].y + t*t*t*cyl_bot[i][j].y;
				cone_bot[i][j].z = (1 - t)*(1 - t)*(1 - t)*cone_bot[i][j].z + 3 * t*(1 - t)*(1 - t)*cone_bot[i][j].z + 3 * t*t*(1 - t)*cone_bot[i][j].z + t*t*t*cyl_bot[i][j].z;
			}
		}
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		draw_figure(cone_top, cone_side, cone_bot);
		//glCallList(1);

		glfwSwapBuffers(WINDOW);
		Sleep(20);
	}
}

void save()
{
	FILE *file;
	file = fopen("save.txt", "w");
	fprintf(file, "%f\n", alpha);
	fprintf(file, "%f\n", beta);
	fprintf(file, "%f\n", gamma);
	fprintf(file, "%f\n", h);
	fprintf(file, "%f\n", r);
	fprintf(file, "%d\n", freq_x);
	fprintf(file, "%d\n", freq_y);
	fprintf(file, "%d\n", type);
	fclose(file);

	printf("-> current state saved\n");
}

void load()
{
	FILE *file;
	file = fopen("save.txt", "r");

	if (!file) {
		printf("-> failed to load!\n");
		return;
	}

	fscanf(file, "%f\n", &alpha);
	fscanf(file, "%f\n", &beta);
	fscanf(file, "%f\n", &gamma);
	fscanf(file, "%f\n", &h);
	fscanf(file, "%f\n", &r);
	fscanf(file, "%d\n", &freq_x);
	fscanf(file, "%d\n", &freq_y);
	fscanf(file, "%d\n", &type);
	fclose(file);

	printf("-> previous state loaded\n");
}

void initial_state()
{
	shift = true;
	tex_flag = true;
	twinkie_flag = false;
	light_flag = true;

	t = 0;
	freq_x = 15;
	freq_y = 3;
	type = -1;
	h = 3;
	r = 1;
	alpha = 0;
	beta = 0;
	gamma = 0;
	light_pos_x = 0.0;
	light_pos_y = 2.0;
	light_pos_z = 5.0;

	load_textures();

	count_cone();
	count_cylinder();

#ifdef X_LIGHTDEC
	light_flag = true;
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	light_enable();
#endif

#ifdef X_TEXSIZE
	tex_flag = true;
	glEnable(GL_TEXTURE_2D);
#endif

}

void controls(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) glfwSetWindowShouldClose(window, GL_TRUE);
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) { shift = !shift; printf("switched\n"); }
	if (key == GLFW_KEY_Z && action == GLFW_PRESS) initial_state();
	if (GLFW_KEY_UP == key && action == GLFW_PRESS) {
		if (h <= MAXHEIGHT) { h += 2; freq_y++; count_cone(); count_cylinder(); }
	}
	if (GLFW_KEY_DOWN == key && action == GLFW_PRESS) {
		if (h > MINHEIGHT) { h -= 2; freq_y--; count_cone(); count_cylinder(); }
	}
	if (GLFW_KEY_RIGHT == key && action == GLFW_PRESS) {
		if (r < MAXRAD) { r++; count_cone(); count_cylinder(); }
	}
	if (GLFW_KEY_LEFT == key && action == GLFW_PRESS) {
		if (r > MINRAD) { r--; count_cone(); count_cylinder(); }
	}
	if (GLFW_KEY_R == key && action == GLFW_PRESS) {
		if (freq_x <= MAXFREQ_X) { freq_x++; count_cone(); count_cylinder(); }
	}
	if (GLFW_KEY_F == key && action == GLFW_PRESS) {
		if (freq_x > MINFREQ_X) { freq_x--; count_cone(); count_cylinder(); }
	}
	if (GLFW_KEY_V == key && action == GLFW_PRESS) {
		type = -type;
			if (type == 1) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			if (type == -1) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	if (GLFW_KEY_T == key && action == GLFW_PRESS) { twinkie_flag = !twinkie_flag; tweenking(); }
	if (GLFW_KEY_Y == key && action == GLFW_PRESS) {
		tex_flag = !tex_flag;
		if (tex_flag) { glEnable(GL_TEXTURE_2D);} else { glDisable(GL_TEXTURE_2D); }
	}
	if (GLFW_KEY_H == key && action == GLFW_PRESS) { save(); }
	if (GLFW_KEY_G == key && action == GLFW_PRESS) { load(); count_cone(); count_cylinder(); }
	if (GLFW_KEY_L == key && action == GLFW_PRESS) {
		light_flag = !light_flag;
		if (light_flag) {
			glEnable(GL_LIGHTING);
			glEnable(GL_LIGHT0);
			light_enable();
		}
		else {
			glDisable(GL_LIGHTING);
			glDisable(GL_LIGHT0);
		}
	}

}

GLfloat cabinet_view_matrix[] = {
	1, 0, 0, 0,
	0, 1, 0, 0,
	-0.5*cos(M_PI / 6), -0.5*sin(M_PI / 6), -1, 0,
	0, 0, 0, 1
};

int test, fps;
FILE *X_FILE;

void display(GLFWwindow* window)
{
	printf("[H]: save state\n[G]: load state\n[T]: tweening aninmation\n[Y]: textures\n[L]: lighting\n[Q E W A S D]: rotation\n[Spacebar]:switch model/light rotation\n[Z]: reset state to default\n");
	
	initial_state();

	fps = 0;

	LARGE_INTEGER timerFrequency, timerStart, timerStop;
	QueryPerformanceFrequency(&timerFrequency);
	QueryPerformanceCounter(&timerStart);

#ifdef X_DISPLIST
	draw_figure(cone_top, cone_side, cone_bot);
#endif

	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.4, 0.4, 0.4, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS && shift) alpha -= ANGLE;
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS && shift) alpha += ANGLE;
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS && shift) beta -= ANGLE;
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS && shift) beta += ANGLE;
		if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS && shift) gamma -= ANGLE;
		if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS && shift) gamma += ANGLE;
		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS && !shift) { light_pos_y++; light_enable(); Sleep(120); }
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS && !shift) { light_pos_y--; light_enable(); Sleep(120); }
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS && !shift) { light_pos_x++; light_enable(); Sleep(120); }
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS && !shift) { light_pos_x--; light_enable(); Sleep(120); }
		if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS && !shift) { light_pos_z++; light_enable(); Sleep(120); }
		if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS && !shift) { light_pos_z--; light_enable(); Sleep(120); }

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glLoadMatrixf(cabinet_view_matrix);
		glOrtho(-10, 10, -10, 10, 10, -10);

#ifdef X_DISPLIST
		glCallList(list);
#endif

#ifndef X_DISPLIST
		draw_figure(cone_top, cone_side, cone_bot);
#endif
		fps++;

		QueryPerformanceCounter(&timerStop);
		double const t(static_cast<double>(timerStop.QuadPart - timerStart.QuadPart) / timerFrequency.QuadPart);

#ifdef PERFORM_TESTS
		if (t >= 1.0) {
			if (fps < 100) printf("are you perfroming tests on conole?\ndrakeface");
			fprintf(X_FILE, "%d\n", fps);
			return;
		}
#endif
		glfwSwapBuffers(window);
		//glfwWaitEvents();
		glfwPollEvents();
	}
}



int main(int argc, char** argv)
{

#ifdef PERFORM_TESTS
	X_FILE = fopen("x_file.txt", "w");
	fps = 0;

	for (test = 0; test < 10; test++) {
#endif

		if (!glfwInit())
		{
			printf("glfwInit failed\n");
			return -1;
		}

	#ifdef X_NOSMOOTH
		glfwWindowHint(GLFW_SAMPLES, CONSOLE_MODE);
	#endif

	#ifndef X_NOSMOOTH
		glfwWindowHint(GLFW_SAMPLES, PC_MODE);
	#endif

		glfwSetErrorCallback(error_callback);

	#ifdef X_SCREENSIZE
		GLFWwindow* window = glfwCreateWindow(850, 850, "Test app", NULL, NULL);
	#endif

	#ifndef X_SCREENSIZE
		GLFWwindow* window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Test app", NULL, NULL);
	#endif


			WINDOW = window;
			if (window == NULL)
			{
				printf("glfwOpenWindow failed.\n");
				glfwTerminate();
				return -2;
			}

			int attrib;
			attrib = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
			attrib = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
			attrib = glfwGetWindowAttrib(window, GLFW_OPENGL_PROFILE);

			glfwMakeContextCurrent(window);
			glfwSetKeyCallback(window, controls);

			//glDepthFunc(GL_LEQUAL);
			glEnable(GL_DEPTH_TEST);

	#ifdef X_CULLFACE
			glEnable(GL_CULL_FACE);
	#endif

			if (NULL != window)
			{
				display(window);
			}

			glfwDestroyWindow(window);
			glfwTerminate();

#ifdef PERFORM_TESTS
	}

	fclose(X_FILE);
	getchar();
#endif

	return 0;
}

