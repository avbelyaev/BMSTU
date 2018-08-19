//фронтальная диметрия (cabinet view)
//куб и усеченный конус
//кнопки: wasd q e r f v стрелочки

// glfw_30.cpp : Defines the entry point for the console application.
//  http://www.glfw.org/docs/latest/quick.html

#include "stdafx.h"
#include <time.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <vector>

#define SCREEN_WIDTH 1000
#define SCREEN_HEIGHT 1000
#define PC_MODE 16
#define CONSOLE_MODE 0

#define MAXRAD 4
#define MINRAD 1
#define MAXHEIGHT 7
#define MINHEIGHT 1
#define MAXFREQ_X 40
#define MINFREQ_X 3

#define ANGLE 2

typedef struct
{
	float x, y, z;
} coord;

std::vector<std::vector<coord>> coord_top;
std::vector<std::vector<coord>> coord_side;

int type, stop_i, stop_j, freq_x, freq_y;
static float alpha, beta, gamma, a, b, c, h, r;

static void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}

void draw_cube(
	int type,
	int size,
	int angle_x, int angle_y, int angle_z
	)
{
	glRotatef(angle_x, 1, 0, 0);
	glRotatef(angle_y, 0, 1, 0);
	glRotatef(angle_z, 0, 0, 1);
	glBegin(GL_QUADS);
	//BOTTOM
	glColor3f(0, 0, 0);
	glVertex3f(-size, -size, size);
	glVertex3f(size, -size, size);
	glVertex3f(size, -size, -size);
	glVertex3f(-size, -size, -size);
	glEnd();
	glBegin(GL_QUADS);
	//TOP
	glColor3f(0, 0, 1);
	glVertex3f(-size, size, size);
	glColor3f(1, 0, 0);
	glVertex3f(size, size, size);
	glColor3f(0, 0, 0);
	glVertex3f(size, size, -size);
	glColor3f(0, 0, 1);
	glVertex3f(-size, size, -size);
	glEnd();
	glBegin(GL_QUADS);
	//BACK
	glColor3f(1, 0, 1);
	glVertex3f(-size, -size, -size);
	glVertex3f(size, -size, -size);
	glVertex3f(size, size, -size);
	glVertex3f(-size, size, -size);
	glEnd();
	glBegin(GL_QUADS);
	//FRONT
	glColor3f(0, 0, 1);
	glVertex3f(-size, -size, size);
	glColor3f(0, 1, 0);
	glVertex3f(size, -size, size);
	glColor3f(1, 0, 0);
	glVertex3f(size, size, size);
	glColor3f(0, 0, 1);
	glVertex3f(-size, size, size);
	glEnd();
	glBegin(GL_QUADS);
	//LEFT
	glColor3f(1, 1, 0);
	glVertex3f(-size, -size, size);
	glVertex3f(-size, -size, -size);
	glVertex3f(-size, size, -size);
	glVertex3f(-size, size, size);
	glEnd();
	glBegin(GL_QUADS);
	//RIGHT
	glColor3f(0, 1, 0);
	glVertex3f(size, -size, size);
	glColor3f(0, 0, 1);
	glVertex3f(size, -size, -size);
	glColor3f(0, 0, 0);
	glVertex3f(size, size, -size);
	glColor3f(1, 0, 0);
	glVertex3f(size, size, size);
	glEnd();
}

void count_vertex(int freq_x, int freq_y)
{
	int i = 0, j = 0, top_rad = r, bot_rad = r * 2;
	float angle = 0, dangle, height, dheight, drad;

	height = h;
	dangle = 2 * M_PI / freq_x;
	dheight = height / freq_y;
	drad = (bot_rad - top_rad) / (float)freq_y;

	//count size of array
	while (angle <= 2 * M_PI + dangle) {
		angle += dangle;
		j++;
	}
	stop_j = j;

	//allocate memory for top_ side_ vertices
	coord_top.resize(freq_y + 1);
	for (i = 0; i < freq_y + 1; i++) coord_top[i].resize(stop_j + 1);

	coord_side.resize(freq_y + 1);
	for (i = 0; i < freq_y + 1; i++) coord_side[i].resize(stop_j + 1);

	//count vert
	for (i = 0; i < freq_y; i += 1) {
		for (j = 0, angle = 0; angle <= 2 * M_PI + dangle; j += 1, angle += dangle) {
			coord_top[i][j].x = (0 + drad*i)*cos(angle);
			coord_top[i][j].z = (0 + drad*i)*sin(angle);

			coord_top[i][j + 1].x = (0 + drad*i)*cos(angle + dangle);
			coord_top[i][j + 1].z = (0 + drad*i)*sin(angle + dangle);

			coord_top[i + 1][j + 1].x = (0 + drad*(i + 1))*cos(angle + dangle);
			coord_top[i + 1][j + 1].z = (0 + drad*(i + 1))*sin(angle + dangle);

			coord_top[i + 1][j].x = (0 + drad*(i + 1))*cos(angle);
			coord_top[i + 1][j].z = (0 + drad*(i + 1))*sin(angle);


			coord_side[i][j].x = (top_rad + drad*i)*cos(angle);
			coord_side[i][j].y = height - dheight*i;
			coord_side[i][j].z = (top_rad + drad*i)*sin(angle);

			coord_side[i][j + 1].x = (top_rad + drad*i)*cos(angle + dangle);
			coord_side[i][j + 1].y = height - dheight*i;
			coord_side[i][j + 1].z = (top_rad + drad*i)*sin(angle + dangle);

			coord_side[i + 1][j + 1].x = (top_rad + drad*(i + 1))*cos(angle + dangle);
			coord_side[i + 1][j + 1].y = height - dheight*(i + 1);
			coord_side[i + 1][j + 1].z = (top_rad + drad*(i + 1))*sin(angle + dangle);

			coord_side[i + 1][j].x = (top_rad + drad*(i + 1))*cos(angle);
			coord_side[i + 1][j].y = height - dheight*(i + 1);
			coord_side[i + 1][j].z = (top_rad + drad*(i + 1))*sin(angle);
		}
	}
	stop_i = i;
	stop_j = j;

	printf("coordinates renewed: h=%.2f, top_r=%d, fr_x=%d, fr_y=%d\n", height, top_rad, freq_x, freq_y);
}

void draw_cone(
	int type,
	int top_rad, int bot_rad,
	int height,
	int freq_x, int freq_y,
	int angle_x, int angle_y, int angle_z
	)
{
	int i, j;

	glRotatef(angle_x, 1, 0, 0);
	glRotatef(angle_y, 0, 1, 0);
	glRotatef(angle_z, 0, 0, 1);

	for (i = 0; i < stop_i; i += 1) {
		for (j = 0; j < stop_j; j += 1) {
			//top
			type > 0 ? glBegin(GL_POLYGON) : glBegin(GL_LINE_LOOP);
			glColor3f(0.5, 0.5, 0.5);
			glVertex3f(coord_top[i][j].x, height, coord_top[i][j].z);
			glColor3f(0, 0, 0);
			glVertex3f(coord_top[i][j + 1].x, height, coord_top[i][j + 1].z);
			glVertex3f(coord_top[i + 1][j + 1].x, height, coord_top[i + 1][j + 1].z);
			glVertex3f(coord_top[i + 1][j].x, height, coord_top[i + 1][j].z);
			glEnd();


			//bot
			type > 0 ? glBegin(GL_POLYGON) : glBegin(GL_LINE_LOOP);
			glColor3f(0.3, 0.3, 0.3);
			glVertex3f(coord_top[i][j].x, 0, coord_top[i][j].z);
			glColor3f(0, 0, 0);
			glVertex3f(coord_top[i][j + 1].x, 0, coord_top[i][j + 1].z);
			glVertex3f(coord_top[i + 1][j + 1].x, 0, coord_top[i + 1][j + 1].z);
			glVertex3f(coord_top[i + 1][j].x, 0, coord_top[i + 1][j].z);
			glEnd();

			type > 0 ? glBegin(GL_POLYGON) : glBegin(GL_LINE_LOOP);
			glColor3f(0.3, 0.3, 0.3);
			glVertex3f(coord_side[i][j].x, 0, coord_side[i][j].z);
			glColor3f(0, 0, 0);
			glVertex3f(coord_side[i][j + 1].x, 0, coord_side[i][j + 1].z);
			glVertex3f(coord_side[i + 1][j + 1].x, 0, coord_side[i + 1][j + 1].z);
			glVertex3f(coord_side[i + 1][j].x, 0, coord_side[i + 1][j].z);
			glEnd();


			//sidelines
			type > 0 ? glBegin(GL_POLYGON) : glBegin(GL_LINE_LOOP);
			glColor3f(0, 0, 0);
			glVertex3f(coord_side[i][j].x, coord_side[i][j].y, coord_side[i][j].z);
			glColor3f(0.8, 0, 0);
			glVertex3f(coord_side[i][j + 1].x, coord_side[i][j + 1].y, coord_side[i][j + 1].z);
			glColor3f(0, 0.8, 0);
			glVertex3f(coord_side[i + 1][j + 1].x, coord_side[i + 1][j + 1].y, coord_side[i + 1][j + 1].z);
			glColor3f(0, 0, 0.8);
			glVertex3f(coord_side[i + 1][j].x, coord_side[i + 1][j].y, coord_side[i + 1][j].z);
			glEnd();
		}
	}
}

void draw_axis(
	int len_x_start, int len_x_stop,
	int len_y_start, int len_y_stop,
	int len_z_start, int len_z_stop
	)
{
	glBegin(GL_LINES);
	//X
	glColor3f(1, 0, 0);
	glVertex3f(len_x_start, 0, 0);
	glVertex3f(len_x_stop, 0, 0);
	//Y
	glColor3f(0, 1, 0);
	glVertex3f(0, len_y_start, 0);
	glVertex3f(0, len_y_stop, 0);
	//Z
	glColor3f(0, 0, 1);
	glVertex3f(0, 0, len_z_start);
	glVertex3f(0, 0, len_z_stop);
	glEnd();
}

void draw_net()
{
	int i;
	float j;

	for (i = 0; i < 20; i++) {
		glColor3f(0, 0.7, 0);
		for (j = -10; j <= 10; j++) {
			glBegin(GL_LINES);
			glVertex3f(j, -3, 10);
			glVertex3f(j, -3, -10);
			glEnd();
		}
		for (j = -10; j <= 10; j += 2) {
			glBegin(GL_LINES);
			glVertex3f(10, -3, j);
			glVertex3f(-10, -3, j);
			glEnd();
		}
	}
}

void controls(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	switch (key) {
	case (GLFW_KEY_SPACE) :
		if (action == GLFW_PRESS) {
		alpha = 0;
		beta = 0;
		gamma = 0;
		h = 2;
		r = 1;
		freq_x = 3;
		freq_y = 3;
		count_vertex(freq_x, freq_y);
		}
		break;
	case (GLFW_KEY_UP) :
		if (action == GLFW_PRESS) {
			if (h <= MAXHEIGHT) {
				h += 2;
				freq_y++;
				count_vertex(freq_x, freq_y);
			}
		}
		break;
	case (GLFW_KEY_DOWN) :
		if (action == GLFW_PRESS) {
			if (h > MINHEIGHT) {
				h -= 2;
				freq_y--;
				count_vertex(freq_x, freq_y);
			}
		}
		break;
	case (GLFW_KEY_RIGHT) :
		if (action == GLFW_PRESS) {
			if (r < MAXRAD) {
				r++;
				count_vertex(freq_x, freq_y);
			}
		}
		break;
	case (GLFW_KEY_LEFT) :
		if (action == GLFW_PRESS) {
			if (r > MINRAD) {
				r--;
				count_vertex(freq_x, freq_y);
			}
		}
		break;
	case (GLFW_KEY_R) :
		if (action == GLFW_PRESS) {
			if (freq_x <= MAXFREQ_X) {
				count_vertex(++freq_x, freq_y);
			}
		}
		break;
	case (GLFW_KEY_F) :
		if (action == GLFW_PRESS) {
			if (freq_x > MINFREQ_X) {
				count_vertex(--freq_x, freq_y);
			}
		}
		break;
	case (GLFW_KEY_V) :
		if (action == GLFW_PRESS) {
			type = -type;
			if (type == 1) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			if (type == -1) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		break;
	}

}

void display(GLFWwindow* window)
{
	int i;

	freq_x = 3;
	freq_y = 3;
	type = -1;
	h = 3;
	r = 1;
	alpha = 0;
	beta = 0;
	gamma = 0;

	GLfloat cabinet_view_matrix[] = {
		1, 0, 0, 0,
		0, 1, 0, 0,
		-0.5*cos(M_PI / 6), -0.5*sin(M_PI / 6), -1, 0,
		0, 0, 0, 1
	};      //0.5*cos(M_PI/6), 0.5*sin(M_PI/6), 0, 0,

	count_vertex(freq_x, freq_y);
	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.4, 0.4, 0.4, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) alpha -= ANGLE;
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) alpha += ANGLE;
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) beta -= ANGLE;
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) beta += ANGLE;
		if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) gamma -= ANGLE;
		if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) gamma += ANGLE;

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glLoadMatrixf(cabinet_view_matrix);

		glOrtho(-10, 10, -10, 10, 10, -10);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		draw_net();
		draw_axis(-100, 100, -3, 100, -1, 100);

		glTranslatef(-4, 0, 0);
		
		draw_cube(type, 2, 0, 0, 0);

		glTranslatef(8, 0, 0);
		
		//draw_cube(type, 2, alpha, beta, gamma);
		draw_cone(type, r, r * 2, h, freq_x, freq_y, alpha, beta, gamma);

		glfwSwapBuffers(window);
		glfwWaitEvents();
		//glfwPollEvents();
	}
	
	//free vectors of top_ side_ vertices
	for (i = 0; i < freq_y + 1; i++) {
		coord_top[i].clear();
		coord_top[i].shrink_to_fit();
	}
	coord_top.clear();
	coord_top.shrink_to_fit();

	for (i = 0; i < freq_y + 1; i++) {
		coord_side[i].clear();
		coord_side[i].shrink_to_fit();
	}
	coord_side.clear();
	coord_side.shrink_to_fit();

}

int main(int argc, char** argv)
{
	// initialise GLFW
	if (!glfwInit())
	{
		printf("glfwInit failed\n");
		return -1;
	}
	glfwWindowHint(GLFW_SAMPLES, PC_MODE);
	glfwSetErrorCallback(error_callback);

	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 1);
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
	GLFWwindow* window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Test app", NULL, NULL);

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

	if (NULL != window)
	{
		display(window);
	}
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}

