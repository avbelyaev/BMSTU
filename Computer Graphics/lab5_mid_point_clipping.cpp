//отсечение средней точкой
//управление: задать 2 точки -> d -> рисовать отрезки

// glfw_30.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <time.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <GL/glew.h>
#include <GL\glut.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <cstdlib>

#define SCREEN_WIDTH		600
#define SCREEN_HEIGHT		400
#define	SCREEN_POS_X		100
#define	SCREEN_POS_Y		100

typedef struct point
{
	int x, y;
	bool p_visibility;

	//original idea by Ivan.M
	union {
		struct {
			char _0bit : 1;		//standalone bit
			char _1bit : 1;
			char _2bit : 1;
			char _3bit : 1;
			char _rest : 4;		//standalone 4bits
		};
		char half : 4;		//first half of char (first 4bits) above presented as number
	};
}point;

typedef struct line
{
	point a, b;
	int len;
	bool l_visibility;
}line;

typedef struct clipper
{
	point a, b, c, d;
	int x_min, x_max;
	int y_min, y_max;
}clipper;

bool draw_flag;

std::vector<point> temp_vect;

clipper my_clipper;

std::vector<line> lines;
std::vector<line> pre_visible;

std::vector<line> final_lines;

int line_len(point a, point b)
{
	int dx = abs(a.x - b.x);
	int dy = abs(a.y - b.y);

	return (int)sqrt(dx * dx + dy*dy);
}

int get_code(point a)
{
	int code = 0;

	a._0bit = 0;
	a._1bit = 0;
	a._2bit = 0;
	a._3bit = 0;
	a._rest = 0;

	if (a.x < my_clipper.x_min) {
		code += 1;
		a._0bit = 1;
	}
	if (a.y > my_clipper.y_max) {
		code += 10;
		a._1bit = 1;
	}
	if (a.x > my_clipper.x_max) {
		code += 100;
		a._2bit = 1;
	}
	if (a.y < my_clipper.y_min) {
		code += 1000;
		a._3bit = 1;
	}

	return code;
}

bool point_visibility(point a)
{
	bool inner_x = false, inner_y = false;

	if (my_clipper.x_max > a.x && my_clipper.x_min < a.x) inner_x = true;
	if (my_clipper.y_max > a.y && my_clipper.y_min < a.y) inner_y = true;

	if (inner_x && inner_y) return true;

	return false;
}

int line_visibility(point a, point b)
{
	int code_a = get_code(a);
	int code_b = get_code(b);

	if (a.half & b.half != 0) return -1;	//triv_invisibility

	if (a.half == 0 && b.half == 0) return 1;	//visibility

	return 0;
}

point mid_point(point a, point b)
{
	point middle;

	middle.x = (a.x + b.x) / 2;
	middle.y = (a.y + b.y) / 2;
	int code = get_code(middle);

	return middle;
}

void clip(line AB)
{
	int code_a = get_code(AB.a);
	int code_b = get_code(AB.b);
	int l_v = line_visibility(AB.a, AB.b);
	int len = line_len(AB.a, AB.b);

	line half_1, half_2;

	if (l_v == -1)
		return;

	if (l_v == 1) {
		pre_visible.push_back(AB);
		return;
	}

	if (len == 1) {
		if (point_visibility(AB.a)) 
			pre_visible.push_back(AB);
		return;
	}

	if (l_v == 0) {
		half_1.a = AB.a;
		half_1.b = mid_point(AB.a, AB.b);
		clip(half_1);

		half_2.a = mid_point(AB.a, AB.b);
		half_2.b = AB.b;
		clip(half_2);
	}
}

void changeSize(int w, int h) 
{
	if (0 == h) h = 1;

	float ratio = w * 1.0 / h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluOrtho2D(0, (GLdouble)w, 0, (GLdouble)h);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}


//==================INPUT===================
void set_clipper()
{
	int x_max, x_min, y_max, y_min;
	point tmp1 = temp_vect[0], tmp2 = temp_vect[1];

	if (tmp1.x > tmp2.x) {
		x_max = tmp1.x;
		x_min = tmp2.x;
	}
	else {
		x_max = tmp2.x;
		x_min = tmp1.x;
	}

	if (tmp1.y > tmp2.y) {
		y_max = tmp1.y;
		y_min = tmp2.y;
	}
	else {
		y_max = tmp2.y;
		y_min = tmp1.y;
	}

	my_clipper.x_max = x_max;
	my_clipper.x_min = x_min;
	my_clipper.y_max = y_max;
	my_clipper.y_min = y_min;

	my_clipper.a.x = x_min;
	my_clipper.a.y = y_max;

	my_clipper.b.x = x_max;
	my_clipper.b.y = y_max;

	my_clipper.c.x = x_max;
	my_clipper.c.y = y_min;

	my_clipper.d.x = x_min;
	my_clipper.d.y = y_min;
}

void draw_clipper()
{
	glBegin(GL_LINE_LOOP);
		glColor3f(1, 0, 0);
		glVertex2f(my_clipper.a.x, my_clipper.a.y);
		glVertex2f(my_clipper.b.x, my_clipper.b.y);
		glVertex2f(my_clipper.c.x, my_clipper.c.y);
		glVertex2f(my_clipper.d.x, my_clipper.d.y);
	glEnd();
}

void draw_original_lines()
{
	for (int i = 0; i < lines.size(); i++) {
		glColor3f(0, 0, 1);
		glBegin(GL_LINES);
			glVertex2f(lines[i].a.x, lines[i].a.y);
			glVertex2f(lines[i].b.x, lines[i].b.y);
		glEnd();
	}
}


//=================OUTPUT===================
void draw_visible_lines()
{
	glLineWidth(1.0);
	for (int i = 0; i < final_lines.size(); i++) {
		glBegin(GL_LINES);
			glColor3f(1, 0, 0);
			glVertex2f(final_lines[i].a.x, final_lines[i].a.y);
			glVertex2f(final_lines[i].b.x, final_lines[i].b.y);
		glEnd();
	}
}


//===================MAIN====================
void renderScene(void) 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (draw_flag && temp_vect.size() > 1) {
		draw_clipper();
		draw_original_lines();
		draw_visible_lines();
	}

	glutSwapBuffers();
}

int m = 0;

void mouseButton(int button, int state, int x, int y)
{
	if (GLUT_LEFT_BUTTON == button) {
		int k, p;
		point a;
		line AB, tmp;
		if (GLUT_DOWN == state) {
			a.x = x;
			a.y = SCREEN_HEIGHT - y;

			temp_vect.push_back(a);
			if (m > 1 && m % 2 != 0) {
				AB.a = temp_vect[m - 1];
				AB.b = temp_vect[m];

				lines.push_back(AB);

				clip(AB);	//main process - mid point segmentation clipping

				if (pre_visible.size() == 0) {
					goto fuck;
				} else {
					p = 0;
					k = pre_visible.size() - 1;
				}

				tmp.a = pre_visible[p].a;
				tmp.b = pre_visible[k].b;

				final_lines.push_back(tmp);

				pre_visible.clear();
				pre_visible.shrink_to_fit();
			}
			
		fuck:
			m++;
			glutPostRedisplay();
		}
		if (GLUT_UP == state) {
			glutPostRedisplay();
		}
	}
}

void keyboardButton(unsigned char key, int x, int y)
{
	if ('d' == key || 'D' == key) {
		draw_flag = true;
		set_clipper();

		printf("clipper\n");

		glutPostRedisplay();
	}
	if ('t' == key || 'T' == key) {
		printf("temp_vect.size=%d vertices\n", temp_vect.size());
		printf("lines.size=%d lines\n", lines.size());
		for (int i = 0; i < lines.size(); i++) {
			printf(" %dx%d -> ", lines[i].a.x, lines[i].a.y);
			printf("%dx%d\n", lines[i].b.x, lines[i].b.y);
		}
		printf("vis_lines.size=%d lines\n", pre_visible.size());
		for (int i = 0; i < pre_visible.size(); i++) {
			printf(" %dx%d -> ", pre_visible[i].a.x, pre_visible[i].a.y);
			printf("%dx%d\n", pre_visible[i].b.x, pre_visible[i].b.y);
		}
		printf("temp_vect: p_X x p_Y (code)\n");
		for (int i = 2; i < temp_vect.size(); i++) {
			printf("%dx%d (%d)\n", temp_vect[i].x, temp_vect[i].y, get_code(temp_vect[i]));
		}
		glutPostRedisplay();
	}
}

void initial_state()
{
	glClearColor(0.9, 0.9, 0.9, 1);

	draw_flag = false;
}

int main(int argc, char **argv) 
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(SCREEN_POS_X, SCREEN_POS_Y);
	glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);
	glutCreateWindow("Anthony's Project: Lab 5(clipping)");

	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutIdleFunc(renderScene);

	glutMouseFunc(mouseButton);
	glutKeyboardFunc(keyboardButton);

	initial_state();

	glutMainLoop();

	return 0;
}

