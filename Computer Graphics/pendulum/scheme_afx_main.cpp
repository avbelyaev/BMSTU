//scheme.txt
5
c
1001 90 380 1000
2001
b
2001 200 700 0
1001 0 2002
b
2002 300 480 10
2001 1002 0
c
1002 400 250 300
2002
p
800 930 600 250
end_scheme


1
p
500 930 600 700
end_scheme



3
c
1001 100 450 100
2000
b
2000 200 800 0
1001 0 1002
c
1002 300 600 10
2000
end_scheme



5
b
2001 150 480 10
0 1001 2002
c
1001 120 250 30
2001
b
2002 250 700 0
2001 0 1002
c
1002 310 280 100
2002
p
800 930 600 250
end_scheme




7
c
1001 90 380 100
2001
b
2001 200 750 0
1001 0 2002
b
2002 300 480 10
2001 1002 2003
c
1002 350 210 300
2002
b
2003 400 600 0
2002 0 1003
c
1003 460 370 1000
2003
p
800 930 600 250
end_scheme




10
c
1001 210 400 100
2001
b
2001 250 700 0
1001 0 2002
b
2002 350 500 5
2001 1002 2003
c
1002 320 230 30
2002
b
2003 450 700 0
2002 0 2004
b
2004 550 500 10
2003 1003 2005
c
1003 570 300 50
2004
b
2005 650 700 0
2004 0 1004
c
1004 760 400 10
2005
p
800 930 600 250
end_scheme



//=======================================================================================
//=======================================================================================
//=======================================================================================
//stdafx.h
// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//
#pragma once

#include "targetver.h"

#define _CRT_SECURE_NO_WARNINGS

//if print is lagging switch this
#include <stdio.h>		
//#include <iostream>
#include <vector>
#include <tchar.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

#include <ctime>
#include <time.h>

#include <Windows.h>

#include <string> 

#define GLEW_STATIC
#include <GL\glew.h>
#include <GLFW\glfw3.h>

typedef struct { float x, y; } coord;


#include "cell.h"
#include "map.h"
#include "pendulum.h"
#include "node.h"
#include "dll.h"



#define SCREEN_WIDTH 1000
#define SCREEN_HEIGHT 1000
#define ROOF_HEIGHT 950
#define FLOOR_HEIGHT 50
#define SPLIT_X 10
#define SPLIT_Y 10
#define G_CONST 9.781665
#define MAXRANDID 99999

#define REVD(COORD_Y) ((SPLIT_Y - 1) - (COORD_Y))

#define PC_MODE 16
#define CONSOLE_MODE 0


//Performance tests:
//#define PERFORM_TEST


//=======================================================================================
//=======================================================================================
//=======================================================================================
//glfw30.cpp


// glfw_30.cpp : Defines the entry point for the console application.
//  http://www.glfw.org/docs/latest/quick.html

#include "stdafx.h"

#include <GLAux.h>
#pragma comment(lib,"glaux.lib")


using namespace std;

GLdouble A, B, C, D;

float t_main = 0;
float t_list_pend = 0, t_list_move = 0, t_pend = 0;
static float alpha;

bool pendulum_exists = false;

int mouse_curr_x, mouse_curr_y;

class cell;
class map;
class node;
class dll;
class pendulum;

map *m;
dll *list;
pendulum *pend;

#ifdef PERFORM_TEST
float test_fps;
float test_time = 0;
int test_n = 0;
FILE *TEST_FPS_FILE;
#endif



void draw_constraints(int h) {
	
	glLineWidth(2.0);
	glColor3f(0, 0, 0);
	int i, n = 20;
	float di = SCREEN_WIDTH / n;
	
	if (h > SCREEN_HEIGHT / 2) {
		glBegin(GL_POLYGON);
			glColor3f(1, 1, 1);
			glVertex2f(0, h);
			glVertex2f(0, SCREEN_HEIGHT);
			glVertex2f(SCREEN_WIDTH, SCREEN_HEIGHT);
			glVertex2f(SCREEN_WIDTH, h);
		glEnd();
		for (i = 0; i < SCREEN_WIDTH; i += (int)di) {
			glBegin(GL_LINES);
				glColor3f(0,0,0);
				glVertex2f(i, h);
				glVertex2f(i + di, h + di);
			glEnd();
		}
	}
	else {
		glBegin(GL_POLYGON);
			glColor3f(1, 1, 1);
			glVertex2f(0, -h);
			glVertex2f(0, h);
			glVertex2f(SCREEN_WIDTH, h);
			glVertex2f(SCREEN_WIDTH, -h);
		glEnd();
		for (i = 0; i < SCREEN_WIDTH; i += (int)di) {
			glBegin(GL_LINES);
				glColor3f(0,0,0);
				glVertex2f(i, h);
				glVertex2f(i - di, h - di);
			glEnd();
		}
	}

	glLineWidth(4.0);
	glBegin(GL_LINES);
		glColor3f(0,0,0);
		glVertex2f(0, h);
		glVertex2f((GLfloat)SCREEN_WIDTH, h);
	glEnd();
}

void draw_gears(float start_angle, float angle, int x, int y, int inner_rad, int outer_rad, int teeth_len, int teeth, float color) {
	
	glPushMatrix();
	glRotatef(angle + start_angle, 0, 0, 1);

	int i;
	GLfloat phi;
	float da = 2 * M_PI / teeth / 4.0;

	glBegin(GL_QUAD_STRIP);
	glColor3f(color, color, color);
	
	for (i = 0; i <= teeth; i++) {
		phi = i * 2.0 * M_PI / teeth;
		glVertex2f(x + inner_rad * cos(phi), y + inner_rad * sin(phi));
		glVertex2f(x + outer_rad * cos(phi), y + outer_rad * sin(phi));
		glVertex2f(x + inner_rad * cos(phi), y + inner_rad * sin(phi));
		glVertex2f(x + outer_rad * cos(phi + 3 * da), y + outer_rad * sin(phi + 3 * da));
	}
	glEnd();

	glBegin(GL_QUADS);
	da = 2.0 * M_PI / teeth / 4.0;
	for (i = 0; i < teeth; i++) {
		phi = i * 2.0 * M_PI / teeth;

		glVertex2f(x + outer_rad * cos(phi), y + outer_rad * sin(phi));
		glVertex2f(x + (outer_rad + teeth_len) * cos(phi + da), y + (outer_rad + teeth_len) * sin(phi + da));
		glVertex2f(x + (outer_rad + teeth_len) * cos(phi + 2 * da), y + (outer_rad + teeth_len) * sin(phi + 2 * da));
		glVertex2f(x + outer_rad * cos(phi + 3 * da), y + outer_rad * sin(phi + 3 * da));
	}
	glEnd();

	glPopMatrix();

}

void draw_permission_signal() {

	GLfloat phi = 0, theta = 0;

	int center_x = 50, center_y = ROOF_HEIGHT - 50;
	int rad = 20;

	//outer
	glBegin(GL_POLYGON);

	glColor3f(0, 0, 0);
	for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
		glVertex2f(center_x + rad*cos(phi), center_y + rad*sin(phi));
	}
	glEnd();

	//inner
	rad = 18;
	glBegin(GL_POLYGON);
	glColor3f(1, 0.2, 0.2);
	for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
		glVertex2f(center_x + rad*cos(phi), center_y + rad*sin(phi));
	}
	glEnd();

	rad = 14;
	glBegin(GL_POLYGON);
	glColor3f(1, 0, 0);
	for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
		glVertex2f(center_x + rad*cos(phi), center_y + rad*sin(phi));
	}
	glEnd();

}

static void cursor_callback(GLFWwindow* window, double x, double y)
{
	mouse_curr_x = (int)x;
	mouse_curr_y = (int)y;

	int tmpx = (int)(x / (SCREEN_WIDTH / SPLIT_X));
	int tmpy = (int)(y / (SCREEN_HEIGHT / SPLIT_Y));

	//printf("%d\t%d\t%c\n", (int)x, (int)y, m->get_cell_type(tmpx, tmpy));

	//m->dispell_map();
	//m->highlight_cell(tmpx, tmpy);
}

static void mouse_callback(GLFWwindow* window, int button, int action, int mods)
{
	int tmpx = (int)(mouse_curr_x / (SCREEN_WIDTH / SPLIT_X));
	int tmpy = (int)(mouse_curr_y / (SCREEN_HEIGHT / SPLIT_Y));

	if (button == GLFW_MOUSE_BUTTON_RIGHT)
	{
		if (action == GLFW_PRESS) {
			

			//printf("%d\t%d\t%c\n", (int)mouse_curr_x, (int)mouse_curr_y, m->get_cell_type(tmpx, tmpy));
		}
		if (action == GLFW_RELEASE) glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	}

	if (button == GLFW_MOUSE_BUTTON_LEFT)
	{
		if (action == GLFW_PRESS) {

			//printf("M: %d\t%d\t%c(before click)\n", mouse_curr_x, mouse_curr_y, m->get_cell_type(tmpx, tmpy));
			/*
			if ('e' == m->get_cell_type(tmpx, tmpy)) {
				m->erase_cell(tmpx, tmpy);
				m->set_cell(1, tmpx, tmpy, 'b');
				return;
			}
			if ('b' == m->get_cell_type(tmpx, tmpy)) {
				m->erase_cell(tmpx, tmpy);
				m->set_cell(1, tmpx, tmpy, 'c');
				return;
			}
			if ('c' == m->get_cell_type(tmpx, tmpy)) {
				m->erase_cell(tmpx, tmpy);
				m->set_cell(1, tmpx, tmpy, 'e');
				return;
			}
			*/
			
		}
		
	}
}

static void resize_callback(GLFWwindow* window, int width, int height)
{
	//      windowWidth = width;
	//      windowHeight = height;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0, (GLdouble)width, 0.0, (GLdouble)height, -1, 1);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	A = width / 4.0;
	B = 0.0;
	C = D = height / 2.0;

	//printf("Reshape occured\n");
}

static void keyboard_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);
	if (GLFW_KEY_S == key && action == GLFW_PRESS) { 
		m->validate_map();
	}
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
		//scheme
		if (list->permit_phys_move) {
			if (!list->phys_move_enabled) {
				list->phys_move_enabled = true;
			}
			else {
				list->phys_move_enabled = false;
			}
		}
		else {
			printf("physics can't be enabled at the moment. check back later.\n");
		}
		
	}
}

static void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}


void display(GLFWwindow* window) 
{
	float LO = 0.2;
	float HI = 0.8;

	float clr1 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
	float clr2 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
	float clr3 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
	float clr4 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
	float clr5 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));

	alpha = 0;


#ifdef PERFORM_TEST

	double lastTime = glfwGetTime();
	int nbFrames = 0;
	float fps_sum = 0;

	
	printf("Test#%d\n", test_fps);
#endif

	while (!glfwWindowShouldClose(window))
	{
		glClearColor(0.6, 0.6, 0.6, 1);
		//glClearColor(1, 1, 1, 1);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		draw_gears(0, alpha*0.7, -1500, 1500, 200, 1500, 20, 300, 0.8);
		draw_gears(0, -alpha*0.8, 200, 900, 200, 500, 30, 30, clr1);
		draw_gears(0, alpha * 0.7, -200, 900, 60, 480, 30, 100, clr2);
		draw_gears(0, -alpha, 260, 120, 40, 90, 30, 10, clr3);
		draw_gears(400, alpha*1.2, 850, 330, 50, 200, 50, 17, clr4);
		//draw_gears(clr4*100, alpha*0.3, 850, 330, 50, 200, clr4*10, 17, clr4);
		draw_gears(500, alpha*1.1, 850, 800, 60, 200, 40, 20, clr5);
		draw_gears(clr5 * 100, alpha*clr5, 850, 800, clr4 * 10, 200, 40, clr5 * 10, clr5);
		
		draw_constraints(ROOF_HEIGHT);
		draw_constraints(FLOOR_HEIGHT);
		
		if (!list->permit_phys_move && cos(t_main) > 0) {
			draw_permission_signal();
		}

		if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
			m->print_map();
		}

		//m->draw_map();
		if (pendulum_exists) {
			pend->draw();
			pend->move(t_pend);
			t_pend += 0.05;
		}

		list->draw_objects();
		list->draw_ropes();

		float dt_list_move = 0.005;

		list->check_attenuation();
		list->enable_phys_pend(t_list_pend);

		list->check_collisions();
		list->enable_phys_move(t_list_move, dt_list_move);

		if (list->phys_move_enabled && list->permit_phys_move) {
			t_list_move += dt_list_move;
		}
		if (list->phys_pend_enabled) {
			t_list_pend += 0.06;
		}

		t_main += 0.02;
		alpha += 0.03;
		
#ifdef PERFORM_TEST

		double currentTime = glfwGetTime();
		nbFrames++;
		if (currentTime - lastTime >= 1.0) { 
			printf("%f fps\n", double(nbFrames));
			fps_sum += nbFrames;
			nbFrames = 0;
			lastTime += 1.0;
		}

		if (currentTime >= 10.0) {
			fprintf(TEST_FPS_FILE, "%f\n", fps_sum / 10.0);
			return;
		}

#endif
		
		glfwSwapBuffers(window);
		glfwPollEvents();
		//glfwWaitEvents();
	}

}


int main(int argc, _TCHAR* argv[])
{
	int num;

#ifdef PERFORM_TEST
	TEST_FPS_FILE = fopen("FPS_TEST_FILE.txt", "w");

	for (test_fps = 0; test_fps < 10; test_fps++)
	{

#endif

		if (!glfwInit()) {
			printf("glfwInit failed\n");
			return -1;
		}

		glfwWindowHint(GLFW_SAMPLES, PC_MODE);

		glfwSetErrorCallback(error_callback);

		GLFWwindow* window;

		//      glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 1);
		//      glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
		//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
		window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Anthony's Porject", NULL, NULL);
		if (window == NULL) {
			printf("glfwOpenWindow failed.\n");
			glfwTerminate();
			return -2;
		}

		int i, attrib;
		attrib = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
		attrib = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
		attrib = glfwGetWindowAttrib(window, GLFW_OPENGL_PROFILE);

		glfwMakeContextCurrent(window);

		glfwSetKeyCallback(window, keyboard_callback);
		glfwSetFramebufferSizeCallback(window, resize_callback);
		glfwSetMouseButtonCallback(window, mouse_callback);
		glfwSetCursorPosCallback(window, cursor_callback);
		resize_callback(window, SCREEN_WIDTH, SCREEN_HEIGHT);

		srand(time(NULL));

		FILE *fp;
		fp = fopen("scheme.txt", "r");
		if (NULL == fp) {
			printf("cannot open file\n");
		}


		list = new dll();


		pend = new pendulum();


		m = new map();
		m->set_map(SPLIT_X, SPLIT_Y);


		printf("reading scheme\n");


		
		fscanf(fp, "%d", &num);

		for (i = 0; i < num; i++) {
			char el_type[6];
			printf("type ");
			fscanf(fp, "%s", el_type);

			if (0 == strcmp(el_type, "end_scheme"))
				break;

			node *curr_node;
			int pos_x, pos_y, mass;
			int body_name, conn_left_name, conn_center_name, conn_right_name;

			if (0 == strcmp(el_type, "block") || 0 == strcmp(el_type, "b")) {

				printf("name pos_x pos_y mass ");
				fscanf(fp, "%d%d%d%d", &body_name, &pos_x, &pos_y, &mass);
				printf("connection[L C R] ");
				fscanf(fp, "%d%d%d", &conn_left_name, &conn_center_name, &conn_right_name);

				list->append_back(body_name);
				curr_node = list->get_by_id(body_name);
				curr_node->set_list_param('b', body_name, i, pos_x, pos_y, mass);
				curr_node->set_cp_ids(conn_left_name, conn_center_name, conn_right_name);

				printf("block:%d[x:%d y:%d m:%d CL:%d CC:%d CR:%d]\n\n", body_name, pos_x, pos_y, mass, conn_left_name, conn_center_name, conn_right_name);
			}

			if (0 == strcmp(el_type, "cargo") || 0 == strcmp(el_type, "c")) {

				printf("name pos_x pos_y mass ");
				fscanf(fp, "%d%d%d%d", &body_name, &pos_x, &pos_y, &mass);
				printf("connection[T B] ");
				fscanf(fp, "%d", &conn_center_name);

				list->append_back(body_name);
				curr_node = list->get_by_id(body_name);
				curr_node->set_list_param('c', body_name, i, pos_x, pos_y, mass);
				curr_node->set_cp_ids(conn_center_name);

				printf("cargo:%d[x:%d y:%d m:%d CC:%d]\n\n", body_name, pos_x, pos_y, mass, conn_center_name);
			}

			if (0 == strcmp(el_type, "pendulum") || 0 == strcmp(el_type, "p")) {
				pendulum_exists = true;

				float carrier_x, carrier_y, ticker_x, ticker_y;
				fscanf(fp, "%f%f%f%f", &carrier_x, &carrier_y, &ticker_x, &ticker_y);

				coord carr, tick;

				carr.x = carrier_x;
				carr.y = carrier_y;
				tick.x = ticker_x;
				tick.y = ticker_y;

				pend->set(carr, tick);

				printf("pendulum [CRR:%f,%f TKR:%f,%f] has been added\n\n", pend->carrier.x, pend->carrier.y, pend->curr_ticker.x, pend->curr_ticker.y);
			}
		}
		fclose(fp);
		printf("file has been closed\n");

		list->set_connection_points();
		list->set_pend_params();


		list->display_objects();

		printf("\n");

		if (NULL != window) {

			display(window);

		}

		glfwDestroyWindow(window);
		glfwTerminate();

#ifdef PERFORM_TEST
	}
	fprintf(TEST_FPS_FILE, "scheme difficulty:%d\n", num);

	fclose(TEST_FPS_FILE);
	getchar();
#endif

	return 0;
}

