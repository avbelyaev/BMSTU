//управление: стрелочки

// glfw_30.cpp : Defines the entry point for the console application.
//  http://www.glfw.org/docs/latest/quick.html

#include "stdafx.h"
#include <time.h>
#define _USE_MATH_DEFINES // for C++
#include <cmath>

#define SCREEN_WIDTH 1500
#define SCREEN_HEIGHT 1000

#define CENTER_X SCREEN_WIDTH/2
#define CENTER_Y SCREEN_HEIGHT/2
//------------IDS--------------
#define SUN_ID 0
#define MERCURY_ID 1
#define VENUS_ID 2
#define EARTH_ID 3
#define MARS_ID 4
#define JUPITER_ID 5
#define SATURN_ID 6
#define UrAnus_ID 7
#define PLANEPTUNE_ID 8
#define PLUTO_ID 9
//----------SIZES--------------
#define SUN_SIZE 70
#define MERCURY_SIZE 8
#define VENUS_SIZE 11
#define EARTH_SIZE 11
#define LUNA_SIZE 3
#define MARS_SIZE 8
#define JUPITER_SIZE 58
#define JUPITER_MOON_SIZE_MIN 2
#define JUPITER_MOON_SIZE_MAX 4
#define SATURN_SIZE 37
#define SATURN_MOON_SIZE_MIN 2
#define SATURN_MOON_SIZE_MAX 4
#define UrAnus_SIZE 25
#define PLANEPTUNE_SIZE 23
#define PLUTO_SIZE 3
//----------DISTS--------------
#define SUN_DIST 0
#define MERCURY_DIST 170
#define VENUS_DIST 210
#define EARTH_DIST 245
#define LUNA_DIST (EARTH_SIZE + 8)	//relatively to earth
#define MARS_DIST 290
#define AST_BELT_START 340
#define AST_BELT_STOP 355
#define JUPITER_DIST 460
#define JUPITER_MOON_DIST (JUPITER_SIZE + 8)	//relatively to saturn
#define SATURN_DIST 600
#define SATURN_MOON_DIST (SATURN_SIZE + 10)	//relatively to saturn
#define UrAnus_DIST 685
#define PLANEPTUNE_DIST 750
#define PLUTO_DIST 795
//---------OFFSETS------------
#define SUN_OFFSET 0
#define MERCURY_OFFSET 0
#define VENUS_OFFSET 0
#define EARTH_OFFSET 0
#define MARS_OFFSET 0
#define JUPITER_OFFSET 0
#define SATURN_OFFSET 0
#define UrAnus_OFFSET 15 
#define PLANEPTUNE_OFFSET -20
#define PLUTO_OFFSET 40
//--------VELOCITIES----------
#define MERCURY_SPEED 1.5
#define VENUS_SPEED 3
#define EARTH_SPEED 6.3
#define MARS_SPEED 9
#define JUPITER_SPEED 25
#define SATURN_SPEED 34
#define UrAnus_SPEED 50 
#define PLANEPTUNE_SPEED 60
#define PLUTO_SPEED 70
//-----------MISC------------
#define JUPITER_MOON_NUM 1
#define SATURN_MOON_NUM 10

GLdouble A, B, C, D;

int flag;

static void cursor_callback(GLFWwindow* window, double x, double y)
{
	
}

static void mouse_callback(GLFWwindow* window, int button, int action, int mods)
{
	if(button == GLFW_MOUSE_BUTTON_RIGHT)
	{
		if(action == GLFW_PRESS) glfwSetInputMode( window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
		if(action == GLFW_RELEASE) glfwSetInputMode( window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	}

	if(button == GLFW_MOUSE_BUTTON_LEFT)
	{
		
	}
}

static void resize_callback(GLFWwindow* window, int width, int height)
{
//	windowWidth = width;
//	windowHeight = height;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho( 0.0, (GLdouble)width, 0.0, (GLdouble)height, -1, 1);   

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
		
	A = width / 4.0;
	B = 0.0;
	C = D = height / 2.0;

	printf("Reshape occured\n");
}

static void keyboard_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);
}

void draw_stars(int num)
{
	int i;

	glPointSize(1.0);
	srand(time(NULL));
	for (i = 0; i < num; i++) {
		glBegin(GL_POINTS);
			glColor3f(1, 1, 1);
			glVertex2i(rand() % SCREEN_WIDTH, rand() % SCREEN_HEIGHT);
		glEnd();
	}
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	if (yoffset > 1) {
		flag++;
		draw_stars(1000);
	}
	else {
		flag--;
	}
}

static void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}



void draw_space_object(int id, int dist, int offset, int size, GLfloat angle)
{
	GLfloat phi = 0, theta = 0;
	int i = 0, moon_size = 5, moon_dist = 5, moon_speed = 10, start_pos;
	
	if (id == SUN_ID) {

	}
	if (id == EARTH_ID) {
		glPointSize(2.0);
		glBegin(GL_POLYGON);	//luna
		for (theta = 0; theta < 2 * M_PI; theta += 0.1) {
			glColor3f(0.8, 0.8, 0.9);
			glVertex2f(CENTER_X + offset + dist*cos(angle) - LUNA_DIST*sin(angle*EARTH_SPEED) + LUNA_SIZE*cos(theta), CENTER_Y + offset + dist*sin(angle) + LUNA_DIST*cos(angle*EARTH_SPEED) + LUNA_SIZE*sin(theta));
		}
		glEnd();
	}

	srand(time(NULL));
	if (id == JUPITER_ID) {
		for (i = 0; i < JUPITER_MOON_NUM; i++) {
			glPointSize(1.0);
			glBegin(GL_POLYGON);	//jupiter rings
			for (theta = 0; theta < 2 * M_PI; theta += 0.1) {
				glColor3f(1, 1, 1);
				moon_size = rand() % JUPITER_MOON_SIZE_MAX + JUPITER_MOON_SIZE_MIN;
				moon_dist = rand() % JUPITER_MOON_NUM * 10 + JUPITER_MOON_DIST + i * 3;
				moon_speed = 10;
				glVertex2f(CENTER_X + offset + dist*cos(angle) + (JUPITER_MOON_DIST + i * 3)*cos(angle*(moon_speed + i) + i) + moon_size*cos(theta), CENTER_Y + offset + dist*sin(angle) + (JUPITER_MOON_DIST + i * 3)*sin(angle*(moon_speed + i) + i) + moon_size*sin(theta));
			}
			glEnd();
		}
	}
	if (id == SATURN_ID) {
		glPointSize(2.0);
		glBegin(GL_POLYGON);	//saturn ring = 1.9*saturn_size
		for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
			glColor3f(0.8, 0.5, 0.2);
			glVertex2f(CENTER_X + offset + dist*cos(angle) + size*1.9*cos(phi), CENTER_Y + offset + dist*sin(angle) + size*1.9*sin(phi));
		}
		glEnd();

		glPointSize(2.0);
		glBegin(GL_POLYGON);	//space between ring and planet is just  black ring
		for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
			glColor3f(0, 0, 0);
			glVertex2f(CENTER_X + offset + dist*cos(angle) + size*1.1*cos(phi), CENTER_Y + offset + dist*sin(angle) + size*1.1*sin(phi));
		}
		glEnd();
		for (i = 0; i < SATURN_MOON_NUM; i++) {
			glPointSize(1.0);
			glBegin(GL_POLYGON);	//saturn rings
			for (theta = 0; theta < 2 * M_PI; theta += 0.1) {
				glColor3f(1, 0.4, 0.2);
				moon_size = rand() % SATURN_MOON_SIZE_MAX + SATURN_MOON_SIZE_MIN;
				moon_dist = rand() % SATURN_MOON_NUM * 10 + SATURN_MOON_DIST + i * 3;
				moon_speed = 15;
				glVertex2f(CENTER_X + offset + dist*cos(angle) + (SATURN_MOON_DIST + i * 3)*cos(angle*(moon_speed + i*0.8) + i) + moon_size*cos(theta), CENTER_Y + offset + dist*sin(angle) + (SATURN_MOON_DIST + i * 3)*sin(angle*(moon_speed + i*0.8) + i) + moon_size*sin(theta));
			}
			glEnd();
		}
	}
	glPointSize(2.0);
	glBegin(GL_POLYGON);
	for (phi = 0; phi < 2*M_PI; phi += 0.1) {
		//glColor3f(1, 1, 1);
		if (id == SUN_ID) glColor3f(253, 160, 0.2);
		if (id == MERCURY_ID) glColor3f(0.5, 0.5, 0.5);
		if (id == VENUS_ID) glColor3f(1, 0.7, 0);
		if (id == EARTH_ID) glColor3f(0, 0.4, 1);
		if (id == MARS_ID) glColor3f(0.5, 0, 0);
		if (id == JUPITER_ID) glColor3f(1, 0.6, 0.3);
		if (id == SATURN_ID) glColor3f(1, 0.5, 0.2);
		if (id == UrAnus_ID) glColor3f(0, 0.7, 0.7);
		if (id == PLANEPTUNE_ID) glColor3f(0, 0, 0.6);
		if (id == PLUTO_ID) glColor3f(1, 0.7, 0);
		glVertex2f(CENTER_X + offset + dist*cos(angle) + size*cos(phi), CENTER_Y + offset + dist*sin(angle) + size*sin(phi));
	}
	glEnd();
}

void draw_asteroid_belt()
{
	GLfloat phi = 0, theta = 0;

	glPointSize(2.0);
	glBegin(GL_POLYGON);
	for (phi = 0; phi < 2 * M_PI; phi++) {

	}
	glEnd();
}

void drawPoint(int x, int y)
{
	glPointSize(1.0);
	glBegin(GL_POINTS);
		glColor3f(0.1f, 0.4f, 0.4f);
		glVertex2i(x,y);
	glEnd();
}

typedef struct GLintPoint
{
	GLint	x;
	GLint	y;
} GLintPoint;

void serpinsky_drawing(void)
{
	GLintPoint T[3]= {{10,10},{210,310},{410,10}};

	int index = (int)floor((3.0 * rand()) / RAND_MAX);
	// 0, 1 or 2 equally likely
	// 0, 1 или 2 равновероятны
	GLintPoint point = T[index];
	// initial point
	// начальная точка
	drawPoint(point.x, point.y);
	// draw initial point
	// рисуем начальную точку
	for(int i = 0; i < 100000; i++)
		// draw 1000 dots
		// рисуем 1000 точек
	{
		index = (int)floor((3.0 * rand()) / RAND_MAX);
		point.x = (point.x + T[index].x) / 2;
		point.y = (point.y + T[index].y) / 2;
		drawPoint(point.x,point.y);
	}
}

void function_drawing(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
	
	glEnable(GL_LINE_STIPPLE);
	glColor3f(1.0f,0.3f,0.3f);
	//glBegin(GL_LINE_STRIP);
	//glBegin(GL_LINE_LOOP);
	glBegin(GL_POLYGON);
		for(GLdouble x = 0; x < 1.0; x += 0.05)
		{
			GLdouble func = exp(-x) * cos(2 * 3.14159265 * x);
			glVertex2d(A * x + B, C * func + D);
		}
	glEnd();
	glDisable(GL_LINE_STIPPLE);
	//glBegin(GL_LINE_STRIP);
	glBegin(GL_LINE_LOOP);
		for(GLdouble x = 1.0; x < 4.0; x += 0.05)
		{
			GLdouble func = exp(-x) * cos(2 * 3.14159265 * x);
			glVertex2d(A * x + B, C * func + D);
		}
	glEnd();


	glPointSize(4.0);
	glColor3b((GLbyte)20,(GLbyte)26,(GLbyte)250);
	glBegin(GL_POINTS);
		for(GLdouble x = 0; x < 4.0; x += 0.05)
		{
			GLdouble func = exp(-x) * cos(2 * 3.14159265 * x);
			glVertex2d(A * x + B, C * func + D);
		}
	glEnd();
}

void test_drawing_A(void)
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glColor3f(1.0,0.0,0.0);

	glPointSize(10);

	glBegin(GL_POINTS);
		glVertex2i(200,200);
		glVertex2i(300,300);
	glEnd();

	glBegin(GL_QUADS);
	{
		glVertex2i(200, 200);
		glVertex2i(100, 200);
		glVertex2i(100, 100);
		glVertex2i(200, 100);
	}
	glEnd();
}

void draw(void)
{
	//test_drawing();
	//test_drawing_A();	
	//serpinsky_drawing();
	//function_drawing();	
}

int main(int argc, _TCHAR* argv[])
{
	// initialise GLFW
    if(!glfwInit())
	{
		printf("glfwInit failed\n");
		return -1;
	}

	glfwSetErrorCallback(error_callback);

	GLFWwindow* window;
	
//	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 1);
//	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Test app", NULL, NULL);
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

	glfwSetKeyCallback(window, keyboard_callback);
	glfwSetFramebufferSizeCallback(window, resize_callback);
	glfwSetMouseButtonCallback(window, mouse_callback);
	glfwSetCursorPosCallback(window, cursor_callback);
	glfwSetScrollCallback(window, scroll_callback);
	
    resize_callback(window, SCREEN_WIDTH, SCREEN_HEIGHT);
	GLfloat phi = 0, center = 0;
	int i = 1;
	while (!glfwWindowShouldClose(window))
	{
		// ESC = exit
		// ^, v = change size
		// <--, --> = rotate left, right
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {
			i += 5;
		}
		if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {
			i -= 5;
		}
		if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) {
			phi += 0.3;
			draw_asteroid_belt();
			draw_stars(1000);
		}
		if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) {
			phi -= 0.3;
			draw_asteroid_belt();
			draw_stars(1000);
		}
		
		draw_space_object(SUN_ID, SUN_DIST, SUN_OFFSET, SUN_SIZE + i, phi + M_PI);
		draw_space_object(MERCURY_ID, MERCURY_DIST + i, MERCURY_OFFSET, MERCURY_SIZE, phi / MERCURY_SPEED + M_PI);
		draw_space_object(VENUS_ID, VENUS_DIST + i, VENUS_OFFSET, VENUS_SIZE, phi / VENUS_SPEED + M_PI);
		draw_space_object(EARTH_ID, EARTH_DIST + i, EARTH_OFFSET, EARTH_SIZE, phi / EARTH_SPEED + M_PI);
		draw_space_object(MARS_ID, MARS_DIST + i, MARS_OFFSET, MARS_SIZE, phi / MARS_SPEED + M_PI);
		
		draw_space_object(JUPITER_ID, JUPITER_DIST + i, JUPITER_OFFSET, JUPITER_SIZE, phi / JUPITER_SPEED);
		draw_space_object(SATURN_ID, SATURN_DIST + i, SATURN_OFFSET, SATURN_SIZE, phi / SATURN_SPEED);
		draw_space_object(UrAnus_ID, UrAnus_DIST + i, UrAnus_OFFSET, UrAnus_SIZE, phi / UrAnus_SPEED);
		draw_space_object(PLANEPTUNE_ID, PLANEPTUNE_DIST + i, PLANEPTUNE_OFFSET, PLANEPTUNE_SIZE, phi / PLANEPTUNE_SPEED - M_PI / 10);
		draw_space_object(PLUTO_ID, PLUTO_DIST + i, PLUTO_OFFSET, PLUTO_SIZE, phi / PLUTO_SPEED);

		glfwSwapBuffers(window);
		//glfwPollEvents();
		glfwWaitEvents();
	}

	glfwDestroyWindow(window);

	// clean up and exit
    glfwTerminate();

	return 0;
}

