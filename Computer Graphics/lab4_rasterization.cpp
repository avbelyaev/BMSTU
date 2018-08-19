//закраска упорядоченным списком ребер и сглаживание усредненной постфильтрацией 3х3
//управление: мышь, q-толщина линий, d-рисовать, a-замкнуть контур, f-заливка, s-постфильтрация, t-распечатать упорядоченный список ребер

// glfw_30.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <time.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <GL/glew.h>
#include <GL\glut.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <algorithm>
#include <cstdlib>

int SCREEN_WIDTH = 800, SCREEN_HEIGHT = 600;

#define	SCREEN_POS_X		100
#define	SCREEN_POS_Y		100

#define DATA_FORMAT			3
#define BUFFER_SIZE			SCREEN_WIDTH*SCREEN_HEIGHT*DATA_FORMAT
#define PXL_INDX(X,Y,COMP) (((Y)*SCREEN_WIDTH+(X))*DATA_FORMAT+(COMP))

typedef struct point{
	int x, y;
	bool extremum;
	point* incident_next;
	point* incident_prev;
} point;

typedef struct edge{
	point a, b;
	edge* incident_e_next;
	edge* incident_e_prev;
}edge;

bool draw_flag, fill_flag, anti_aliacing_flag, update_flag, thick_flag, clear_flag;
int x_max, x_min, y_max, y_min;

GLfloat *contour;
GLfloat *intrails;
GLfloat *blur;

std::vector<point> vertices;
std::vector<edge> edges;
std::vector<std::vector<int>> edge_table;



bool extremum(int x, int y)
{
	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i].x == x && vertices[i].y == y) {
			return vertices[i].extremum;
		}
	}

	return false;
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

void clear_vector(std::vector<point>& vector)
{
	vector.clear();
	vector.shrink_to_fit();
}

void clear_buffers()
{
	delete[] contour;
	delete[] intrails;
	delete[] blur;
}

void copy_buffer(GLfloat *buff_source, GLfloat *buff_dest)
{
	for (int i = 0; i < BUFFER_SIZE; i++) {
		buff_dest[i] = buff_source[i];
	}
}

void create_buffers()
{
	contour = new GLfloat[BUFFER_SIZE];
	intrails = new GLfloat[BUFFER_SIZE];
	blur = new GLfloat[BUFFER_SIZE];
}

void figure_pos()
{
	int max_x = 0, min_x = SCREEN_WIDTH;
	int max_y = 0, min_y = SCREEN_HEIGHT;

	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i].x > max_x) max_x = vertices[i].x;
		if (vertices[i].x < min_x) min_x = vertices[i].x;
		if (vertices[i].y > max_y) max_y = vertices[i].y;
		if (vertices[i].y < min_y) min_y = vertices[i].y;
	}

	x_max = max_x;
	x_min = min_x;
	y_max = max_y;
	y_min = min_y;
}


//===========================DRAWING===========================
void Bresenham(int x1, int y1, int const x2, int const y2)
{
	int delta_x(x2 - x1);
	// if x1 == x2, then it does not matter what we set here
	signed char const ix((delta_x > 0) - (delta_x < 0));
	delta_x = std::abs(delta_x) << 1;

	int delta_y(y2 - y1);
	// if y1 == y2, then it does not matter what we set here
	signed char const iy((delta_y > 0) - (delta_y < 0));
	delta_y = std::abs(delta_y) << 1;



	//plot(x1, y1);
	for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1, y1, i)] = 1;
	//========================================
	if (thick_flag) {
		for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 - 1, y1 + 1, i)] = 1;
		for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1, y1 + 1, i)] = 1;
		for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 + 1, y1 + 1, i)] = 1;
		for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 - 1, y1, i)] = 1;
		for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 + 1, y1, i)] = 1;
		for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 - 1, y1 - 1, i)] = 1;
		for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1, y1 - 1, i)] = 1;
		for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 + 1, y1 - 1, i)] = 1;
	}
	//========================================



	if (delta_x >= delta_y)
	{
		// error may go below zero
		int error(delta_y - (delta_x >> 1));

		while (x1 != x2)
		{
			if ((error >= 0) && (error || (ix > 0)))
			{
				error -= delta_x;
				y1 += iy;
			}
			// else do nothing

			error += delta_y;
			x1 += ix;


			//plot(x1, y1);
			for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1, y1, i)] = 1;
			//========================================
			if (thick_flag) {
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 - 1, y1 + 1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1, y1 + 1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 + 1, y1 + 1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 - 1, y1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 + 1, y1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 - 1, y1 - 1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1, y1 - 1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 + 1, y1 - 1, i)] = 1;
			}
			//========================================
		}
	}
	else
	{
		// error may go below zero
		int error(delta_x - (delta_y >> 1));

		while (y1 != y2)
		{
			if ((error >= 0) && (error || (iy > 0)))
			{
				error -= delta_y;
				x1 += ix;
			}
			// else do nothing

			error += delta_x;
			y1 += iy;



			//plot(x1, y1);
			for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1, y1, i)] = 1;
			//========================================
			if (thick_flag) {
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 - 1, y1 + 1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1, y1 + 1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 + 1, y1 + 1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 - 1, y1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 + 1, y1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 - 1, y1 - 1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1, y1 - 1, i)] = 1;
				for (int i = 0; i < DATA_FORMAT; i++) contour[PXL_INDX(x1 + 1, y1 - 1, i)] = 1;
			}
			//========================================
		}
	}
}

void set_contour()//D
{
	for (int i = 0; i < vertices.size() - 1; i++) 
		Bresenham(vertices[i].x, vertices[i].y, vertices[i + 1].x, vertices[i + 1].y); 
}

void draw_contour()
{
	glDrawPixels(SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGB, GL_FLOAT, (void*)contour); 
}


//===========================FILLING===========================
void intersect_counter(int x1, int y1, int const x2, int const y2)
{
	int prev, curr;
	prev = y1;
	curr = y1;

	int delta_x(x2 - x1);
	signed char const ix((delta_x > 0) - (delta_x < 0));
	delta_x = std::abs(delta_x) << 1;

	int delta_y(y2 - y1);
	signed char const iy((delta_y > 0) - (delta_y < 0));
	delta_y = std::abs(delta_y) << 1;

	if (extremum(x1, y1)) {
		//edge_table[y1 - y_min].push_back(x1);
		//edge_table[y1 - y_min].push_back(x1);
	}
	else {
		//edge_table[y1 - y_min].push_back(x1);
	}

	if (delta_x >= delta_y)
	{
		int error(delta_y - (delta_x >> 1));
		while (x1 != x2)
		{
			if ((error >= 0) && (error || (ix > 0)))
			{
				error -= delta_x;
				y1 += iy;
				curr = y1;
			}
			error += delta_y;
			x1 += ix;

			if (curr != prev) {
				if (extremum(x1, y1)) {
					edge_table[y1 - y_min].push_back(x1);
					edge_table[y1 - y_min].push_back(x1);
				}
				else {
					edge_table[y1 - y_min].push_back(x1);
				}
			}
			prev = curr;
		}
	}
	else
	{
		int error(delta_x - (delta_y >> 1));
		while (y1 != y2)
		{
			if ((error >= 0) && (error || (iy > 0)))
			{
				error -= delta_y;
				x1 += ix;
			}
			error += delta_x;
			y1 += iy;

			if (extremum(x1, y1)) {
				edge_table[y1 - y_min].push_back(x1);
				edge_table[y1 - y_min].push_back(x1);
			}
			else {
				edge_table[y1 - y_min].push_back(x1);
			}
		}
	}
}

void set_vertices()
{
	int i;

	//for 0th vertex
	vertices[0].extremum = false;
	vertices[0].incident_next = &vertices[1];
	vertices[0].incident_prev = &vertices[vertices.size() - 1];

	//for other vertices
	for (i = 1; i < vertices.size(); i++) {
		vertices[i].incident_next = &vertices[(i + 1) % vertices.size()];
		vertices[i].incident_prev = &vertices[i - 1];
		vertices[i].extremum = false;
	}

	//find extermums
	for (i = 0; i < vertices.size(); i++) {
		if (
			((*vertices[i].incident_next).y > vertices[i].y && (*vertices[i].incident_prev).y > vertices[i].y) ||
			((*vertices[i].incident_next).y < vertices[i].y && (*vertices[i].incident_prev).y < vertices[i].y)
			) {
			vertices[i].extremum = true;
		}

	}
}

void init_edges()
{
	edge temp_edge;

	for (int i = 0; i < vertices.size(); i++) {
		temp_edge.a = vertices[i];
		temp_edge.b = vertices[(i + 1) % vertices.size()];

		edges.push_back(temp_edge);
	}
}

void set_edges()//F
{
	set_vertices();

	init_edges();

	int i;
	
	edge_table.resize(y_max);

	for (i = 0; i < edge_table.size(); i++) 
		edge_table[i].resize(0);
	
	for (i = 0; i < vertices.size(); i++)
		intersect_counter(vertices[i].x, vertices[i].y, vertices[(i + 1) % vertices.size()].x, vertices[(i + 1) % vertices.size()].y);

	for (i = 0; i < edge_table.size(); i++)
		std::sort(edge_table[i].begin(), edge_table[i].end());
}

void set_intrails_pixels()
{
	int i, j, z, a;

	for (i = 0; i < edge_table.size(); i++) {
		if (edge_table[i].size() % 2 == 0 && 0 != edge_table[i].size()) {
			for (j = 0; j < edge_table[i].size() - 1; j += 2) {
				int start_filling = edge_table[i][j];
				int stop_filling = edge_table[i][j + 1];
				for (z = start_filling; z < stop_filling; z++) {
					for (a = 0; a < 3; a++) {
						intrails[PXL_INDX(z, i + y_min, a)] = 1;
					}
				}
			}
		}
	}
}

void draw_intrails()
{
	glDrawPixels(SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGB, GL_FLOAT, (void*)intrails);
}


//===========================BLUR===========================
void postfiltration()//S
{
	int i, j, k, x, y, num;
	GLfloat r, g, b;

	//copy_buffer(contour, blur);

	for (y = 0; y < SCREEN_HEIGHT; y++) {
		for (x = 0; x < SCREEN_WIDTH; x++) {
			r = 0;
			g = 0;
			b = 0;
			num = 0;

			if (x + 1 < SCREEN_WIDTH && y + 1 < SCREEN_HEIGHT) {
				num++;
				r += intrails[PXL_INDX(x + 1, y + 1, 0)];
				g += intrails[PXL_INDX(x + 1, y + 1, 1)];
				b += intrails[PXL_INDX(x + 1, y + 1, 2)];
			}
			if (y + 1 < SCREEN_HEIGHT) {
				num++;
				r += intrails[PXL_INDX(x, y + 1, 0)];
				g += intrails[PXL_INDX(x, y + 1, 1)];
				b += intrails[PXL_INDX(x, y + 1, 2)];
			}
			if (x - 1 > 0 && y + 1 < SCREEN_HEIGHT) {
				num++;
				r += intrails[PXL_INDX(x - 1, y + 1, 0)];
				g += intrails[PXL_INDX(x - 1, y + 1, 1)];
				b += intrails[PXL_INDX(x - 1, y + 1, 2)];
			}

			if (x + 1 < SCREEN_WIDTH) {
				num++;
				r += intrails[PXL_INDX(x + 1, y, 0)];
				g += intrails[PXL_INDX(x + 1, y, 1)];
				b += intrails[PXL_INDX(x + 1, y, 2)];
			}
			if (1) {
				num++;
				r += intrails[PXL_INDX(x, y, 0)];
				g += intrails[PXL_INDX(x, y, 1)];
				b += intrails[PXL_INDX(x, y, 2)];
			}
			if (x - 1 > 0) {
				num++;
				r += intrails[PXL_INDX(x - 1, y, 0)];
				g += intrails[PXL_INDX(x - 1, y, 1)];
				b += intrails[PXL_INDX(x - 1, y, 2)];
			}

			if (x + 1 < SCREEN_WIDTH && y - 1 > 0) {
				num++;
				r += intrails[PXL_INDX(x + 1, y - 1, 0)];
				g += intrails[PXL_INDX(x + 1, y - 1, 1)];
				b += intrails[PXL_INDX(x + 1, y - 1, 2)];
			}
			if (y - 1 > 0) {
				num++;
				r += intrails[PXL_INDX(x, y - 1, 0)];
				g += intrails[PXL_INDX(x, y - 1, 1)];
				b += intrails[PXL_INDX(x, y - 1, 2)];
			}
			if (x - 1 > 0 && y - 1 > 0) {
				num++;
				r += intrails[PXL_INDX(x - 1, y - 1, 0)];
				g += intrails[PXL_INDX(x - 1, y - 1, 1)];
				b += intrails[PXL_INDX(x - 1, y - 1, 2)];
			}

			blur[PXL_INDX(x, y, 0)] = r / (GLfloat)num;
			blur[PXL_INDX(x, y, 1)] = g / (GLfloat)num;
			blur[PXL_INDX(x, y, 2)] = b / (GLfloat)num;
			
		}

	}


//radiance goes here. its not necessary for lab
	/*
	for (y = 0; y < SCREEN_HEIGHT; y++) {
		for (x = 0; x < SCREEN_WIDTH; x++) {
			if (contour[PXL_INDX(x, y, 0)] == 1 && 
				contour[PXL_INDX(x, y, 1)] == 1 && 
				contour[PXL_INDX(x, y, 2)] == 1
				) {
				for (i = 0; i < 3; i++) {
					blur[PXL_INDX(x, y, i)] = 1;
				}
			}
		}
	}*/
}

void draw_blur()//S
{
	glDrawPixels(SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGB, GL_FLOAT, (void*)blur);
}


//===========================MAIN===========================
void renderScene(void) 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (draw_flag) {
		if (0 != vertices.size()) {
			set_contour();
			draw_contour();
		}
	}

	if (fill_flag) {
		set_edges();
		set_intrails_pixels();
		draw_intrails();
	}

	if (anti_aliacing_flag) {
		postfiltration();
		draw_blur();
	}
	if (clear_flag) {
		clear_buffers();
		clear_vector(vertices);
	}

	glutSwapBuffers();
}

void mouseButton(int button, int state, int x, int y)
{
	if (GLUT_LEFT_BUTTON == button) {
		point a;
		if (GLUT_DOWN == state) {
			a.x = x;
			a.y = SCREEN_HEIGHT - y;
			printf("%d %d\n", a.x, a.y);
			vertices.push_back(a);
			figure_pos();
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
		draw_flag = !draw_flag;
		figure_pos();
		printf("x_min=%d x_max=%d\ny_min=%d y_max=%d\n", x_min, x_max, y_min, y_max);
		draw_flag ? printf("drawing:ON\n") : printf("drawing:OFF\n");
		glutPostRedisplay();
	}
	if ('s' == key || 'S' == key) {
		anti_aliacing_flag = !anti_aliacing_flag;
		anti_aliacing_flag ? printf("PS4 / Xbox ONE mode\n") : printf("PC mode\n");
		glutPostRedisplay();
	}
	if ('f' == key || 'F' == key) {
		fill_flag = !fill_flag;
		figure_pos();
		fill_flag ? printf("filling:ON\n") : printf("filling:OFF\n");
		glutPostRedisplay();
	}
	/*if ('z' == key || 'Z' == key) {
		if (0 < vertices.size()) {
			vertices.pop_back();
			figure_pos();
			glutPostRedisplay();
		}
	}*/
	if ('a' == key || 'A' == key) {
		Bresenham(vertices[vertices.size() - 1].x, vertices[vertices.size() - 1].y, vertices[0].x, vertices[0].y);
		glutPostRedisplay();
	}
	if ('q' == key || 'Q' == key) {
		thick_flag = !thick_flag;
		thick_flag ? printf("thick lines:ON\n") : printf("thick lines:OFF\n");
		glutPostRedisplay();
	}
	if ('w' == key || 'W' == key) {
		clear_flag = !clear_flag;
		printf("cleared\n");
		glutPostRedisplay();
	}
	if ('t' == key || 'T' == key) {
		printf("table_size=%d\n", edge_table.size());
		for (int i = 0; i < edge_table.size(); i++) {
			printf("[%d](%d):", i, edge_table[i].size());
			for (int j = 0; j < edge_table[i].size(); j++) {
				printf("%d ", edge_table[i][j]);
			}
			printf(" -> ");
			if (edge_table[i].size() % 2 == 0 && 0 != edge_table[i].size()) {
				for (int j = 0; j < edge_table[i].size() - 1; j += 2) {
					int start = edge_table[i][j];
					int stop = edge_table[i][j + 1];
					printf("%d..%d ", start, stop);
				}
				
			}
			else {
				printf("MISTAKE");
			}
			printf("\n");
		}

		printf("vertices:\n");
		for (int i = 0; i < vertices.size(); i++) {
			printf("(%d)", i);
			if (vertices[i].extremum == true) {
				printf("extremum\n");
			}
			else {
				printf("normal\n");
			}
		}
		glutPostRedisplay();
	}
}

void initial_state()
{
	glClearColor(0, 0, 0, 1);

	draw_flag = false;
	fill_flag = false;
	anti_aliacing_flag = false;
	update_flag = false;
	thick_flag = false;

	thick_flag ? printf("thick lines:ON\n") : printf("thick lines:OFF\n");

	x_max = 0;
	x_min = SCREEN_WIDTH;
	y_max = 0;
	y_min = SCREEN_HEIGHT;

	create_buffers();
	clear_buffers();
}

int main(int argc, char **argv) 
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(SCREEN_POS_X, SCREEN_POS_Y);
	glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);
	glutCreateWindow("Anthony's Project: Lab 4");

	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutIdleFunc(renderScene);

	glutMouseFunc(mouseButton);
	glutKeyboardFunc(keyboardButton);

	initial_state();

	glutMainLoop();

	return 0;
}

