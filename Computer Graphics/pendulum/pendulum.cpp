//pendulum.h
#pragma once
class pendulum
{
public:
	coord carrier;
	coord init_ticker;
	coord curr_ticker;

	float rope_len;
	float curr_angle, angle_degree, start_angle;
	float amplitude, start_amplitude;
	int mass;

	void set(coord carrier, coord ticker);
	void set_physics(bool k);

	void count_rope_len();
	void count_angle();

	void draw();
	void move(float t);

	bool get_physics();

	pendulum();
	~pendulum();

private:
	bool physics;
	int ticker_rad, ticker_border;
};





//pendulum.cpp
#include "stdafx.h"

using namespace std;

pendulum::pendulum() {
	physics = false;
}

pendulum::~pendulum()
{
}

void pendulum::set_physics(bool ka) {
	this->physics = ka;
}

bool pendulum::get_physics() {
	return this->physics;
}

void pendulum::set(coord carrier, coord ticker) {

	this->carrier = carrier;
	
	this->init_ticker = ticker;
	this->curr_ticker = ticker;

	this->start_angle = 0;

	this->count_angle();

	this->amplitude = this->curr_angle;
	this->start_amplitude = this->amplitude;

	this->physics = false;

	this->ticker_rad = 40;
	this->ticker_border = 5;
}

void pendulum::count_angle() {

	this->count_rope_len();

	float TM = abs(this->carrier.x - this->curr_ticker.x);
	float CT = this->rope_len;

	float angle = asin(TM / CT);
	
	this->curr_angle = angle;
	this->angle_degree = angle * 180 / M_PI;
}

void pendulum::count_rope_len() {
	coord carr = this->carrier;
	coord tick = this->curr_ticker;

	float CM = abs(carr.y - tick.y);
	float TM = abs(carr.x - tick.x);
	float CT = sqrt(TM*TM + CM*CM);

	this->rope_len = CT;
}

void draw_pendulum_bracket(int x, int y, int size) {

	glLineWidth(3.0);
	glColor3f(0, 0, 0);

	glBegin(GL_LINES);
	glVertex2f(x, y);
	glVertex2f(x, y - size);
	glEnd();

	//outer
	GLfloat phi = 0, theta = 0;
	glBegin(GL_POLYGON);
	glColor3f(0, 0, 0);
	for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
		glVertex2f(x + 8 * cos(phi), (y - size) + 8 * sin(phi));
	}
	glEnd();

	//inner
	glBegin(GL_POLYGON);
	glColor3f(0.6, 0.6, 0.6);
	for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
		glVertex2f(x + 5 * cos(phi), (y - size) + 5 * sin(phi));
	}
	glEnd();

}

void pendulum::draw() {

	coord carrier = this->carrier;
	coord ticker = this->curr_ticker;

	
	//carrier
	draw_pendulum_bracket(carrier.x, carrier.y + 20, 15);

	//ticker outer
	GLfloat phi = 0, theta = 0;
	glBegin(GL_POLYGON);
	glColor3f(0, 0, 0);
	for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
		glVertex2f(ticker.x + this->ticker_rad * cos(phi), ticker.y + this->ticker_rad * sin(phi));
	}
	glEnd();

	//ticker inner
	glBegin(GL_POLYGON);
	glColor3f(1, 1, 1);
	for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
		glVertex2f(ticker.x + (this->ticker_rad - this->ticker_border) * cos(phi), ticker.y + (this->ticker_rad - this->ticker_border) * sin(phi));
	}
	glEnd();

	//ticker center
	glPointSize(7.0);
	glColor3f(0, 0, 0);
	glBegin(GL_POINTS);
		glVertex2f(ticker.x, ticker.y);
	glEnd();

	//rope
	glLineWidth(3.0);
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
		glVertex2f(carrier.x, carrier.y);
		glVertex2f(ticker.x, ticker.y);
	glEnd();
}

void pendulum::move(float t) {

	coord carr = this->carrier;
	coord tick = this->curr_ticker;

	if (tick.y - this->ticker_rad > FLOOR_HEIGHT && tick.y + this->ticker_rad < ROOF_HEIGHT) {
		this->count_angle();

		float omega = sqrt(G_CONST / this->rope_len);
		//float attenuation = pow(M_E, -log(this->amplitude / (this->amplitude * 0.5)));
		float attenuation_k = 0.99999;
		float phi = this->amplitude * cos(omega * t + this->start_angle);
		this->amplitude *= attenuation_k;

		this->curr_angle = phi;

		float TM = this->rope_len * sin(phi);
		float CM = this->rope_len * sin(90*M_PI/180 - phi);	//1,57 rad

		coord T;
		T.x = carr.x - TM;
		T.y = carr.y - CM;

		this->curr_ticker = T;
	}
}

