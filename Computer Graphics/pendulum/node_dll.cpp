//===================================================================================================
//===================================================================================================
//===================================================================================================
//node.h
#pragma once
class node
{
public:
	int id;
	char type;
	int listpos;
	int mass;
	float velo, accel;
	bool halt_move, halt_pend, explicit_accel;

	node *next_node;
	node *prev_node;

	node *cp_left_node, *cp_center_node, *cp_right_node;
	int cp_left_id, cp_center_id, cp_right_id;

	coord init;
	coord curr;

	float pend_angle, pend_rope_len, pend_amplitude, pend_start_amplitude;
	coord ticker, carrier;

	void count_pend_param();

	void set_list_param(char type, int id, int listpos, float pos_x, float pos_y, int mass);
	void set_cp_ids(int left, int center, int right);
	void set_cp_ids(int center);
	void set_list_size(int list_size, int cargos_amount, int blocks_amount);

	void draw_body();
	void draw_ropes();

	void phys_pend(float t);
	void phys_move(float t);
	void deliver_motion_to_list();
	void deliver_motion_by_id(int id);

	node *get_by_id(int id);
	node *get_next_by_type(char type);
	node *get_prev_by_type(char type);
	coord get_cp_left_pos();
	coord get_cp_center_pos();
	coord get_cp_right_pos();

	node(int id);
	~node();

private:
	int list_size, cargos_amount, blocks_amount;
	int size_x, size_y, size_border, size_rope, size_bracket;
	coord cp_left_pos, cp_center_pos, cp_right_pos;

};




//===================================================================================================
//===================================================================================================
//===================================================================================================
//node.cpp
#include "stdafx.h"

using namespace std;

node::node(int x)
{
	id = x; 
	next_node = NULL;
	prev_node = NULL;
	halt_move = false;
	halt_pend = false;
}

node::~node()
{
}

void node::set_list_param(char type, int id, int listpos, float pos_x, float pos_y, int mass) {
	this->init.x = pos_x;
	this->init.y = pos_y;
	this->curr.x = pos_x;
	this->curr.y = pos_y;

	this->type = type;

	this->id = id;

	this->listpos = listpos;

	this->mass = mass;

	this->size_x = 100;
	this->size_y = 100;
	this->size_border = 5;
	this->size_rope = 3;
	this->size_bracket = 15;

	this->velo = 0;
	this->accel = (-1)*G_CONST;

	this->halt_move = false;
	this->halt_pend = false;
	this->explicit_accel = false;
}

void node::set_cp_ids(int left, int center, int right) {
	this->cp_left_id = left;
	this->cp_center_id = center;
	this->cp_right_id = right;
}

void node::set_cp_ids(int center) {
	this->cp_center_id = center;
}

void node::set_list_size(int size, int cargos, int blocks) {
	//for every block to know are there another blocks
	this->list_size = size;
	this->cargos_amount = cargos;
	this->blocks_amount = blocks;
}

node* node::get_by_id(int id) {

	node *tmp = this;
	while (tmp != NULL) {
		if (tmp->prev_node == NULL) 
			break;
		tmp = tmp->prev_node;
	}

	while (tmp != NULL) {
		if (tmp->id == id) {
			return tmp;
		}
		tmp = tmp->next_node;
	}
	return NULL;
}

node* node::get_next_by_type(char type) {

	node *tmp = this;
	tmp = tmp->next_node;

	while (tmp != NULL) {
		if (tmp->type == type) {
			return tmp;
		}
		tmp = tmp->next_node;
	}
	return NULL;
}

node* node::get_prev_by_type(char type) {

	node *tmp = this;
	tmp = tmp->prev_node;

	while (tmp != NULL) {
		if (tmp->type == type) {
			return tmp;
		}
		tmp = tmp->prev_node;
	}
	return NULL;
}

coord node::get_cp_left_pos() {
	coord retval;
	retval.x = this->curr.x - this->size_x / 2;
	retval.y = this->curr.y;

	this->cp_left_pos.x = retval.x;
	this->cp_left_pos.y = retval.y;

	return retval;
}

coord node::get_cp_center_pos() {
	
	if ('b' == this->type) {
		cp_center_pos = curr;
		return this->curr;
	}

	coord retval;

	retval.x = this->curr.x;
	retval.y = this->curr.y + this->size_y / 2;

	this->cp_center_pos.x = retval.x;
	this->cp_center_pos.y = retval.y;

	return retval;
}

coord node::get_cp_right_pos() {
	coord retval;
	retval.x = this->curr.x + this->size_x / 2;
	retval.y = this->curr.y;

	this->cp_right_pos.x = retval.x;
	this->cp_right_pos.y = retval.y;

	return retval;
}

void draw_line(float from_x, float from_y, float to_x, float to_y) {
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
		glVertex2f(from_x, from_y);
		glVertex2f(to_x, to_y);
	glEnd();
}

void draw_line(coord from, coord to) {
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
		glVertex2f(from.x, from.y);
		glVertex2f(to.x, to.y);
	glEnd();
}

void draw_bracket(int x, int y, int size) {

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

	glBegin(GL_POLYGON);
	glColor3f(0.6, 0.6, 0.6);
	glVertex2f(x - 8, y - size);
	glVertex2f(x - 8, y - 7);
	glVertex2f(x - 2, y - 7);
	glVertex2f(x - 2, y - size);
	glEnd();
}

void node::draw_body() {
	float center_x = this->curr.x;
	float center_y = this->curr.y;

	int size_x = this->size_x;
	int size_y = this->size_y;
	int rad = size_x / 2;

	int border = this->size_border;

	if ('c' == this->type) {

		coord lb, lt, rt, rb;

		lb.x = center_x - size_x / 2;
		lb.y = center_y - size_y / 2;

		lt.x = center_x - size_x / 2;
		lt.y = center_y + size_y / 2;

		rt.x = center_x + size_x / 2;
		rt.y = center_y + size_y / 2;

		rb.x = center_x + size_x / 2;
		rb.y = center_y - size_y / 2;

		//outer black square
		glBegin(GL_POLYGON);
			glColor3f(0, 0, 0);
			glVertex2f(lb.x, lb.y);
			glColor3f(0.4, 0.4, 0.4);
			glVertex2f(lt.x, lt.y);
			glColor3f(0, 0, 0);
			glVertex2f(rt.x, rt.y);
			glVertex2f(rb.x, rb.y);
		glEnd();

		//inner white square
		glBegin(GL_POLYGON);
			glColor3f(0.7, 0.6, 0.6);
			glVertex2f(lb.x + border, lb.y + border);
			glColor3f(1, 1, 1);
			glVertex2f(lt.x + border, lt.y - border);
			glColor3f(0.8, 0.9, 0.8);
			glVertex2f(rt.x - border, rt.y - border);
			glColor3f(0.4, 0.4, 0.6);
			glVertex2f(rb.x - border, rb.y + border);
		glEnd();
	}

	if ('b' == this->type) {

		GLfloat phi = 0, theta = 0;

		//outer black circle
		glBegin(GL_POLYGON);
		glColor3f(0, 0, 0);
		for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
			glVertex2f(center_x + rad*cos(phi), center_y + rad*sin(phi));
		}
		glEnd();

		//inner white circle
		glBegin(GL_POLYGON);
		
		for (phi = 0; phi < 2 * M_PI; phi += 0.1) {
			glColor3f(0.7, 0.7, 0.8);
			if (cos(phi) < 0 && sin(phi) > 0) {
				glColor3f(0.75, 0.7, 0.75);
			}
			glVertex2f(center_x + (rad - border)*cos(phi), center_y + (rad - border)*sin(phi));
		}
		glEnd();

	}
}

void node::draw_ropes() {

	node *next = this->next_node;
	node *prev = this->prev_node;

	if (NULL == next) {
		return;
	}

	glLineWidth((GLfloat)this->size_rope);
	glColor3f(0, 0, 0);

	if ('b' == this->type) {

		//suspension to the ceiling
		if (0 == this->cp_left_id) {
			coord this_block_cp_left = this->get_cp_left_pos();
			draw_bracket(this_block_cp_left.x, ROOF_HEIGHT, this->size_bracket);
			draw_line(this_block_cp_left.x, this_block_cp_left.y, this_block_cp_left.x, ROOF_HEIGHT - this->size_bracket - 2);
		}
		if (0 == this->cp_center_id) {
			coord this_block_cp_center = this->get_cp_center_pos();
			draw_bracket(this_block_cp_center.x, ROOF_HEIGHT, this->size_bracket);
			draw_line(this_block_cp_center.x, this_block_cp_center.y, this_block_cp_center.x, ROOF_HEIGHT - this->size_bracket -2);
		}
		if (0 == this->cp_right_id) {
			coord this_block_cp_right = this->get_cp_right_pos();
			draw_bracket(this_block_cp_right.x, ROOF_HEIGHT, this->size_bracket);
			draw_line(this_block_cp_right.x, this_block_cp_right.y, this_block_cp_right.x, ROOF_HEIGHT - this->size_bracket -2);
		}

		//mark central point
		coord this_block_cp_center = this->get_cp_center_pos();
		glPointSize(7.0);
		glColor3f(0, 0, 0);
		glBegin(GL_POINTS);
			glVertex2f(this_block_cp_center.x, this_block_cp_center.y);
		glEnd();
		
	}



	if ('c' == this->type && 'b' == next->type) {
		if (NULL != this->prev_node && 'b' == this->prev_node->type && 'b' == this->next_node->type) return;
		coord this_cargo_cp = this->get_cp_center_pos();
		coord next_block_cp_left = next->get_cp_left_pos();

		draw_line(this_cargo_cp, next_block_cp_left);
	}

	if ('b' == this->type && 'b' == next->type) {
		coord this_block_cp_right = this->get_cp_right_pos();
		coord next_block_cp_left = next->get_cp_left_pos();

		draw_line(this_block_cp_right, next_block_cp_left);
	}

	
	if ('b' == this->type && 'c' == next->type) {
		//cargo hanging on a block
		if (0 != this->cp_center_id && -1 != this->cp_center_id) {
			coord this_block_cp_center = this->get_cp_center_pos();
			coord next_cargo_cp_center = this->next_node->get_cp_center_pos();

			draw_line(this_block_cp_center, next_cargo_cp_center);
			
		}
		//crago is twined
		if (0 == this->cp_center_id && 0 != this->cp_right_id) {
			coord this_block_cp_right = this->get_cp_right_pos();
			coord next_cargo_cp = next->get_cp_center_pos();

			draw_line(this_block_cp_right, next_cargo_cp);
			
		}
		if (NULL != next->next_node && 'b' == next->next_node->type) {
			coord this_block_cp_right = this->get_cp_right_pos();
			coord nearest_block_cp_left = next->next_node->get_cp_left_pos();

			draw_line(this_block_cp_right, nearest_block_cp_left);
		}
	}

	
}



void node::phys_move(float t) {

	if (0 != this->mass && false == this->halt_move) {
		float a, v, y0;
		float all_mass;
		float bot_limit = (float)FLOOR_HEIGHT;
		float top_limit = (float)ROOF_HEIGHT;
		node *next_cargo = this->get_next_by_type('c');
		node *prev_cargo = this->get_prev_by_type('c');
		node *next_block = this->get_next_by_type('b');
		node *prev_block = this->get_prev_by_type('b');

		if ('c' == this->type) {
			if (NULL != this->next_node) 
				top_limit = this->next_node->curr.y - this->next_node->size_y / 2;
			if (NULL != this->prev_node)
				top_limit = this->prev_node->curr.y - this->prev_node->size_y / 2;
		}

		if (NULL == prev_cargo && ((NULL != next_cargo) || (NULL != prev_block))) {

			all_mass = next_cargo->mass;
			node *curr = this;
			for (int i = 0; i < this->cargos_amount - this->listpos; i++) {
				node *tmp = curr->next_node;
				if ('c' == tmp->type) {
					all_mass += tmp->mass;
					i++;
				}
			}
			
			a = G_CONST *(all_mass - this->mass) / (all_mass + this->mass);
			v = this->velo;
			y0 = this->curr.y;
			
			this->accel = a;
			this->velo = v;
			this->curr.y = y0;
			
			this->deliver_motion_to_list();
			
		}
		count_coord:
		if (bot_limit <= this->curr.y - this->size_y / 2 && top_limit >= this->curr.y + this->size_y / 2) {

			this->curr.x = this->curr.x;

			this->curr.y = this->curr.y + this->velo*t + (this->accel*t*t) / 2;

		}
		else {
			this->halt_move = true;
		}
		

	}

}

void node::deliver_motion_to_list() {

	float accel = this->accel;
	node *curr = this->next_node;

	int i = 0;
	while (NULL != curr) {
		
		if ('c' == curr->type) {
			int his_blocks_id = curr->cp_center_id;
			node *his_block = this->get_by_id(his_blocks_id);

			if (NULL != his_block) {

				if (his_block->cp_center_id == curr->id) {
					//this cargo is hanging on block
					//curr->accel = -accel / 2;
					curr->deliver_motion_by_id(curr->cp_center_id);
				} else {
					//its twined
					//its also the last cargo
					curr->accel = -accel;
				}

			} else {
				printf("[node::deliver_motion_to_list]: couldnt find block[id:%d] for cargo[id:%d]\n", his_blocks_id, curr->id);
			}
		}

		if ('b' == curr->type) {
			if (0 != curr->cp_center_id && -1 != curr->cp_center_id) {
				curr->accel = -accel / 2;
				curr->deliver_motion_by_id(curr->cp_center_id);
			}
			i++;
		}
		curr = curr->next_node;
	}

}

void node::deliver_motion_by_id(int id) {

	node *tmp = this->get_by_id(id);
	tmp->accel = this->accel;
	tmp->explicit_accel = true;

}

void node::count_pend_param() {

	//find pos of carrier (it's cargo cp)
	if (NULL != this->prev_node) {
		if (this->prev_node->cp_center_id == this->id) {
			this->carrier = this->prev_node->get_cp_center_pos();
		}
		if (this->prev_node->cp_right_id == this->id) {
			this->carrier = this->prev_node->get_cp_right_pos();
		}
	}
	if (NULL == this->prev_node && NULL != this->next_node) {
		if (this->next_node->cp_center_id == this->id) {
			this->carrier = this->next_node->get_cp_center_pos();
		}
		if (this->next_node->cp_left_id == this->id) {
			this->carrier = this->next_node->get_cp_left_pos();
		}
	}
	this->ticker = curr;

	float CM = abs(this->carrier.y - this->ticker.y);
	float TM = abs(this->carrier.x - this->ticker.x);
	float CT = sqrt(TM*TM + CM*CM);

	this->pend_rope_len = CT;
	this->pend_angle = asin(TM / CT);
	this->pend_amplitude = this->pend_angle;
	this->pend_start_amplitude = this->pend_angle;

}

void node::phys_pend(float t) {

	coord carr = this->carrier;
	coord tick = this->ticker;

	if (!halt_pend && tick.y - this->size_y > FLOOR_HEIGHT && tick.y + this->size_y < ROOF_HEIGHT) {

		float omega = sqrt(G_CONST / this->pend_rope_len);
		//float attenuation = pow(M_E, -log(this->amplitude / (this->amplitude * 0.5)));
		float attenuation_k = 0.999;
		float phi = this->pend_amplitude * cos(omega * t + this->pend_start_amplitude);
		
		if (this->pend_amplitude < 0.07) {
			attenuation_k = 0.9999;
		}
		this->pend_amplitude *= attenuation_k;

		this->pend_angle = phi;

		float TM = this->pend_rope_len * sin(phi);
		float CM = this->pend_rope_len * sin(90 * M_PI / 180 - phi);	//1,57 rad

		coord T;
		T.x = carr.x - TM;
		T.y = carr.y - CM;

		this->ticker = T;

		//this->curr.x = T.x;
	}

}





//===================================================================================================
//===================================================================================================
//===================================================================================================
//dll.h
class dll {
public:
	bool permit_phys_move;
	bool phys_move_enabled;
	bool phys_pend_enabled;

	void append_front(int id);
	void append_back(int id);
	void destroy_list();

	int get_length();
	node* get_front();
	node* get_back();
	node* get_by_id(int id);
	node* get_by_pos(int pos);

	void display_forward();
	void display_backward();
	void display_objects();

	void draw_objects();
	void draw_ropes();
	void set_connection_points();
	void set_pend_params();

	void enable_phys_pend(float time);
	void enable_phys_move(float time, float delta_time);
	void check_collisions();
	void check_attenuation();

	dll() { front = NULL; back = NULL; length = 0; 
		phys_move_enabled = false; phys_pend_enabled = true; 
		permit_phys_move = false; }
	~dll() { destroy_list(); }


private:
	int length;
	int cargos_amount, blocks_amount;

	node *front;
	node *back;

};


//===================================================================================================
//===================================================================================================
//===================================================================================================
//dll.cpp
#include "stdafx.h"

using namespace std;

void dll::append_front(int x) {

	node *n = new node(x);

	if (front == NULL) {
		front = n;
		back = n;
	}
	else {
		front->prev_node = n;
		n->next_node = front;
		front = n;
	}

	this->length++;

	printf("fornt appended\n");
}

void dll::append_back(int x) {

	node *n = new node(x);

	if (back == NULL) {
		front = n;
		back = n;
	}
	else {
		back->next_node = n;
		n->prev_node = back;
		back = n;
	}

	this->length++;

	printf("back appended\n");
}

int dll::get_length() {
	return this->length;
}

node* dll::get_front() {
	return this->front;
}

node* dll::get_back() {
	return this->back;
}

node* dll::get_by_id(int id) {

	node *tmp = front;
	node *retval;

	while (tmp != NULL) {
		if (id == tmp->id) {
			retval = tmp;
			return retval;
		}
		tmp = tmp->next_node;
	}
	return NULL;
}

node* dll::get_by_pos(int pos) {

	node *tmp = front;
	node *retval;

	int i = 0;
	while (tmp != NULL) {
		if (i == pos) {
			retval = tmp;
			return retval;
		}
		i++;
		tmp = tmp->next_node;
	}
	return NULL;
}

void dll::display_forward() {
	node *tmp = front;

	printf("\n\nnodes(%d) in forward order:\n", this->length);

	while (tmp != NULL) {
		printf("%d   ", tmp->id);
		tmp = tmp->next_node;
	}
	printf("\n");
}

void dll::display_backward() {
	node *tmp = back;

	printf("\n\nnodes(%d) in reverse order:\n", this->length);

	while (tmp != NULL) {
		printf("%d   ", tmp->id);
		tmp = tmp->prev_node;
	}
	printf("\n");
}

void dll::display_objects() {
	node *tmp = front;

	printf("\n\nobjects(%d) in forward order:\n", this->length);

	while (NULL != tmp) {
		printf("(id:%d", tmp->id);
		if ('c' == tmp->type)
			printf("C, ");
		if ('b' == tmp->type)
			printf("B, ");
		printf("listpos:%d, mass:%d)   ", tmp->listpos, tmp->mass);
		tmp = tmp->next_node;
	}
	printf("\n");
}

void dll::draw_objects() {
	node *tmp = front;

	while (NULL != tmp) {
		tmp->draw_body();
		tmp = tmp->next_node;
	}
}

void dll::draw_ropes() {

	node *tmp = front;
	int i;

	for (i = 0; i < this->length; i++) {
		tmp->draw_ropes();
		tmp = tmp->next_node;

	}


}

void dll::destroy_list() {
	node *T = back;

	while (T != NULL) {
		node *T2 = T;
		T = T->prev_node;
		delete T2;
	}

	front = NULL;
	back = NULL;
}

void dll::enable_phys_pend(float t) {
	if (this->phys_pend_enabled) {
		node *tmp = front;

		while (NULL != tmp) {
			if ('c' == tmp->type) {
				tmp->phys_pend(t);
			}
			tmp = tmp->next_node;
		}
	}

	
}


void dll::enable_phys_move(float t, float dt) {
	if (this->phys_move_enabled && this->permit_phys_move) {
		
		node *tmp = front;

		float force = 0;
		float mass = 0;


		while (NULL != tmp) {
			tmp->phys_move(t);
			if (this->cargos_amount == 3 && tmp->type == 'c')
				tmp->explicit_accel = true;
			tmp = tmp->next_node;
		}

		/*while (NULL != tmp) {
			tmp->accel = -5;
			tmp->curr.y = RungeKutta4(tmp->curr.y, tmp->accel, t, dt);

			tmp = tmp->next_node;
		}*/

	}
}

void dll::check_collisions() {

	if (this->permit_phys_move) {
		node *tmp = front;

		while (NULL != tmp) {
			if (tmp->halt_move) {
				this->phys_move_enabled = false;
				tmp->halt_pend = true;
			}

			tmp = tmp->next_node;
		}
	}
}

void dll::check_attenuation() {

	node *tmp = front;
	int mark = 0;
	while (NULL != tmp) {
		if ('c' == tmp->type) {
			if (tmp->pend_amplitude < 0.07) {

				//tmp->halt_pend = true;
				tmp->curr.x = tmp->ticker.x;
				mark++;
			}
			else {
				tmp->curr = tmp->ticker;
			}
		}
		tmp = tmp->next_node;
	}

	if (mark == this->cargos_amount) {
		this->permit_phys_move = true;
	}
}

void dll::set_connection_points() {
	node *tmp = front;

	//count by type
	while (NULL != tmp) {
		if ('c' == tmp->type) this->cargos_amount++;
		if ('b' == tmp->type) this->blocks_amount++;
		tmp = tmp->next_node;
	}

	tmp = front;

	while (NULL != tmp) {
		tmp->set_list_size(this->length, this->cargos_amount, this->blocks_amount);
		if ('c' == tmp->type) {
			int cp_center_id = tmp->cp_center_id;
			tmp->cp_center_node = this->get_by_id(cp_center_id);
		}
		if ('b' == tmp->type) {
			tmp->cp_left_node = this->get_by_id(tmp->cp_left_id);
			tmp->cp_center_node = this->get_by_id(tmp->cp_center_id);
			tmp->cp_right_node = this->get_by_id(tmp->cp_right_id);
		}
		tmp = tmp->next_node;
	}
}

void dll::set_pend_params() {
	node *tmp = front;

	while (NULL != tmp) {
		if ('c' == tmp->type) {
			tmp->count_pend_param();
		}
		tmp = tmp->next_node;
	}
}

