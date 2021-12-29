#include "kernel2.h"

//inline float dist2_cl(individual_container* list, output_type_cl* now_x, int neighbour) {
//	float out = 0;
//	for (int i = 0; i < 3; i++) {
//		float d = list->list_cl[neighbour].position[i] - now_x->position[i];
//		out += d * d;
//	}
//	for (int i = 0; i < 4; i++) {
//		float d = list->list_cl[neighbour].orientation[i] - now_x->orientation[i];
//		out += d * d;
//	}
//	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++) {
//		float d = list->list_cl[neighbour].lig_torsion[i] - now_x->lig_torsion[i];
//		out += d * d;
//	}
//	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
//		float d = list->list_cl[neighbour].flex_torsion[i] - now_x->flex_torsion[i];
//		out += d * d;
//	}
//	return out;
//}

inline float dist2_g(ele_cl list, output_type_cl* now_x) {
	float out = 0;
	for (int i = 0; i < 3; i++) {
		float d = list.position[i] - now_x->position[i];
		out += d * d;
	}
	for (int i = 0; i < 4; i++) {
		float d = list.orientation[i] - now_x->orientation[i];
		out += d * d;
	}
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++) {
		float d = list.lig_torsion[i] - now_x->lig_torsion[i];
		out += d * d;
	}
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
		float d = list.flex_torsion[i] - now_x->flex_torsion[i];
		out += d * d;
	}
	return out;
}

inline float dist2_3D_cl(ele_cl list_cl, output_type_cl* now_x) {
	float d, out = 0;
	for (int i = 0; i < 3; i++) {
		d = list_cl.position[i] - now_x->position[i];
		out += d * d;
	}
	return out;
}

int get_n(int thread, int search_depth) {
	for (int i = 0;; i += 3) {
		if (thread * search_depth < (2<<i))
			return i;
	}
}

inline ele_cl ele_init_cl(output_type_cl* x_cl, float f_cl, change_cl* d_cl) {
	ele_cl present_ele;
	for (int i = 0; i < 3; i++)present_ele.position[i] = x_cl->position[i];
	for (int i = 0; i < 4; i++)present_ele.orientation[i] = x_cl->orientation[i];
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)present_ele.lig_torsion[i] = x_cl->lig_torsion[i];
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)present_ele.flex_torsion[i] = x_cl->flex_torsion[i];
	present_ele.energy = f_cl;
	present_ele.d_zero = 0;
	present_ele.d_positive = 0;

	const long ONE  = 1;
	long bitMask    = 0;

	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if (d_cl->position[i] == 0)
			present_ele.d_zero |= bitMask;
		else if (d_cl->position[i] > 0)
			present_ele.d_positive |= bitMask;
	}
	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if (d_cl->orientation[i] == 0)
			present_ele.d_zero |= bitMask;
		else if (d_cl->orientation[i] > 0)
			present_ele.d_positive |= bitMask;
	}
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++) {
		bitMask = ONE << i;
		if (d_cl->lig_torsion[i] == 0)
			present_ele.d_zero |= bitMask;
		else if (d_cl->lig_torsion[i] > 0)
			present_ele.d_positive |= bitMask;
	}
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
		bitMask = ONE << i;
		if (d_cl->flex_torsion[i] == 0)
			present_ele.d_zero |= bitMask;
		else if (d_cl->flex_torsion[i] > 0)
			present_ele.d_positive |= bitMask;
	}
	return present_ele;
}


inline vec3_cl origin_init(float center_x, float center_y, float center_z) {
	vec3_cl origin;
	origin.x = center_x;
	origin.y = center_y;
	origin.z = center_z;
	return origin;
}

inline vec3_cl boxsize_init(float size_x, float size_y, float size_z) {
	vec3_cl boxsize;
	boxsize.x = size_x;
	boxsize.y = size_y;
	boxsize.z = size_z;
	return boxsize;
}


inline int set_Binary_Id_cl(int thread, int search_depth, output_type_cl* x, vec3_cl origin, vec3_cl dimension) {
	int n = get_n(thread, search_depth) - 3;
	int binary_id = 0;
	int oct[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	for (int i = 0; i < n / 3; i++) {
		if (x->position[0] >= origin.x) oct[i] |= 1;
		if (x->position[1] >= origin.y) oct[i] |= 2;
		if (x->position[2] >= origin.z) oct[i] |= 4;
		switch (oct[i])
		{
		case 0: origin.x = origin.x - dimension.x / (4 << (2 * i));
			origin.y = origin.y - dimension.y / (4 << (2 * i));
			origin.z = origin.z - dimension.z / (4 << (2 * i));
			break;
		case 1: origin.x = origin.x - dimension.x / (4 << (2 * i));
			origin.y = origin.y - dimension.y / (4 << (2 * i));
			origin.z = origin.z + dimension.z / (4 << (2 * i));
			break;
		case 2: origin.x = origin.x - dimension.x / (4 << (2 * i));
			origin.y = origin.y + dimension.y / (4 << (2 * i));
			origin.z = origin.z - dimension.z / (4 << (2 * i));
			break;
		case 3: origin.x = origin.x - dimension.x / (4 << (2 * i));
			origin.y = origin.y + dimension.y / (4 << (2 * i));
			origin.z = origin.z + dimension.z / (4 << (2 * i));
			break;
		case 4: origin.x = origin.x + dimension.x / (4 << (2 * i));
			origin.y = origin.y - dimension.y / (4 << (2 * i));
			origin.z = origin.z - dimension.z / (4 << (2 * i));
			break;
		case 5: origin.x = origin.x + dimension.x / (4 << (2 * i));
			origin.y = origin.y - dimension.y / (4 << (2 * i));
			origin.z = origin.z + dimension.z / (4 << (2 * i));
			break;
		case 6: origin.x = origin.x + dimension.x / (4 << (2 * i));
			origin.y = origin.y + dimension.y / (4 << (2 * i));
			origin.z = origin.z - dimension.z / (4 << (2 * i));
			break;
		case 7: origin.x = origin.x + dimension.x / (4 << (2 * i));
			origin.y = origin.y + dimension.y / (4 << (2 * i));
			origin.z = origin.z + dimension.z / (4 << (2 * i));
			break;
		}
	}
	for (int i = 0; i < n / 3; i++) {
		binary_id += (oct[i] << (n - 3 - 3 * i));
	}
	return binary_id;
}


inline void circularvisited_init_cl(individual_container* list, output_type_cl* conf_v) {
	list->n_variable = conf_v->lig_torsion_size;
	list->p = 0;
	list->full = false;
	list->counter = 0;

}

inline void global_init_cl(global_container* g_container, output_type_cl* conf_v) {
	g_container->n_variable = conf_v->lig_torsion_size;
	g_container->binary_ID = 0;
}


inline void add_to_individual_buffer(individual_container* list, output_type_cl* conf_v, float f, change_cl* change_v) {
	//circularvisited_init_cl(list, conf_v);
	list->tempf = f;
	ele_cl e = ele_init_cl(conf_v, f, change_v);	
	if (!list->full) {
		list->list_cl[list->counter] = e;
		list->counter++;
		if (list->counter >= 5 * list->n_variable) {
			list->full = true;
			list->p = 0;
		}
	}
	else {
		list->list_cl[list->p] = e;
		list->p = (list->p + 1) % (5 * list->n_variable);
	}
}

inline void add_to_global_buffer(global_container* global_c, int thread, int search_depth, 
								 output_type_cl* conf_v, float f, change_cl* change_v, 
								 vec3_cl origin, vec3_cl dimension, ele_cl* ptr) {
	//global_init_cl(global_c, conf_v);
	global_c->tempf = f;
	ele_cl e = ele_init_cl(conf_v, f, change_v);;
	global_c->binary_ID = set_Binary_Id_cl(thread, search_depth, conf_v, origin, dimension);

	//add e to global buffer according to binary_ID
	*(ptr + global_c->binary_ID * 8) = e;
	global_c->counter[global_c->binary_ID]++;
	if (global_c->counter[global_c->binary_ID] >= 8) {
		global_c->counter[global_c->binary_ID] = 7;
	}

}

inline void getPointsWithinCutoff_cl(float cutoff, global_container* global_c, int thread, int search_depth, output_type_cl* point, ele_cl* results, float* distances, vec3_cl origin, vec3_cl dimension, ele_cl* ptr) {
	global_c->binary_ID = set_Binary_Id_cl(thread, search_depth, point, origin, dimension);
	if (dimension.x <= cutoff && dimension.y <= cutoff && dimension.z <= cutoff) {
		//取global_c->binary_ID对应的global buffer里面的点，结果存放在result数组里面
		for (int i = 0; i < global_c->counter[global_c->binary_ID]; i++) {
			results[i] = *(ptr + global_c->binary_ID * 8 + i);
		}
		//并把这些点与待检测点的n维欧式距离添加到distances数组里面
		for (int i = 0; i < global_c->counter[global_c->binary_ID]; i++) {
			distances[i] = dist2_g(results[i], point);
		}
	}
	else {
		for (int i = 0; i < global_c->counter[global_c->binary_ID]; i++) {
			float dist2_3D = dist2_3D_cl(*(ptr + global_c->binary_ID * 8 + i),point);
			if (dist2_3D <= cutoff) {
				results[i] = *(ptr + global_c->binary_ID * 8 + i);
				distances[i] = dist2_g(results[i], point);
			}
		}
		
	}
	
}


inline bool check_cl(ele_cl* list_cl, output_type_cl* now_x, float now_f, change_cl* now_d,int neighbour) {
	bool newXBigger, newYBigger;
	const long ONE = 1;
	long bitMask;
	newYBigger = (now_f - list_cl[neighbour].energy) > 0;

	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if ((list_cl[neighbour].d_zero & bitMask) || !(now_d->position[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->position[i] > 0;
			const bool dPositive = list_cl[neighbour].d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->position[i] - list_cl[neighbour].position[i]) > 0;
				if (!(nowPositive ^ newXBigger ^ newYBigger)) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if ((list_cl[neighbour].d_zero & bitMask) || !(now_d->orientation[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->orientation[i] > 0;
			const bool dPositive = list_cl[neighbour].d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->orientation[i] - list_cl[neighbour].orientation[i]) > 0;
				if (!(nowPositive ^ newXBigger ^ newYBigger)) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++) {
		bitMask = ONE << i;
		if ((list_cl[neighbour].d_zero & bitMask) || !(now_d->lig_torsion[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->lig_torsion[i] > 0;
			const bool dPositive = list_cl[neighbour].d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->lig_torsion[i] - list_cl[neighbour].lig_torsion[i]) > 0;
				if (!(nowPositive ^ newXBigger ^ newYBigger)) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
		if ((list_cl[neighbour].d_zero & bitMask) || !(now_d->flex_torsion[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->flex_torsion[i] > 0;
			const bool dPositive = list_cl[neighbour].d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->flex_torsion[i] - list_cl[neighbour].flex_torsion[i]) > 0;
				if (!(nowPositive ^ newXBigger ^ newYBigger)) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}
}



inline int global_interesting_cl(global_container* global_mem, output_type_cl* now_x, float now_f, change_cl* now_d, int excluded, int thread, int search_depth, vec3_cl origin, vec3_cl dimension, ele_cl* ptr) {
	//ele_cl nearbyPoints[MAX_SIZE_OF_LIST];
	//float distances[MAX_SIZE_OF_LIST];
	//float dist[MAX_SIZE_OF_LIST];
	//int len = global_mem->counter;
	bool notYetChecked[MAX_NUM_OF_RESULTS];
	for (int i = 0; i < MAX_NUM_OF_RESULTS; i++) {
		notYetChecked[i] = true;
	}
	const float CUTOFF = 5.0;
	
	getPointsWithinCutoff_cl(CUTOFF, global_mem, thread, search_depth, now_x, global_mem->nearbyPoints, global_mem->distances, origin, dimension, ptr);

	int nearbyPoints_size = sizeof(global_mem->nearbyPoints) / sizeof(ele_cl);
	const int grandMaxCheck = 1 * global_mem->n_variable; //1N in this case
	const int maxCheck = (nearbyPoints_size <= grandMaxCheck) ? nearbyPoints_size : grandMaxCheck;

	double min = 1e10;
	int i = 0; //counts checked done so far
	int p; //pointer to current nearest point
	for (; i < maxCheck; i++) {
		min = 1e10;
		for (int j = 0; j < nearbyPoints_size; j++) {
			if (notYetChecked[j] && (global_mem->distances[j] <= min)) {
				p = j;
				min = global_mem->distances[p];
			}
		}
		notYetChecked[p] = false;

		if (check_cl(global_mem->nearbyPoints, now_x, now_f, now_d, p)) {
			return -1; //i.e. return success
		}
	}
	return i;

}

inline int individual_interesting_cl(individual_container* list, output_type_cl* now_x, float now_f, change_cl* now_d, int excluded) {
	int len = list->counter;
	float dist[MAX_SIZE_OF_LIST];
	bool notPicked[MAX_SIZE_OF_LIST];
	
	if (len == 0) {
		return -1; //i.e. interesting
	}
	else {
		if (!list->full) {
			return -1; //i.e. interesting
		}
		for (int i = 0; i < len; i++) {
			notPicked[i] = true;
		}
		for (int i = 0; i < len; i++) {
			dist[i] = dist2_g(list->list_cl[i], now_x);
		}
		double min = 1e10;
		int p = 0;
		const int maxCheck = 4 * list->n_variable - excluded;
		int i = 0;
		for (; i < maxCheck; i++) {
			min = 1e10;
			for (int j = 0; j < len; j++) {
				if (notPicked[j] && (dist[j] < min)) {
					p = j;
					min = dist[j];
				}
			}
			notPicked[p] = false;
			if (check_cl(list->list_cl, now_x, now_f, now_d, p))
				return -1; //i.e. interesting
		}
		return i;
	}
}
