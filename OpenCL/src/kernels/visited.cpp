#include "kernel2.h"

inline void ele_init_cl(ele_cl* present_ele, output_type_cl* x_cl, float f_cl, change_cl* d_cl) {
	for (int i = 0; i < 3; i++)present_ele->position[i] = x_cl->position[i];
	for (int i = 0; i < 4; i++)present_ele->orientation[i] = x_cl->orientation[i];
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)present_ele->lig_torsion[i] = x_cl->lig_torsion[i];
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)present_ele->flex_torsion[i] = x_cl->flex_torsion[i];
	present_ele->energy = f_cl;
	present_ele->d_zero = 0;
	present_ele->d_positive = 0;

	const long ONE  = 1;
	long bitMask    = 0;

	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if (d_cl->position[i] == 0)
			present_ele->d_zero |= bitMask;
		else if (d_cl->position[i] > 0)
			present_ele->d_positive |= bitMask;
	}
	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if (d_cl->orientation[i] == 0)
			present_ele->d_zero |= bitMask;
		else if (d_cl->orientation[i] > 0)
			present_ele->d_positive |= bitMask;
	}
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++) {
		bitMask = ONE << i;
		if (d_cl->lig_torsion[i] == 0)
			present_ele->d_zero |= bitMask;
		else if (d_cl->lig_torsion[i] > 0)
			present_ele->d_positive |= bitMask;
	}
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
		bitMask = ONE << i;
		if (d_cl->flex_torsion[i] == 0)
			present_ele->d_zero |= bitMask;
		else if (d_cl->flex_torsion[i] > 0)
			present_ele->d_positive |= bitMask;
	}
}


inline void circularvisited_init_cl(individual_container* list, output_type_cl* conf_v) {
	list->n_variable = conf_v->lig_torsion_size;
	list->p = 0;
	list->full = false;	
	list->counter = 0;
	
}

inline float dist2_cl(individual_container* list, output_type_cl* now_x, int neighbour) {
	float out = 0;
	for (int i = 0; i < 3; i++) {
		double d = list->list_cl[neighbour]->position[i] - now_x->position[i];
		out += d * d;
	}
	for (int i = 0; i < 4; i++) {
		double d = list->list_cl[neighbour]->orientation[i] - now_x->orientation[i];
		out += d * d;
	}
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++) {
		double d = list->list_cl[neighbour]->lig_torsion[i] - now_x->lig_torsion[i];
		out += d * d;
	}
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++) {
		double d = list->list_cl[neighbour]->flex_torsion[i] - now_x->flex_torsion[i];
		out += d * d;
	}
	return out;
}

inline float dist2_3D_cl(individual_container* list, output_type_cl* now_x, int neighbour) {
	double d, out = 0;
	for (int i = 0; i < 3; i++) {
		d = list->list_cl[neighbour]->position[i] - now_x->position[i];
		out += d * d;
	}
	return out;
}

inline bool add_to_individual_buffer(individual_container* list, output_type_cl* conf_v, float f, change_cl* change_v) {
	list->tempf = f;
	ele_cl* e;
	ele_init_cl(e, conf_v, f, change_v);
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
	return true;
}

inline bool check_cl(individual_container* list, output_type_cl* now_x, float now_f, change_cl* now_d,int neighbour) {
	bool newXBigger, newYBigger;
	const long ONE = 1;
	long bitMask;
	newYBigger = (now_f - list->list_cl[neighbour]->energy) > 0;

	for (int i = 0; i < 3; i++) {
		bitMask = ONE << i;
		if ((list->list_cl[neighbour]->d_zero & bitMask) || !(now_d->position[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->position[i] > 0;
			const bool dPositive = list->list_cl[neighbour]->d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->position[i] - list->list_cl[neighbour]->position[i]) > 0;
				if (!(nowPositive ^ newXBigger ^ newYBigger)) {
					return true;
				}
				else {
					return false;
				}
			}
		}
	}
	for (int i = 0; i < 4; i++) {
		bitMask = ONE << i;
		if ((list->list_cl[neighbour]->d_zero & bitMask) || !(now_d->orientation[i - 1])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->orientation[i - 1] > 0;
			const bool dPositive = list->list_cl[neighbour]->d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->orientation[i] - list->list_cl[neighbour]->orientation[i]) > 0;
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
		if ((list->list_cl[neighbour]->d_zero & bitMask) || !(now_d->lig_torsion[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->lig_torsion[i] > 0;
			const bool dPositive = list->list_cl[neighbour]->d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->lig_torsion[i] - list->list_cl[neighbour]->lig_torsion[i]) > 0;
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
		if ((list->list_cl[neighbour]->d_zero & bitMask) || !(now_d->flex_torsion[i])) {
			continue;
		}
		else {
			const bool nowPositive = now_d->flex_torsion[i] > 0;
			const bool dPositive = list->list_cl[neighbour]->d_positive & bitMask;
			if (nowPositive ^ dPositive) {
				continue;
			}
			else {
				newXBigger = (now_x->flex_torsion[i] - list->list_cl[neighbour]->flex_torsion[i]) > 0;
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

inline int global_interesting_cl(individual_container* list, output_type_cl* now_x, float now_f, change_cl* now_d, int excluded) {
	ele_cl* nearbyPoints[MAX_SIZE_OF_LIST];
	float distances[MAX_SIZE_OF_LIST];
	float dist[MAX_SIZE_OF_LIST];
	bool notPicked[MAX_SIZE_OF_LIST];


	

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
			dist[i] = dist2_cl(list, now_x, i);
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
			if (check_cl(list, now_x, now_f, now_d, p))
				return -1; //i.e. interesting
		}
		return i;
	}
}
