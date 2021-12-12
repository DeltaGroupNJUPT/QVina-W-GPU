#include "kernel2.h"

inline void ele_init_cl(std_vector_cl* x_cl, float f_cl, change_cl* d_cl) {
	ele_cl* present_ele;
	for (int i = 0; i < 3; i++)present_ele->position[i] = x_cl->position[i];
	for (int i = 0; i < 3; i++)present_ele->orientation[i] = x_cl->orientation[i];
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)present_ele->lig_torsion[i] = x_cl->lig_torsion[i];
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)present_ele->flex_torsion[i] = x_cl->flex_torsion[i];
	present_ele->energy = f_cl;

	present_ele->d_zero     = 0x0000000000000000;
	present_ele->d_positive = 0x0000000000000000;
	const long ONE  = 0x0000000000000001;
	long bitMask    = 0x0000000000000000;

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

inline void circularvisited_cl_init(circularvisited_cl* visited) {
	visited->n_variable = 0;
	visited->p = 0;
	visited->full = false;
}

