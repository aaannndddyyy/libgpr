#include "globals.h"

/* convert RGB to HSL */
void rgb_to_hsl(unsigned char R,
				unsigned char G,
				unsigned char B,
				float * H,
				float * S,
				float * L)
{
	float var_R = R / 255.0f;
	float var_G = G / 255.0f;
	float var_B = B / 255.0f;
	float var_Min, var_Max, del_Max, del_R, del_G, del_B;

	var_Min = var_R;
	if (var_G < var_Min) var_Min = var_G;
	if (var_B < var_Min) var_Min = var_B;

	var_Max = var_R;
	if (var_G > var_Max) var_Max = var_G;
	if (var_B > var_Max) var_Max = var_B;

	del_Max = var_Max - var_Min;

	*L = (var_Max + var_Min) * 0.5f;

	if (del_Max == 0) {
		*H = 0;
		*S = 0;
	}
	else {
		if ( *L < 0.5f ) {
			*S = del_Max / (var_Max + var_Min);
		}
		else {
			*S = del_Max / (2.0f - var_Max - var_Min);
		}

		del_R = (((var_Max - var_R) / 6.0f) + (del_Max / 2.0f)) / del_Max;
		del_G = (((var_Max - var_G) / 6.0f) + (del_Max / 2.0f)) / del_Max;
		del_B = (((var_Max - var_B) / 6.0f) + (del_Max / 2.0f)) / del_Max;

		if (var_R == var_Max) {
			*H = del_B - del_G;
		}
		else if ( var_G == var_Max ) {
			*H = (1 / 3.0f) + del_R - del_B;
		}
		else if (var_B == var_Max) {
			*H = (2.0f / 3.0f) + del_G - del_R;
		}

		if (*H < 0) *H += 1;
		if (*H > 1) *H -= 1;
	}
}

static float hue_to_rgb(float v1, float v2, float vH )
{
	if (vH < 0) vH += 1;
	if (vH > 1) vH -= 1;
	if ((6 * vH) < 1) return v1 + ( v2 - v1 ) * 6 * vH;
	if ((2 * vH) < 1) return v2;
	if ((3 * vH) < 2) return v1 + ( v2 - v1 ) * ( ( 2 / 3 ) - vH ) * 6;
	return v1;
}

void hsl_to_rgb(float H,
				float S,
				float L,
				unsigned char * R,
				unsigned char * G,
				unsigned char * B)
{
	float var_1, var_2;

	if (S == 0) {
		*R = (unsigned char)(L * 255);
		*G = (unsigned char)(L * 255);
		*B = (unsigned char)(L * 255);
	}
	else {
		if (L < 0.5f) {
			var_2 = L * ( 1 + S );
		}
		else {
			var_2 = (L + S) - (S * L);
		}

		var_1 = 2 * L - var_2;

		*R = (unsigned char)(255 * hue_to_rgb(var_1, var_2, H + (1.0f / 3.0f)));
		*G = (unsigned char)(255 * hue_to_rgb(var_1, var_2, H));
		*B = (unsigned char)(255 * hue_to_rgb(var_1, var_2, H - (1.0f / 3.0f)));
	}
}
