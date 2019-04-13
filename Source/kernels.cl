#pragma OPENCL EXTENSION cl_khr_fp64 : enable

constant float focal_length = 500.0;
constant float3 indirect_light = (float3)(0.5f, 0.5f, 0.5f);
#define SCREEN_WIDTH 1400
#define SCREEN_HEIGHT 1400

typedef struct
{
  int4 v0;
  int4 v1;
  int4 v2;
  float4 normal;
  float3 color;
} triangle;

kernel void draw(global uint  *screen_buffer, global triangle *triangles, 
				 global float *rot_matrix,    global float4 camera_pos)
{         /* accumulated magnitudes of velocity for each cell */
  const short x = get_global_id(0);
  const short y = get_global_id(1);


}