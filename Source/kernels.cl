#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include <glm/glm.hpp>
typedef struct
{
  glm::vec4 v0;
  glm::vec4 v1;
  glm::vec4 v2;
  glm::vec4 normal;
  glm::vec3 color;
} triangle;

kernel void draw(global uint *screen_buffer, global triangle *triangle, global float *r)
{         /* accumulated magnitudes of velocity for each cell */

  const short x = get_global_id(0);
  const short y = get_global_id(1);


}