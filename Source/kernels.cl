#pragma OPENCL EXTENSION cl_khr_fp64 : enable

constant float3 indirect_light = (float3)(0.5f, 0.5f, 0.5f);
constant float3 light_color    = (float3) (14.0f, 14.0f, 14.0f);
#define SCREEN_WIDTH 4608
#define SCREEN_HEIGHT 4608

/////READ ONLY BUFFERS

typedef struct  {
  float3 position;
  float distance;
  int triangle_index;
} Intersection;

inline float det(float3 M[3]) {
	return M[0].x * (M[1].y * M[2].z - M[1].z * M[2].y) -
		     M[0].y * (M[1].x * M[2].z - M[1].z * M[2].x) +
		     M[0].z * (M[1].x * M[2].y - M[1].y * M[2].x);
}

inline void PutPixelSDL(global uint *screen_buffer, int x, int y, float3 colour) {
  // if(x<0 || x>=SCREEN_WIDTH || y<0 || y>=SCREEN_HEIGHT)  {
  //   printf("apa\n");
  //   return;
  // }
  uint3 rgb = convert_uint3(min(max(255*colour, 0.f), 255.f));
  screen_buffer[y*SCREEN_WIDTH+x] = (128<<24) + (rgb.x<<16) + (rgb.y<<8) + rgb.z;
}


bool closest_intersection(float3 start, float3 d, local float3 *triangle_vertexes, private Intersection* closest_intersection, int triangle_n) {
  // Set closest intersection to be the max float value
  float current_t = MAXFLOAT;
  // Make 4D ray into 3D ray
  for (uint i = 0; i < triangle_n; i++) {
    // Define two corners of triangle relative to the other corner
    const float3 v0 = triangle_vertexes[i*3];
    const float3 v1 = triangle_vertexes[i*3+1];
    const float3 v2 = triangle_vertexes[i*3+2];

    const float3 e1 = (float3) (v1.x-v0.x,    v1.y-v0.y,    v1.z-v0.z);
    const float3 e2 = (float3) (v2.x-v0.x,    v2.y-v0.y,    v2.z-v0.z);
    const float3 b  = (float3) (start.x-v0.x, start.y-v0.y, start.z-v0.z);

    // Cramers, might be det repeated computation..?
    const float3 A[3]  = {-d, e1, e2};
    const float3 A0[3] = {b,  e1, e2};
    const float3 A1[3] = {-d, b,  e2};
    const float3 A2[3] = {-d, e1, b};

    const float detA  = det(A);
    const float detA0 = det(A0);
    const float detA1 = det(A1);
    const float detA2 = det(A2); 

    float3 x = (float3) (detA0/detA, detA1/detA, detA2/detA);

    // If ray goes through triangle, and is the closest triangle
    if (x.x >= 0 && x.y >= 0 && x.z >= 0 && (x.y + x.z) <= 1 && x.x < current_t) {
      float3 position = ((float3) (v0.x, v0.y, v0.z)) + (x.y * e1) + (x.z * e2);

      closest_intersection->position       = (float3) (position.x, position.y, position.z);
      float3 dist_vec                      = x.x*d;
      closest_intersection->distance       = native_sqrt(dist_vec.x*dist_vec.x + dist_vec.y*dist_vec.y + dist_vec.z*dist_vec.z);
      closest_intersection->triangle_index = i;
      current_t                            = x.x;
    }
  }
  if (current_t == MAXFLOAT) return false;
  return true;
}


float3 direct_light(const Intersection intersection, local float3 *triangle_vertexes, local float3 *triangle_normals, float3 light_pos, int triangle_n) {

  // Vector from the light to the point of intersection
  float3 r = light_pos - intersection.position;
  // Distance of the checked point to the light source
  float radius = native_sqrt(r.x*r.x + r.y*r.y + r.z*r.z);;

  Intersection obstacle_intersection;
  float threshold = 0.001f;
  float3 intersect_pos = intersection.position + (float3) (r.x * threshold, r.y * threshold, r.z * threshold);

  if (closest_intersection(intersect_pos, r, triangle_vertexes, &obstacle_intersection, triangle_n)) {
    if (obstacle_intersection.distance < radius) return (float3)(0.0f, 0.0f, 0.0f);
  }

  // Get the normal of the triangle that the light has hit
  float3 n = triangle_normals[intersection.triangle_index];
  // Intensity of the colour, based on the distance from the light
  float3 D = (light_color * max(dot(r, n) , 0.0f)) / (4 * ((float)M_PI) * radius * radius);

  return D;
}


kernel void draw(global uint  *screen_buffer,    global float3 *triangle_vertexes,   global float3 *triangle_normals,
				 global float3 *triangle_colors, global float3 *rot_matrix,           float3 camera_pos, float3 light_pos, 
				 int triangle_n, float focal_length, local float3 *LOC_triangle_vertexes,  local float3 *LOC_triangle_normals,
				 local float3 *LOC_triangle_colors)
{         /* accumulated magnitudes of velocity for each cell */
  const short x = get_global_id(0);
  const short y = get_global_id(1);


  event_t e = async_work_group_copy(LOC_triangle_vertexes, triangle_vertexes, triangle_n*3, 0);
  e         = async_work_group_copy(LOC_triangle_normals,  triangle_normals,  triangle_n,  0);
  e         = async_work_group_copy(LOC_triangle_colors,   triangle_colors,   triangle_n,  0);
  wait_group_events(3, &e);


  // Declare ray for given position on the screen. Rotate ray by current view angle
  float3 d = (float3) (x - SCREEN_WIDTH/2.0, y - SCREEN_HEIGHT/2.0, focal_length);
  d        = (float3) (dot(rot_matrix[0], d), dot(rot_matrix[1], d), dot(rot_matrix[2], d));

  // Find intersection point with closest geometry. If no intersection, paint the abyss
  Intersection intersection;
  if (closest_intersection(camera_pos, d, LOC_triangle_vertexes, &intersection, triangle_n)) {
    const float3 p = LOC_triangle_colors[intersection.triangle_index];
    const float3 final_color = p*(direct_light(intersection, LOC_triangle_vertexes, LOC_triangle_normals, light_pos, triangle_n) + indirect_light);
  	PutPixelSDL(screen_buffer, x, y, final_color);

  } else {
    // Otherwise draw black
  	PutPixelSDL(screen_buffer, x, y, (float3) (0.0f, 0.0f, 0.0f));
  }
}

uint3 getRGB(uint pixel)  {
  return (uint3) ((uint)((pixel >> 16) & 128), (uint)((pixel >> 8) & 128), (uint)(pixel & 128));
}

kernel void average_pixels(global uint *screen_buffer)  {
  const short x = get_global_id(0);
  const short y = get_global_id(1);

  const short nx = get_global_size(1);
  
  uint3 surrounding_cell_total = (uint3) (0, 0, 0);
  // surrounding_cell_total  += getRGB(screen_buffer[(y*3)*SCREEN_WIDTH+(x*3)]);
  // surrounding_cell_total  += getRGB(screen_buffer[(y*3)*SCREEN_WIDTH+(x*3+1)]);
  // surrounding_cell_total  += getRGB(screen_buffer[(y*3)*SCREEN_WIDTH+(x*3+2)]);

  // surrounding_cell_total  += getRGB(screen_buffer[(y*3+1)*SCREEN_WIDTH+(x*3)]);
  surrounding_cell_total  += getRGB(screen_buffer[(y*3+1)*SCREEN_WIDTH+(x*3+1)]);
  // surrounding_cell_total  += getRGB(screen_buffer[(y*3+1)*SCREEN_WIDTH+(x*3+2)]);

  // surrounding_cell_total  += getRGB(screen_buffer[(y*3+2)*SCREEN_WIDTH+(x*3)]);
  // surrounding_cell_total  += getRGB(screen_buffer[(y*3+2)*SCREEN_WIDTH+(x*3+1)]);
  // surrounding_cell_total  += getRGB(screen_buffer[(y*3+2)*SCREEN_WIDTH+(x*3+2)]);
  // printf("Averaged pixel Val = %d\n", surrounding_cell_total/9);

  if (x == 300 && y == 300)
    printf("screen_buffer %ud\n", screen_buffer[(y*3+1)*SCREEN_WIDTH+(x*3+1)]);

  surrounding_cell_total /= 1;

  // surrounding_cell_total = min(max(surrounding_cell_total, (uint3)0), (uint3)255);
  // screen_buffer[y*nx+x] = screen_buffer[(y*3+1)*SCREEN_WIDTH+(x*3+1)];
  screen_buffer[y*nx+x] = (128<<24) + (surrounding_cell_total.x<<16) + (surrounding_cell_total.y<<8)
                                                                      + surrounding_cell_total.z;
  if (x == 300 && y == 300)
    printf("screen_buffer %ud\n", (128<<24) + (surrounding_cell_total.x<<16) + (surrounding_cell_total.y<<8)
                                                                      + surrounding_cell_total.z);

}


