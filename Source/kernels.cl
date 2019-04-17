#pragma OPENCL EXTENSION cl_khr_fp64 : enable

constant float3 indirect_light = (float3)(0.5f, 0.5f, 0.5f);
constant float3 light_color    = (float3) (14.0f, 14.0f, 14.0f);
#define SCREEN_WIDTH 1536.0f
#define SCREEN_HEIGHT 1536.0f

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
  screen_buffer[y*(short)SCREEN_WIDTH+x] = (128<<24) + (rgb.x<<16) + (rgb.y<<8) + rgb.z;
}
inline float3 rnd(float3 seed, float3 range) {
  // return seed;
  return fmod(7*seed, range)-range/2.0f;
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

    const float detA  = det(A);
    const float detA0 = det(A0);

    float t = detA0/detA;

    // If ray goes through triangle, and is the closest triangle
    if ( t < current_t && t >= 0 ) {
      const float3 A1[3] = {-d, b,  e2};
      const float3 A2[3] = {-d, e1, b};

      const float detA1 = det(A1);
      const float detA2 = det(A2);
      float u = detA1/detA;
      float v = detA2/detA; 

      if (u >= 0 && v >= 0 && (u+v) <= 1) {
        float3 position = ((float3) (v0.x, v0.y, v0.z)) + (u * e1) + (v * e2);

        closest_intersection->position       = (float3) (position.x, position.y, position.z);
        float3 dist_vec                      = t*d;
        closest_intersection->distance       = native_sqrt(dist_vec.x*dist_vec.x + dist_vec.y*dist_vec.y + dist_vec.z*dist_vec.z);
        closest_intersection->triangle_index = i;
        current_t                            = t;
      }

    }
  }
  if (current_t == MAXFLOAT) return false;
  return true;
}

bool in_shadow(float3 start, float3 d, local float3 *triangle_vertexes, float radius_sq, int triangle_n) {
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

    const float detA  = det(A);
    const float detA0 = det(A0);

    float t = detA0/detA;

    // If ray goes through triangle, and is the closest triangle
    if (t >= 0 ) {
      const float3 A1[3] = {-d, b,  e2};
      const float3 A2[3] = {-d, e1, b};

      const float detA1 = det(A1);
      const float detA2 = det(A2);
      float u = detA1/detA;
      float v = detA2/detA; 
      float3 dist_vec            = t*d;
      float intersect_dist       = dist_vec.x*dist_vec.x + dist_vec.y*dist_vec.y + dist_vec.z*dist_vec.z;


      if (u >= 0 && v >= 0 && (u+v) <= 1 && intersect_dist < radius_sq) {
        return true;
      }

    }
  }
  return false;
}

float3 direct_light(const Intersection intersection, local float3 *triangle_vertexes, local float3 *triangle_normals, 
                    float3 light_pos, int triangle_n, float3 intersect_normal) {

  //Declare colour for point to be 0
  float3 total_colour = (float3) 0.0f;

  //Get vector from intersection point to light position, and its magnitude
  float3 dir = light_pos - intersection.position;
  float radius_sq = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;

  //Declare threshold to get intersection position that is not going to intersect with own triangle
  const float threshold = 0.001f;
  float3 start = intersection.position + threshold*dir;

  const float soft_shadows = 10.0f;
  const float3 soft_shadow_color_step = (float3)(0.9f/soft_shadows);

  // Check parallel ghost surfaces for soft triangles
  for (int i = 0; i < soft_shadows; i++)  {

    float3 ghost_dir = dir + rnd(intersection.position, 0.05)+threshold*dir;
    float ghost_radius_sq = ghost_dir.x*ghost_dir.x + ghost_dir.y*ghost_dir.y + ghost_dir.z*ghost_dir.z;
    
    if (in_shadow(start, ghost_dir, triangle_vertexes, ghost_radius_sq, triangle_n)) {
      total_colour -= soft_shadow_color_step;
    }
  }
  
  total_colour += (light_color * max(dot(dir, intersect_normal), 0.0f)) / (4 * ((float)M_PI) * radius_sq);

  
  return total_colour;
}




kernel void draw(global uint  *screen_buffer,    global float3 *triangle_vertexes,   global float3 *triangle_normals,
				 global float3 *triangle_colors, global float3 *rot_matrix,           float3 camera_pos, float3 light_pos, 
				 int triangle_n, float focal_length, local float3 *LOC_triangle_vertexes,  local float3 *LOC_triangle_normals,
				 local float3 *LOC_triangle_colors)
{         /* accumulated magnitudes of velocity for each cell */
  const short x = get_global_id(0);
  const short y = get_global_id(1);


  event_t e = async_work_group_copy(LOC_triangle_vertexes, triangle_vertexes, triangle_n*3, 0);
  e         = async_work_group_copy(LOC_triangle_normals,  triangle_normals,  triangle_n,   0);
  e         = async_work_group_copy(LOC_triangle_colors,   triangle_colors,   triangle_n,   0);
  wait_group_events(3, &e);

  const char rays_x = 2;
  const char rays_y = 2;


  float3 final_color_total = (float3) (0.0f);

  // const float rndx = rnd(x/y);
  // const float rndy = rnd(y/x);

  const float rndx = 0.0f;
  const float rndy = 0.0f;

  for (float dy = y*rays_y; dy < (y+1)*rays_y; dy+=1)  {

    for (float dx = x*rays_y; dx < (x+1)*rays_x; dx+=1)  {
    // Declare ray for given position on the screen. Rotate ray by current view angle
        float3 d = (float3) (dx - (SCREEN_WIDTH*rays_x)/2.0f + rndx, dy - (SCREEN_HEIGHT*rays_y)/2.0f + rndy, focal_length);
        d        = (float3) (dot(rot_matrix[0], d), dot(rot_matrix[1], d), dot(rot_matrix[2], d));

        // Find intersection point with closest geometry. If no intersection, paint the abyss
        Intersection intersect;
        if (closest_intersection(camera_pos, d, LOC_triangle_vertexes, &intersect, triangle_n)) {
          const float3 p = LOC_triangle_colors[intersect.triangle_index];
          const float3 final_color = p*(indirect_light + direct_light(intersect, LOC_triangle_vertexes, LOC_triangle_normals, 
                                                                      light_pos, triangle_n, LOC_triangle_normals[intersect.triangle_index]));
          final_color_total += final_color;
        }
    }
  }
  PutPixelSDL(screen_buffer, (short)x, (short)y, final_color_total/(rays_x*rays_y));
}  


// uint3 getRGB(uint pixel)  {
//   return (uint3) ((uint)((pixel >> 16) & 255), (uint)((pixel >> 8) & 255), (uint)(pixel & 255));
// }

// kernel void average_pixels(global uint *screen_buffer)  {
//   const short x = get_global_id(0);
//   const short y = get_global_id(1);

//   const short nx = get_global_size(1);
  
//   uint3 surrounding_cell_total = (uint3) (0, 0, 0);
//   surrounding_cell_total  += getRGB(screen_buffer[(y*2)*SCREEN_WIDTH+(x*2)]);
//   surrounding_cell_total  += getRGB(screen_buffer[(y*2)*SCREEN_WIDTH+(x*2+1)]);
// // 
//   surrounding_cell_total  += getRGB(screen_buffer[(y*2+1)*SCREEN_WIDTH+(x*2)]);
//   surrounding_cell_total  += getRGB(screen_buffer[(y*2+1)*SCREEN_WIDTH+(x*2+1)]);


//   surrounding_cell_total /= 4;

//   screen_buffer[y*nx+x] = (128<<24) + (surrounding_cell_total.x<<16) + (surrounding_cell_total.y<<8)
//                                                                       + surrounding_cell_total.z;
// }


