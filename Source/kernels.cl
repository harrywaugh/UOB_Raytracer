#pragma OPENCL EXTENSION cl_khr_fp64 : enable

constant float3 indirect_light = (float3)(0.35f, 0.35f, 0.35f);
constant float3 light_color    = (float3) (15.0f, 15.0f, 15.0f);
constant float3 bias           = (float3) (0.00001f, 0.00001f, 0.00001f);
#define SCREEN_WIDTH 1024.0f
#define SCREEN_HEIGHT 1024.0f

/////READ ONLY BUFFERS

typedef struct  {
  float3 position;
  // float distance;
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
  screen_buffer[y*(short)SCREEN_WIDTH+x] = (rgb.x<<16) + (rgb.y<<8) + rgb.z;
}
inline float rnd(uint seed, float range) { // XORSHIFT Alg.
  seed ^= seed << 13;
  seed ^= seed >> 17;
  seed ^= seed << 5;  
  const float fl_seed = range*((float)seed)/((float)UINT_MAX);
  return fl_seed - range/2.f;
}

bool closest_intersection(float3 start, float3 dir, local float3 *triangle_vertexes, local float3 *triangle_normals, private Intersection* closest_intersection, int triangle_n) {
  // Set closest intersection to be the max float value
  float current_t = MAXFLOAT;
  // Make 4D ray into 3D ray
  for (uint i = 0; i < triangle_n; i++) {
    // Define two corners of triangle relative to the other corner
    if(dot(dir, triangle_normals[i]) < 0.f)  {
      const float3 v0 = triangle_vertexes[i*3];

      const float3 e1 = triangle_vertexes[i*3+1] - v0;
      const float3 e2 = triangle_vertexes[i*3+2] - v0;
      const float3 b  = start - v0;

      // Cramers, might be det repeated computation..?
      const float3 A[3]  = {-dir, e1, e2};
      const float3 A0[3] = {b,  e1, e2};

      const float detA_recip  = native_recip(det(A));
      const float t = det(A0)*detA_recip;

      // If ray goes through triangle, and is the closest triangle
      if ( t >= 0 && t < current_t  ) {
        const float3 A1[3] = {-dir, b,  e2};
        const float3 A2[3] = {-dir, e1, b};
        const float u = det(A1)*detA_recip;
        const float v = det(A2)*detA_recip; 

        if (u >= 0 && v >= 0 && (u+v) <= 1) {
          // float3 dist_vec                      = t*dir;
          // closest_intersection->distance       = native_sqrt(dist_vec.x*dist_vec.x + dist_vec.y*dist_vec.y + dist_vec.z*dist_vec.z);
          closest_intersection->position       = v0 + (u * e1) + (v * e2);
          closest_intersection->triangle_index = i;
          current_t                            = t;
        }

      }
    }
    
  }
  if (current_t == MAXFLOAT) return false;
  return true;
}

bool in_shadow(float3 start, float3 dir, local float3 *triangle_vertexes, local float3 *triangle_normals,  float radius_sq, int triangle_n) {
  // Make 4D ray into 3D ray

  for (uint i = 0; i < triangle_n; i++) {
    // Define two corners of triangle relative to the other corner
    if(dot(dir, triangle_normals[i]) < 0.f)  {
      const float3 v0 = triangle_vertexes[i*3];

      const float3 e1 = triangle_vertexes[i*3+1] - v0;
      const float3 e2 = triangle_vertexes[i*3+2] - v0;
      const float3 b  = start - v0;

      const float3 A[3]  = {-dir, e1, e2};
      const float3 A0[3] = {b,  e1, e2};
      const float detA_recip  = native_recip(det(A));
      const float t = det(A0)*detA_recip;

      const float3 dist_vec            = t*dir;
      const float intersect_dist       = dist_vec.x*dist_vec.x + dist_vec.y*dist_vec.y + dist_vec.z*dist_vec.z;

      if (t >= 0 && intersect_dist < radius_sq ) {
        const float3 A1[3] = {-dir, b,  e2};
        const float3 A2[3] = {-dir, e1, b};
        const float u = det(A1)*detA_recip;
        const float v = det(A2)*detA_recip; 

        if (u >= 0 && v >= 0 && (u+v) <= 1 ) {
          return true;
        }
      }
    }
  }
  return false;
}

float3 direct_light(const Intersection intersection, local float3 *triangle_vertexes, local float3 *triangle_normals, float3 light_pos, int triangle_n, const float3 intersect_normal, const int global_id) {

  //Declare colour for point to be 0
  float3 total_colour = (float3) 0.0f;

  //Get vector from intersection point to light position, and its magnitude
  float3 dir = light_pos - intersection.position;
  float radius_sq = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;

  //Declare threshold to get intersection position that is not going to intersect with own triangle
  float3 start = intersection.position + bias*dir;

  const short light_sources = 20;
  short light_count = light_sources;

  const float light_spread = 0.01f;

  // Check parallel ghost surfaces for soft triangles
  for (int i = 0; i < light_sources; i++)  {
    float3 ghost_dir = dir + (float3) (rnd(i*(global_id), light_spread), rnd(i*(global_id), light_spread), rnd(i*(global_id), light_spread));
    float ghost_radius_sq = ghost_dir.x*ghost_dir.x + ghost_dir.y*ghost_dir.y + ghost_dir.z*ghost_dir.z;

    if (!in_shadow(start, ghost_dir, triangle_vertexes, triangle_normals, ghost_radius_sq, triangle_n)) {
      total_colour += (light_color * max(dot(dir, intersect_normal), 0.0f)) / (light_sources * 4 * ((float)M_PI) * radius_sq);
    }
  }

  return total_colour;
}

inline float3 reflect_ray(float3 ray, float3 normal)  {
  return ray-2*(dot(ray, normal)*normal);
}

float3 secondary_light(float3 start, float3 dir, local float3 *triangle_vertexes, local float3 *triangle_normals, local float3 *triangle_colors, int triangle_n, float3 light_pos, const int global_id)  {
  
  Intersection intersect;
  float3 light_accumulator = (float3) 0.f;

  if (closest_intersection(start + bias*dir, dir, triangle_vertexes, triangle_normals, &intersect, triangle_n)) {
    light_accumulator = direct_light(intersect, triangle_vertexes, triangle_normals, light_pos, triangle_n, triangle_normals[intersect.triangle_index], global_id);
    light_accumulator *= triangle_colors[intersect.triangle_index];
  }
  return light_accumulator;
}

kernel void draw(global uint  *screen_buffer,    global float3 *triangle_vertexes,   global float3 *triangle_normals,
				 global float3 *triangle_colors, global float3 *rot_matrix,           float3 camera_pos, float3 light_pos, 
				 int triangle_n, float focal_length, local float3 *LOC_triangle_vertexes,  local float3 *LOC_triangle_normals,
				 local float3 *LOC_triangle_colors)
{         /* accumulated magnitudes of velocity for each cell */
  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const int global_id = y*SCREEN_WIDTH+x;

  event_t e = async_work_group_copy(LOC_triangle_vertexes, triangle_vertexes, triangle_n*3, 0);
  e         = async_work_group_copy(LOC_triangle_normals,  triangle_normals,  triangle_n,   0);
  e         = async_work_group_copy(LOC_triangle_colors,   triangle_colors,   triangle_n,   0);
  wait_group_events(3, &e);

  const char rays_x = 3;
  const char rays_y = 3;

  float3 final_color_total = (float3) (0.0f);

  for (float dy = y*rays_y; dy < (y+1)*rays_y; dy+=1)  {

    for (float dx = x*rays_y; dx < (x+1)*rays_x; dx+=1)  {
    // Declare ray for given position on the screen. Rotate ray by current view angle
        float3 dir = (float3) (dx - (SCREEN_WIDTH*rays_x)/2., dy - (SCREEN_HEIGHT*rays_y)/2., focal_length);
        dir        = (float3) (dot(rot_matrix[0], dir), dot(rot_matrix[1], dir), dot(rot_matrix[2], dir));
        


        // Find intersection point with closest geometry. If no intersection, paint the abyss
        Intersection intersect;
        if (closest_intersection(camera_pos, dir, LOC_triangle_vertexes, LOC_triangle_normals, &intersect, triangle_n)) {
          const float3 p = LOC_triangle_colors[intersect.triangle_index];
          if (p.x == -1.0f)  {
            const float3 outgoing_dir    = reflect_ray(dir, LOC_triangle_normals[intersect.triangle_index]);
            final_color_total += 0.8f*secondary_light(intersect.position, outgoing_dir, LOC_triangle_vertexes, LOC_triangle_normals, LOC_triangle_colors, triangle_n, light_pos, global_id);

          }  else {
            const float3 light = direct_light(intersect, LOC_triangle_vertexes, LOC_triangle_normals, light_pos, triangle_n, LOC_triangle_normals[intersect.triangle_index], global_id);        
            final_color_total += p*(indirect_light+light); //Add indirect back!
            
          }
        }
    }
  }
  PutPixelSDL(screen_buffer, (short)x, (short)y, final_color_total/(rays_x*rays_y));
}



