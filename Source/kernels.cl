#pragma OPENCL EXTENSION cl_khr_fp64 : enable

constant float3 indirect_light = (float3)(0.35f, 0.35f, 0.35f);
constant float3 light_color    = (float3) (15.0f, 15.0f, 15.0f);
constant float3 bias           = (float3) (0.00001f, 0.00001f, 0.00001f);
constant char rays_x = 3;
constant char rays_y = 3;
constant char aa_rays = 9;

#define SCREEN_WIDTH 1024.0f
#define SCREEN_HEIGHT 1024.0f

/////READ ONLY BUFFERS

typedef struct  {
  float3 position;
  // float distance;
  int triangle_index;
} Intersection;

typedef struct  {
  float3 start;
  float3 direction;
  float3 intersect;
  int intersect_triangle;
} Ray;

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

inline uint3 random(uint3 seed) { // XORSHIFT Alg.
  seed ^= seed << 13;
  seed ^= seed >> 17;
  seed ^= seed << 5;  
  // seed ^= seed << 13; //Not needed as chaining random now
  // seed ^= seed >> 17;
  // seed ^= seed << 5; 
  // seed ^= seed << 13;
  // seed ^= seed >> 17;
  // seed ^= seed << 5; 
  return seed;
  // const float3 fl_seed = range*((float3)seed)/((float3)UINT_MAX);
  // return fl_seed - range/2.f;
}

inline float3 crush(uint3 vector, float range) { // XORSHIFT Alg.

  const float3 fl_seed = range*(convert_float3(vector))/((float3)UINT_MAX);
  return fl_seed - range/2.f;
}

bool closest_intersection(Ray *ray, local float3 *triangle_vertexes, local float3 *triangle_normals, int triangle_n) {
  // Set closest intersection to be the max float value
  float current_t = MAXFLOAT;
  // Make 4D ray into 3D ray
  for (uint i = 0; i < triangle_n; i++) {
    // Define two corners of triangle relative to the other corner
      const float3 v0 = triangle_vertexes[i*3];

      const float3 e1 = triangle_vertexes[i*3+1] - v0;
      const float3 e2 = triangle_vertexes[i*3+2] - v0;
      const float3 b  = ray->start - v0;

      const float3 A[3]  = {-ray->direction, e1, e2};
      const float3 A0[3] = {b,  e1, e2};

      const float detA_recip  = native_recip(det(A));
      const float t           = det(A0)*detA_recip;

      // If ray goes through triangle, and is the closest triangle
      if (  t < current_t  ) {
        const float3 A1[3] = {-ray->direction, b,  e2};
        const float3 A2[3] = {-ray->direction, e1, b};
        const float u = det(A1)*detA_recip;
        const float v = det(A2)*detA_recip; 

        if (u >= 0 && v >= 0 && (u+v) <= 1) {
          // float3 dist_vec                      = t*dir;
          // closest_intersection->distance       = native_sqrt(dist_vec.x*dist_vec.x + dist_vec.y*dist_vec.y + dist_vec.z*dist_vec.z);
          ray->intersect          =  v0 + (u * e1) + (v * e2);
          ray->intersect_triangle = i;
          current_t                            = t;
        }

      }
  }
  if (current_t == MAXFLOAT) return false;
  return true;
}

// bool closest_intersection2(float3 start, float3 dir, local float3 *triangle_vertexes, local float3 *triangle_normals, private Intersection* closest_intersection, int triangle_n) {
//   // Set closest intersection to be the max float value
//   float current_t = MAXFLOAT;
//   // Make 4D ray into 3D ray
//   for (uint i = 0; i < triangle_n; i++) {
//     // Define two corners of triangle relative to the other corner
//     if(dot(dir, triangle_normals[i]) < 0.f)  {
//       const float3 v0 = triangle_vertexes[i*3];

//       const float3 e1 = triangle_vertexes[i*3+1] - v0;
//       const float3 e2 = triangle_vertexes[i*3+2] - v0;
//       const float3 b  = start - v0;

//       const float3 A[3]  = {-dir, e1, e2};
//       const float3 A0[3] = {b,  e1, e2};

//       const float detA_recip  = native_recip(det(A));
//       const float t           = det(A0)*detA_recip;

//       // If ray goes through triangle, and is the closest triangle
//       if (  t < current_t  ) {
//         const float3 A1[3] = {-dir, b,  e2};
//         const float3 A2[3] = {-dir, e1, b};
//         const float u = det(A1)*detA_recip;
//         const float v = det(A2)*detA_recip; 

//         if (u >= 0 && v >= 0 && (u+v) <= 1) {
//           // float3 dist_vec                      = t*dir;
//           // closest_intersection->distance       = native_sqrt(dist_vec.x*dist_vec.x + dist_vec.y*dist_vec.y + dist_vec.z*dist_vec.z);
//           closest_intersection->position       = v0 + (u * e1) + (v * e2);
//           closest_intersection->triangle_index = i;
//           current_t                            = t;
//         }

//       }
//     }
//   }
//   if (current_t == MAXFLOAT) return false;
//   return true;
// }


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

float3 direct_light(const Ray ray, local float3 *triangle_vertexes, local float3 *triangle_normals, float3 light_pos, int triangle_n, const float3 intersect_normal, const int global_id) {

  //Declare colour for point to be 0
  const short light_sources = 10;
  const float light_spread = 0.05f;
  float3 total_light = (float3) 0.0f;
  uint3 rand_vec = random((uint3) (global_id, global_id*91.0f, global_id*19.0f));


  //Get vector from intersection point to light position, and its magnitude
  Ray shadow_ray;
  shadow_ray.direction = light_pos - ray.intersect;
  shadow_ray.start     = ray.intersect + bias*shadow_ray.direction;
  const float radius_sq = shadow_ray.direction.x*shadow_ray.direction.x + shadow_ray.direction.y*shadow_ray.direction.y + shadow_ray.direction.z*shadow_ray.direction.z;


  // Check parallel ghost surfaces for soft triangles
  for (int i = 0; i < light_sources; i++)  {
    rand_vec = random(rand_vec);

    if (!in_shadow(shadow_ray.start, shadow_ray.direction + crush(rand_vec, light_spread), triangle_vertexes, triangle_normals, radius_sq, triangle_n)) {
      total_light += (light_color * max(dot(shadow_ray.direction , intersect_normal), 0.0f)) / ( 4.0f * ((float)M_PI) * radius_sq);
      // total_light += (light_color * max(dot(dir, intersect_normal+5*rand_vec), 0.0f)) / ( 4 * ((float)M_PI) * radius_sq);
    }
  }

  return total_light/light_sources;
}

inline float3 reflect_ray(float3 ray, float3 normal)  {
  return ray-2*(dot(ray, normal)*normal);
}

// float3 secondary_light(float3 start, float3 dir, local float3 *triangle_vertexes, local float3 *triangle_normals, local float4 *triangle_colors, int triangle_n, float3 light_pos, const int global_id, const int depth, float diffusity)  {
  
//   Intersection intersect;
//   float3 light_accumulator_total = (float3) 0.f;
//   float fraction = 0.8;
//   for (int d = 0; d < depth; d++) {
//       if (closest_intersection2(start + bias*dir, dir, triangle_vertexes, triangle_normals, &intersect, triangle_n)) {
//       float3 light_accumulator = indirect_light+direct_light(intersect, triangle_vertexes, triangle_normals, light_pos, triangle_n, triangle_normals[intersect.triangle_index], global_id);
//       light_accumulator *= fraction*(1-diffusity)*triangle_colors[intersect.triangle_index].xyz;
//       light_accumulator_total += light_accumulator;
      
//       dir = reflect_ray(dir, triangle_normals[intersect.triangle_index]);
//       start = intersect.position;
//       diffusity = triangle_colors[intersect.triangle_index].w;
//       fraction *= 0.8f;

//     }
//   }

//   return light_accumulator_total;
// }

float3 norm(float3 vec) {
  float len = native_rsqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
  return (float3) { vec.x * len, vec.y * len, vec.z * len };
}

kernel void draw(global uint  *screen_buffer,    global float3 *triangle_vertexes,   global float3 *triangle_normals,
				 global float4 *triangle_colors, global float3 *rot_matrix,           float3 camera_pos, float3 light_pos, 
				 int triangle_n, float focal_length, local float3 *LOC_triangle_vertexes,  local float3 *LOC_triangle_normals,
				 local float4 *LOC_triangle_colors)
{         /* accumulated magnitudes of velocity for each cell */
  event_t e = async_work_group_copy(LOC_triangle_vertexes, triangle_vertexes, triangle_n*3, 0);
  e         = async_work_group_copy(LOC_triangle_normals,  triangle_normals,  triangle_n,   0);
  e         = async_work_group_copy(LOC_triangle_colors,   triangle_colors,   triangle_n,   0);

  const int x = get_global_id(0);
  const int y = get_global_id(1);
  const int global_id = y*SCREEN_WIDTH+x;
  float3 final_color_total = (float3) (0.0f);

  const float3 base_dir = (float3) (x*rays_x - (SCREEN_WIDTH*rays_x)/2.0f, y*rays_y - (SCREEN_HEIGHT*rays_y)/2.0f, focal_length);



  Ray rays[9];
  const float3 r0 = rot_matrix[0];
  const float3 r1 = rot_matrix[1];
  const float3 r2 = rot_matrix[2];

  for (char dy = 0; dy < rays_y; dy++)  {
    for (char dx = 0; dx < rays_x; dx++)  {
      rays[r].start = camera_pos;

      rays[dy*rays_y + dx].direction = base_dir + (float3) (dx, dy, 0.0f);
      rays[dy*rays_y + dx].direction = (float3) (dot(r0, rays[dy*rays_y + dx].direction), 
                                                 dot(r1, rays[dy*rays_y + dx].direction), 
                                                 dot(r2, rays[dy*rays_y + dx].direction));
    }
  }

  wait_group_events(3, &e);
  for (char r = 0; r < rays_y*rays_x; dy++)  {
    // Find intersection point with closest geometry. If no intersection, paint the abyss
    if (closest_intersection(&rays[r], LOC_triangle_vertexes, LOC_triangle_normals, triangle_n)) {
      const float4 color = LOC_triangle_colors[rays[r].intersect_triangle];
      const float3 first_light = direct_light(rays[r], LOC_triangle_vertexes, LOC_triangle_normals, light_pos, triangle_n, LOC_triangle_normals[rays[r].intersect_triangle], global_id);

      // const float diffusity = color.w; // Diffuse = 1, mirror = 0        
      // const float3 outgoing_dir = reflect_ray(dir, LOC_triangle_normals[intersect.triangle_index]);
      // const float3 second_light = secondary_light(intersect.position, outgoing_dir, LOC_triangle_vertexes, LOC_triangle_normals, LOC_triangle_colors, triangle_n, light_pos, global_id, 1, diffusity);

      final_color_total += color.xyz*(indirect_light+first_light); //Add indirect back! 
    }
  }
  PutPixelSDL(screen_buffer, (short)x, (short)y, final_color_total/(rays_x*rays_y));
}



