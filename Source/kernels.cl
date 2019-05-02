#pragma OPENCL EXTENSION cl_khr_fp64 : enable

constant float3 indirect_light = (float3)(0.5f, 0.5f, 0.5f);
constant float3 light_color    = (float3) (16.0f, 16.0f, 16.0f);
constant float3 bias           = (float3) (0.0001f, 0.0001f, 0.0001f);
// constant float4 ball_color     = (float4) (0.5f, 0.0f, 0.0f, -1.0f);
// constant float3 circle_center  = (float3) (0.0f, 0.6f, 0.0f);
constant float4 sphere_centers[2]   = {(float4) (-0.5f, 0.7f, 0.0f, 0.0f), (float4) (0.5f, 0.7f, 0.0f, 0.0f)};
constant float sphere_radius_sqs[2] = {0.1f, 0.1f};
// constant float circle_radius_sq = 0.1f;
constant float4 sphere_colors[2]     = {(float4) (0.0f, 0.0f, 0.0f, 0.0f), (float4) (0.0f, 0.0f, 0.0f, -1.0f)};

constant char rays_x = 3;
constant char rays_y = 3;
#define aa_rays 9

#define SCREEN_WIDTH 1024.0f
#define SCREEN_HEIGHT 1024.0f
#define GLASS 2.0f
#define AIR 1.0f
#define SPHERES 2

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
  float3 intersect_normal;
  float4 intersect_color;
  float medium;
  int intersect_triangle;
} Ray;

inline float det(const float3 *M) {
	return M[0].x * (M[1].y * M[2].z - M[1].z * M[2].y) -
		     M[0].y * (M[1].x * M[2].z - M[1].z * M[2].x) +
		     M[0].z * (M[1].x * M[2].y - M[1].y * M[2].x);
}

inline void color_pixel(global uint *screen_buffer, int x, int y, float3 colour) {
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

inline Ray reflect_ray(Ray ray)  {
  Ray reflected_ray;
  reflected_ray.intersect_triangle = -1;
  reflected_ray.intersect_color.w = 1.0f;

  reflected_ray.direction = (ray.direction-2*(dot(ray.direction, ray.intersect_normal)*ray.intersect_normal)); // MINUS here???
  reflected_ray.start     = ray.intersect + bias*reflected_ray.direction;
  reflected_ray.medium    = AIR;
  reflected_ray.direction = normalize(reflected_ray.direction);

  return reflected_ray;
}

inline Ray refract_ray(Ray ray)  {
  Ray refracted_ray;
  float3 normal          = ray.intersect_normal;
  const char medium_mask = (ray.medium == AIR); 
  const float n1         = medium_mask*AIR   + (medium_mask^1)*GLASS; 
  const float n2         = medium_mask*GLASS + (medium_mask^1)*AIR; 
  float c1               = (dot(normal, ray.direction));
  if (c1 < 0.0f) normal *= -1.0f;   
  c1 = fabs(c1);
  const float n          = native_divide(n1, n2);
  const float c2         = native_sqrt(1-(n*n)*(1-(c1*c1)));
  if(c2 < 0.0f)  {
    return reflect_ray(ray);
  }  
  refracted_ray.intersect_triangle = -1;
  refracted_ray.intersect_color = (float4) (1.0f, 0.0f, 0.0f, 1.0f);
  refracted_ray.direction = (n*ray.direction + (n*c1-c2)*(-normal));
  refracted_ray.start = ray.intersect + bias*refracted_ray.direction;
  refracted_ray.medium = n2;
  refracted_ray.direction = normalize(refracted_ray.direction);
  return refracted_ray;
}



void batch_ray_intersections(Ray *rays, local float3 *triangle_vertexes, local float3 *triangle_normals, local float4 *triangle_colors, int triangle_n) {
  // Set closest intersection to be the max float value
  // Make 4D ray into 3D ray

  for (uint r = 0; r < aa_rays; r++)  {
    float current_t = MAXFLOAT;
    

    for (uint i = 0; i < triangle_n; i++) {
        // Define two corners of triangle relative to the other corner
      const float3 v0 = triangle_vertexes[i*3];
      const float3 e1 = triangle_vertexes[i*3+1] - v0;
      const float3 e2 = triangle_vertexes[i*3+2] - v0;


      const float3 b  = rays[r].start - v0;


      const float3 A[3]  = {-rays[r].direction, e1, e2};
      const float3 A0[3] = {b,  e1, e2};


      const float detA_recip  = native_recip(det(A));
      const float t           = det(A0)*detA_recip;


      const float3 A1[3] = {-rays[r].direction, b,  e2};
      const float3 A2[3] = {-rays[r].direction, e1, b};
      const float u = det(A1)*detA_recip;
      const float v = det(A2)*detA_recip; 


      const int intersected = (t < current_t && u >= 0 && v >= 0 && (u+v) <= 1 && t >=0);
      // const int not_intersected = intersected ^ 1;

      // If ray goes through triangle, and is the closest triangle
      if (  intersected  ) {
        // ray->intersect          = intersected*(v0 + (u*e1) + (v*e2)) + not_intersected*ray->intersect;
        // current_t[r]               = intersected*t + not_intersected*current_t[r] ;
        // ray->intersect_triangle = intersected*i + not_intersected*ray->intersect_triangle;
        rays[r].intersect_triangle = i;
        rays[r].intersect          = (v0 + (u*e1) + (v*e2));
        rays[r].intersect_normal   = triangle_normals[i];
        rays[r].intersect_color    = triangle_colors[i];
        // printf("box %f %f %f\n", rays[r].intersect_color.x, rays[r].intersect_color.y, rays[r].intersect_color.z);
        current_t                  = t;
      }
    }
    for (int i=0; i < SPHERES; i++)  {
      
      

      const float3 L = (rays[r].start - sphere_centers[i].xyz);
      const float a = dot(rays[r].direction, rays[r].direction);
      const float b = 2*dot(rays[r].direction, L);
      const float c = dot(L, L) - sphere_radius_sqs[i];

      const float discriminant = b*b - 4.0f*a*c;
      if (discriminant < 0.0f) continue;
      const float q = (b > 0) ? -0.5 * (b + sqrt(discriminant)) : -0.5 * (b - sqrt(discriminant)); 
      const float x0 = q / a;
      const float x1 = c / q;

      const float x_min = min(x0, x1);
      const float x_max = max(x0, x1);

      if(x_min >= 0.0f && x_min < current_t)  {
        rays[r].intersect_triangle = -2;
        rays[r].intersect          = rays[r].start + rays[r].direction*x_min;
        rays[r].intersect_normal   = normalize(rays[r].intersect - sphere_centers[i].xyz);
        rays[r].intersect_color    = sphere_colors[i];

        current_t                  = x_min;
      } else if(x_max >= 0.0f && x_max < current_t)  {
        rays[r].intersect_triangle = -2;
        rays[r].intersect          = rays[r].start + rays[r].direction*x_max;
        rays[r].intersect_normal   = normalize(rays[r].intersect - sphere_centers[i].xyz);
        rays[r].intersect_color    = sphere_colors[i];

        current_t                  = x_max;
      }
    }
    
  }
}

void single_ray_intersections(Ray *ray, local float3 *triangle_vertexes, local float3 *triangle_normals, local float4 *triangle_colors, int triangle_n) {
  // Set closest intersection to be the max float value
  // Make 4D ray into 3D ray

  float current_t = MAXFLOAT;

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


    const float3 A1[3] = {-ray->direction, b,  e2};
    const float3 A2[3] = {-ray->direction, e1, b};
    const float u = det(A1)*detA_recip;
    const float v = det(A2)*detA_recip; 


    const int intersected = (t < current_t && u >= 0 && v >= 0 && (u+v) <= 1 && t >=0);
    // const int not_intersected = intersected ^ 1;

    // If ray goes through triangle, and is the closest triangle
    if (  intersected  ) {
      // ray->intersect          = intersected*(v0 + (u*e1) + (v*e2)) + not_intersected*ray->intersect;
      // current_t[r]               = intersected*t + not_intersected*current_t[r] ;
      // ray->intersect_triangle = intersected*i + not_intersected*ray->intersect_triangle;
      ray->intersect_triangle = i;
      ray->intersect          = (v0 + (u*e1) + (v*e2));
      ray->intersect_normal   = triangle_normals[i];
      ray->intersect_color    = triangle_colors[i];
      current_t               = t;
    }
  }
  for (int i=0; i < SPHERES; i++)  {
    const float3 L = (ray->start - sphere_centers[i].xyz);
    const float a = dot(ray->direction, ray->direction);
    const float b = 2*dot(ray->direction, L);
    const float c = dot(L, L) - sphere_radius_sqs[i];

    const float discriminant = b*b - 4.0f*a*c;
    if (discriminant < 0.0f) continue;
    const float q = (b > 0) ? -0.5 * (b + sqrt(discriminant)) : -0.5 * (b - sqrt(discriminant)); 
    const float x0 = q / a;
    const float x1 = c / q;

    const float x_min = min(x0, x1);
    const float x_max = max(x0, x1);

    if(x_min >= 0.0f && x_min < current_t)  {
      ray->intersect_triangle = -2;
      ray->intersect          = ray->start + ray->direction*x_min;
      ray->intersect_normal   = normalize(ray->intersect - sphere_centers[i].xyz);
      ray->intersect_color    = sphere_colors[i];
      // printf("ball %f %f %f\n", ray->intersect_color.x, ray->intersect_color.y, ray->intersect_color.z);

      current_t                  = x_min;
    } else if(x_max >= 0.0f && x_max < current_t)  {
      ray->intersect_triangle = -2;
      ray->intersect          = ray->start + ray->direction*x_max;
      ray->intersect_normal   = normalize(ray->intersect - sphere_centers[i].xyz);
      ray->intersect_color    = sphere_colors[i];
      // printf("ball %f %f %f\n", rays[r].intersect_color.x, rays[r].intersect_color.y, rays[r].intersect_color.z);

      current_t                  = x_max;
    }
  }
  
}


bool in_shadow(float3 start, float3 dir, local float3 *triangle_vertexes, local float4 *triangle_colors, float radius_sq, int triangle_n) {
  // Make 4D ray into 3D ray

  for (uint i = 0; i < triangle_n; i++) {
    if(triangle_colors[i].w == -1.0f)
      i+=10;
    // Define two corners of triangle relative to the other corner
    const float3 v0 = triangle_vertexes[i*3];

    const float3 e1 = triangle_vertexes[i*3+1] - v0;
    const float3 e2 = triangle_vertexes[i*3+2] - v0;
    const float3 b  = start - v0;

    const float3 A[3]  = {-dir, e1, e2};
    const float3 A0[3] = {b,  e1, e2};
    const float detA_recip  = native_recip(det(A));
    const float t = det(A0)*detA_recip;

    const float3 dist_vec            = t*dir;
    const float intersect_dist       = dist_vec.x*dist_vec.x + dist_vec.y*dist_vec.y + dist_vec.z*dist_vec.z; //// May be just t??

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
  for (int i=0; i<SPHERES; i++)  {
    if(sphere_colors[i].w == -1.0f) continue;
    const float3 L = (start - sphere_centers[i].xyz);
    const float a = dot(dir, dir);
    const float b = 2*dot(dir, L);
    const float c = dot(L, L) - sphere_radius_sqs[i];

    const float discriminant = b*b - 4.0f*a*c;
    if(discriminant < 0.0f)  return false;
    const float q = (b > 0) ? -0.5 * (b + sqrt(discriminant)) : -0.5 * (b - sqrt(discriminant)); 
    const float x0 = q / a;
    const float x1 = c / q;

    const float x_min = min(x0, x1);
    const float x_max = max(x0, x1);

    const float3 min_dir = x_min*dir;
    const float3 max_dir = x_max*dir;
    const float min_dist = dot(min_dir, min_dir);
    const float max_dist = dot(max_dir, max_dir);

    if(x_min >= 0.0f && min_dist < radius_sq)  {
      return true;
    } else if(x_max >= 0.0f && max_dist < radius_sq)  {
      return true;
    }
  }
  

  return false;
}

float3 direct_light(const Ray ray, local float3 *triangle_vertexes, local float4 *triangle_colors, float3 light_pos, int triangle_n, const float3 intersect_normal, const int global_id) {

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

    int mask = (!in_shadow(shadow_ray.start, shadow_ray.direction + crush(rand_vec, light_spread), triangle_vertexes, triangle_colors, radius_sq, triangle_n));
    total_light += mask*(light_color * max(dot(shadow_ray.direction , intersect_normal), 0.0f)) / ( 4.0f * ((float)M_PI) * radius_sq);
    // total_light += mask*(light_color * max(dot(shadow_ray.direction , intersect_normal+5*crush(rand_vec, 2*light_spread)), 0.0f)) / ( 4.0f * ((float)M_PI) * radius_sq);
  }

  return total_light/light_sources;
}

float3 secondary_light(Ray ray, local float3 *triangle_vertexes, local float3 *triangle_normals, local float4 *triangle_colors, int triangle_n, float3 light_pos, const int global_id)  {
  const int bounces = 10;

  Ray perturbed_ray = ray;

  for (int b = 0; b < bounces; b++ )  {
    perturbed_ray = (perturbed_ray.intersect_color.w == 0.0f) ? reflect_ray(perturbed_ray) : refract_ray(perturbed_ray);
    single_ray_intersections(&perturbed_ray, triangle_vertexes, triangle_normals, triangle_colors, triangle_n);

    int intersected = (perturbed_ray.intersect_triangle != -1  && perturbed_ray.intersect_color.w > 0.0f);
    if (intersected)  {
      const float3 light = indirect_light + direct_light(perturbed_ray, triangle_vertexes, triangle_colors,  light_pos, triangle_n, perturbed_ray.intersect_normal, global_id);
      return 0.75f*light*perturbed_ray.intersect_color.xyz;
    }
  }
  // if(light_accumulator.x == 0.0f)  {
  //   printf("(%f %f %f %f)\n", perturbed_ray.intersect_color.x, perturbed_ray.intersect_color.y, perturbed_ray.intersect_color.z, perturbed_ray.intersect_color.w);

  // }
  return (float3) 0.0f;
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



  Ray rays[aa_rays];
  const float3 r0 = rot_matrix[0];
  const float3 r1 = rot_matrix[1];
  const float3 r2 = rot_matrix[2];

  for (char dy = 0; dy < rays_y; dy++)  {
    for (char dx = 0; dx < rays_x; dx++)  {
      rays[dy*rays_y + dx].start = camera_pos;

      rays[dy*rays_y + dx].direction = base_dir + (float3) (dx, dy, 0.0f);
      rays[dy*rays_y + dx].direction = (float3) (dot(r0, rays[dy*rays_y + dx].direction), 
                                                 dot(r1, rays[dy*rays_y + dx].direction), 
                                                 dot(r2, rays[dy*rays_y + dx].direction));
      rays[dy*rays_y + dx].intersect_triangle = -1;
      rays[dy*rays_y + dx].intersect = (float3) 0.0f;
      rays[dy*rays_y + dx].medium = AIR;
      rays[dy*rays_y + dx].intersect_color = (float4) (0.0f, 0.0f, 0.0f, 1.0f);
      rays[dy*rays_y + dx].direction = normalize(rays[dy*rays_y + dx].direction);
    }
  }

  wait_group_events(3, &e);

  batch_ray_intersections(&rays, LOC_triangle_vertexes, LOC_triangle_normals, LOC_triangle_colors, triangle_n);

  for (char r = 0; r < 9; r++)  {
    // Find intersection point with closest geometry. If no intersection, paint the abyss
    if (rays[r].intersect_triangle != -1) {

      if(rays[r].intersect_color.w <= 0.0f)  { // Mirror or glass
        final_color_total += secondary_light(rays[r], LOC_triangle_vertexes, LOC_triangle_normals, LOC_triangle_colors, triangle_n, light_pos, global_id);
      } else { 
        const float3 first_light  = direct_light(rays[r], LOC_triangle_vertexes, LOC_triangle_colors,  light_pos, triangle_n, rays[r].intersect_normal, global_id);
        final_color_total += rays[r].intersect_color.xyz*(indirect_light + first_light);
      }
    }
  }
  color_pixel(screen_buffer, (short)x, (short)y, final_color_total/((float)aa_rays));
}



