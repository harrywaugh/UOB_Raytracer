#pragma OPENCL EXTENSION cl_khr_fp64 : enable

constant float3 indirect_light = (float3)(0.35f, 0.35f, 0.35f);
constant float3 light_color    = (float3) (15.0f, 15.0f, 15.0f);
constant float3 bias           = (float3) (0.0000001f, 0.0000001f, 0.0000001f);
constant char rays_x = 3;
constant char rays_y = 3;
#define aa_rays 9

#define SCREEN_WIDTH 1024.0f
#define SCREEN_HEIGHT 1024.0f
#define GLASS 1.52f
#define AIR 1.00f

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

inline Ray reflect_ray(Ray ray)  {
  Ray reflected_ray;
  reflected_ray.intersect_triangle = -1;
  reflected_ray.intersect_color.w = 1.0f;

  reflected_ray.direction = (ray.direction-2*(dot(ray.direction, ray.intersect_normal)*ray.intersect_normal)); // MINUS here???
  reflected_ray.start = ray.intersect + bias*reflected_ray.direction;
  reflected_ray.medium = AIR;
  return reflected_ray;
}
 



inline Ray refract_ray(Ray ray)  {


  
  Ray refracted_ray;
  float cosi = max(-1.0f, min(1.0f, (dot(ray.intersect_normal, ray.direction))));

  float etai = AIR, etat = GLASS;
  float3 n = ray.intersect_normal;
  if (cosi < 0) { cosi = -cosi; } else { etai = GLASS; etat = AIR; n *= -1.0f; }

  float eta = etai / etat;
  float k = 1.0f - eta*eta*(1.0f - cosi*cosi);


  refracted_ray.direction = (k < 0) ? (float3)(0.0f) : eta*ray.direction + (eta*cosi-sqrtf(k))*n;

  refracted_ray.intersect_triangle = -1;
  refracted_ray.intersect_color.w = 1.0f;
  refracted_ray.start = ray.intersect + bias*refracted_ray.direction;
  refracted_ray.medium = GLASS;

  // printf(" (%f %f %f) (%f %f %f)\n", ray.direction.x, ray.direction.y, ray.direction.z, refracted_ray.direction.x, refracted_ray.direction.y, refracted_ray.direction.z);

  return refracted_ray;
}



void batch_ray_intersections(Ray *rays, local float3 *triangle_vertexes, local float3 *triangle_normals, local float4 *triangle_colors, int triangle_n) {
  // Set closest intersection to be the max float value
  // Make 4D ray into 3D ray

  for (uint r = 1; r < aa_rays-1; r+=2)  {
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
        current_t                  = t;
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
    // total_light += (light_color * max(dot(dir, intersect_normal+5*rand_vec), 0.0f)) / ( 4 * ((float)M_PI) * radius_sq);
  }

  return total_light/light_sources;
}

float3 secondary_light(Ray ray, local float3 *triangle_vertexes, local float3 *triangle_normals, local float4 *triangle_colors, int triangle_n, float3 light_pos, const int global_id)  {
  const int bounces = 5;

  Ray perturbed_ray = ray;
  float3 light_accumulator = (float3) 0.f;

  for (int b = 0; b < bounces && perturbed_ray.intersect_color.w <= 0.0f; b++ )  {
    perturbed_ray = (perturbed_ray.intersect_color.w == 0.0f) ? reflect_ray(perturbed_ray) : refract_ray(perturbed_ray);
    single_ray_intersections(&perturbed_ray, triangle_vertexes, triangle_normals, triangle_colors, triangle_n);



    int intersected = (perturbed_ray.intersect_triangle != -1  && perturbed_ray.intersect_color.w > 0.0f);

    if (intersected)  {
      light_accumulator += indirect_light + direct_light(perturbed_ray, triangle_vertexes, triangle_colors,  light_pos, triangle_n, perturbed_ray.intersect_normal, global_id);
    }

  }
  
  // dir = reflect_ray(dir, triangle_normals[intersect.triangle_index]);
  // start = intersect.position;
  // diffusity = triangle_colors[intersect.triangle_index].w;

  return 0.8f*light_accumulator*perturbed_ray.intersect_color.xyz;
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
    }
  }

  wait_group_events(3, &e);

  batch_ray_intersections(&rays, LOC_triangle_vertexes, LOC_triangle_normals, LOC_triangle_colors, triangle_n);

  for (char r = 1; r < aa_rays-1; r+=2)  {
    // Find intersection point with closest geometry. If no intersection, paint the abyss
    if (rays[r].intersect_triangle != -1) {
      // const float diffusity = color.w; // Diffuse = 1, mirror = 0        
      float3 second_light = (float3) 0.0f;
      if(rays[r].intersect_color.w <= 0.0f)  { // Mirror or glass

        final_color_total += secondary_light(rays[r], LOC_triangle_vertexes, LOC_triangle_normals, LOC_triangle_colors, triangle_n, light_pos, global_id);
      } else { 
        const float3 first_light  = direct_light(rays[r], LOC_triangle_vertexes, LOC_triangle_colors,  light_pos, triangle_n, rays[r].intersect_normal, global_id);
        final_color_total += rays[r].intersect_color.xyz*(indirect_light + first_light);
      }
    }
  }
  PutPixelSDL(screen_buffer, (short)x, (short)y, final_color_total/4);
}



