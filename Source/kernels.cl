#pragma OPENCL EXTENSION cl_khr_fp64 : enable

constant float3 indirect_light = (float3)(0.5f, 0.5f, 0.5f); 
constant float3 light_color    = (float3) (16.0f, 16.0f, 16.0f);
constant float3 bias           = (float3) (0.0001f, 0.0001f, 0.0001f);

#define SPHERES 2
constant float4 sphere_centers[SPHERES]   = {(float4) (0.3f, 0.1f, -0.5f, 0.0f), (float4) (-0.4f, 0.8f, -0.5f, 0.0f) ,(float4) (0.0f, 0.0f, -0.8f, 0.0f) };
constant float4 sphere_colors[SPHERES]    = {(float4) (0.0f, 0.f, 0.f, -1.0f),(float4) (0.0f, 0.f, 0.f, 0.0f), (float4) (0.6f, 0.0f, 0.0f, -1.0f) };
constant float sphere_radius_sqs[SPHERES] = {0.075f, 0.05f, 0.1f};

constant char rays_x = 2;
constant char rays_y = 2;
#define aa_rays 4

#define SCREEN_WIDTH 1024.0f
#define SCREEN_HEIGHT 1024.0f
#define GLASS 1.52f
#define AIR 1.0f

typedef struct  {
  float3 start;              // Origin
  float3 direction;          // Direction
  float3 intersect;          // Vector of intersect position
  float3 intersect_normal;   // Normal of intersected surface
  float4 intersect_color;    // Color of intersected surface
  float medium;              // Current refractive index, air or glass
  int intersect_triangle;    // ID of triangle that ray has intersected, -1 nothing, -2 = sphere
} Ray;

inline float det(const float3 *M) {
	return M[0].x * (M[1].y * M[2].z - M[1].z * M[2].y) -
		     M[0].y * (M[1].x * M[2].z - M[1].z * M[2].x) +
		     M[0].z * (M[1].x * M[2].y - M[1].y * M[2].x);
}

inline void color_pixel(global uint *screen_buffer, int x, int y, float3 colour) {
  uint3 rgb = convert_uint3(min(max(255*colour, 0.f), 255.f));
  screen_buffer[y*(short)SCREEN_WIDTH+x] = (255<<24) + (rgb.x<<16) + (rgb.y<<8) + rgb.z;
}

inline uint3 random(uint3 seed) { // XORSHIFT Alg. <- Wikipedia
  seed ^= seed << 13;
  seed ^= seed >> 17;
  seed ^= seed << 5;
  return seed;
}

inline float3 crush(uint3 vector, float range) { // Crushes random uint3 to be float3 within (-range/2 < random_vector < range/2)
  const float3 fl_seed = range*(convert_float3(vector))/((float3)UINT_MAX);
  return fl_seed - range/2.f;
}

inline Ray reflect_ray(Ray ray)  { // Takes ray and computes reflected ray, direction is normalised
  Ray reflected_ray;
  reflected_ray.intersect_triangle = -1;
  reflected_ray.intersect_color.w = 1.0f;

  reflected_ray.direction = (ray.direction-2*(dot(ray.direction, ray.intersect_normal)*ray.intersect_normal)); // MINUS here???
  reflected_ray.start     = ray.intersect + bias*reflected_ray.direction;
  reflected_ray.medium    = AIR;
  reflected_ray.direction = normalize(reflected_ray.direction);

  return reflected_ray;
}

inline Ray refract_ray(Ray ray)  { // Takes ray and computes refracted ray, direction is normalised
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
    return reflect_ray(ray);  // Total internal reflection
  }  
  refracted_ray.intersect_triangle = -1;
  refracted_ray.intersect_color = (float4) (1.0f, 0.0f, 0.0f, 1.0f);
  refracted_ray.direction = (n*ray.direction + (n*c1-c2)*(-normal));
  refracted_ray.start = ray.intersect + bias*refracted_ray.direction;
  refracted_ray.medium = n2;
  refracted_ray.direction = normalize(refracted_ray.direction);
  return refracted_ray;
}

void single_ray_intersections(Ray *ray, local float3 *triangle_vertexes, local float3 *triangle_normals, local float4 *triangle_colors, int triangle_n) {
  // Set closest intersection to be the max float value
  // Make 4D ray into 3D ray

  float current_t = MAXFLOAT;
  // Check triangle intersection
  for (uint i = 0; i < triangle_n; i++) {
    // Define two corners of triangle relative to the other corner
    const float3 v0 = triangle_vertexes[i*3];
    const float3 e1 = triangle_vertexes[i*3+1] - v0;
    const float3 e2 = triangle_vertexes[i*3+2] - v0;

    // Vector from start to triangle
    const float3 b  = ray->start - v0;

    // Cramers rule for inverse
    const float3 A[3]  = {-ray->direction, e1, e2};
    const float3 A0[3] = {b,  e1, e2};
    const float detA_recip  = native_recip(det(A));
    const float t           = det(A0)*detA_recip;
    const float3 A1[3] = {-ray->direction, b,  e2};
    const float3 A2[3] = {-ray->direction, e1, b};
    const float u = det(A1)*detA_recip;
    const float v = det(A2)*detA_recip; 

    // If object closer than last found object, and passes through triangle
    const int intersected = (t < current_t && u >= 0 && v >= 0 && (u+v) <= 1 && t >=0);

    if (  intersected  ) {
      // Store collision information in ray
      ray->intersect_triangle = i;
      ray->intersect          = (v0 + (u*e1) + (v*e2));
      ray->intersect_normal   = triangle_normals[i];
      ray->intersect_color    = triangle_colors[i];
      current_t               = t;
    }
  }


  // Check circle intersections
  for (int i=0; i < SPHERES; i++)  {
    // Form quadratic equation, to check if ray intersects with the sphere
    const float3 L = (ray->start - sphere_centers[i].xyz);
    const float a = dot(ray->direction, ray->direction);
    const float b = 2*dot(ray->direction, L);
    const float c = dot(L, L) - sphere_radius_sqs[i];

    const float discriminant = b*b - 4.0f*a*c;
    if (discriminant < 0.0f) continue; // No solution to quadratic, so check next sphere// No solution to quadratic, so check next sphere
    const float q = (b > 0) ? -0.5 * (b + sqrt(discriminant)) : -0.5 * (b - sqrt(discriminant)); 
    const float x0 = q / a;
    const float x1 = c / q;
    // As per scratchapixel, this solution avoids catastrophic cancellation with similar terms
    const float x_min = min(x0, x1);
    const float x_max = max(x0, x1);
    // If intersection is in front of camera and closest object
    if(x_min >= 0.0f && x_min < current_t)  {
      // Store collision information in ray
      ray->intersect_triangle = -2;
      ray->intersect          = ray->start + ray->direction*x_min;
      ray->intersect_normal   = normalize(ray->intersect - sphere_centers[i].xyz);
      ray->intersect_color    = sphere_colors[i];
      current_t                  = x_min;
    } else if(x_max >= 0.0f && x_max < current_t)  {
      // Store collision information in ray
      ray->intersect_triangle = -2;
      ray->intersect          = ray->start + ray->direction*x_max;
      ray->intersect_normal   = normalize(ray->intersect - sphere_centers[i].xyz);
      ray->intersect_color    = sphere_colors[i];
      current_t                  = x_max;
    }
  }
  
}

bool in_shadow(float3 start, float3 dir, local float3 *triangle_vertexes, local float4 *triangle_colors, float radius_sq, int triangle_n) {
  // Make 4D ray into 3D ray
  // Check triangle intersections
  for (uint i = 0; i < triangle_n; i++) {
    if(triangle_colors[i].w == -1.0f) continue; // No glass shadows!

    // Define two corners of triangle relative to the other corner
    const float3 v0 = triangle_vertexes[i*3];
    const float3 e1 = triangle_vertexes[i*3+1] - v0;
    const float3 e2 = triangle_vertexes[i*3+2] - v0;

    // Vector from start to triangle
    const float3 b  = start - v0;
    //CRAMERS
    const float3 A[3]  = {-dir, e1, e2};
    const float3 A0[3] = {b,  e1, e2};
    const float detA_recip  = native_recip(det(A));
    const float t = det(A0)*detA_recip;

    // Check that object is closer than light source
    const float3 dist_vec            = t*dir;
    const float intersect_dist       = dist_vec.x*dist_vec.x + dist_vec.y*dist_vec.y + dist_vec.z*dist_vec.z;

    if (t >= 0 && intersect_dist < radius_sq ) {
      const float3 A1[3] = {-dir, b,  e2};
      const float3 A2[3] = {-dir, e1, b};
      const float u = det(A1)*detA_recip;
      const float v = det(A2)*detA_recip; 

      if (u >= 0 && v >= 0 && (u+v) <= 1 ) {
        return true; // In shadow, no need to check rest
      }
    }
  }
  // Check circle intersections
  for (int i=0; i<SPHERES; i++)  {
    if(sphere_colors[i].w == -1.0f) continue; // Glass, thus no shadow

    // Form quadratic to solve solutions to sphere intersections
    const float3 L = (start - sphere_centers[i].xyz);
    const float a = dot(dir, dir);
    const float b = 2*dot(dir, L);
    const float c = dot(L, L) - sphere_radius_sqs[i];

    const float discriminant = b*b - 4.0f*a*c;
    if(discriminant < 0.0f )  continue;  // No intersection
    const float q = (b > 0) ? -0.5 * (b + sqrt(discriminant)) : -0.5 * (b - sqrt(discriminant)); 
    const float x0 = q / a;
    const float x1 = c / q;

    const float x_min = min(x0, x1);
    const float x_max = max(x0, x1);

    const float3 min_dir = x_min*dir;
    const float3 max_dir = x_max*dir;
    const float min_dist = dot(min_dir, min_dir);
    const float max_dist = dot(max_dir, max_dir);

    // Light intersections and closer than the light?
    if(x_min >= 0.0f && min_dist < radius_sq)  {
      return true;  // In shadow, no need to check rest
    } else if(x_max >= 0.0f && max_dist < radius_sq)  {
      return true;  // In shadow, no need to check rest
    }
  }
  

  return false;
}

float3 direct_light(const Ray ray, local float3 *triangle_vertexes, local float4 *triangle_colors, float3 light_pos, int triangle_n, const float3 intersect_normal, const int global_id) {

  // Number of fake light sources to generate some cuddly shadows, spread is max seperation of light sources 
  const short light_sources = 10;
  const float light_spread = 0.05f;
  float3 total_light = (float3) 0.0f;
  uint3 rand_vec = random((uint3) (global_id, global_id*91.0f, global_id*19.0f)); // Random position


  // Get shadow ray to fake light source
  Ray shadow_ray;
  shadow_ray.direction = light_pos - ray.intersect;  
  shadow_ray.start     = ray.intersect + bias*shadow_ray.direction;
  const float radius_sq = shadow_ray.direction.x*shadow_ray.direction.x + shadow_ray.direction.y*shadow_ray.direction.y + shadow_ray.direction.z*shadow_ray.direction.z;


  // Check parallel ghost surfaces for soft triangles
  for (int i = 0; i < light_sources; i++)  {
    rand_vec = random(rand_vec); // Turn the randomness up another notch for the next shadow
    // In shadow?
    int mask = (!in_shadow(shadow_ray.start, shadow_ray.direction + crush(rand_vec, light_spread), triangle_vertexes, triangle_colors, radius_sq, triangle_n));
    // How much light?
    total_light += mask*(light_color * max(dot(shadow_ray.direction , intersect_normal), 0.0f)) / ( 4.0f * ((float)M_PI) * radius_sq);
    // total_light += mask*(light_color * max(dot(shadow_ray.direction , intersect_normal+5*crush(rand_vec, 1*light_spread)), 0.0f)) / ( 4.0f * ((float)M_PI) * radius_sq); // For noisy surfaces
  }

  return total_light/light_sources;
}

float3 secondary_light(Ray ray, local float3 *triangle_vertexes, local float3 *triangle_normals, local float4 *triangle_colors, int triangle_n, float3 light_pos, const int global_id)  {
  const int bounces = 10; // Max number of light bounces off mirrors and glass 

  Ray perturbed_ray = ray; // Current refracted ray or reflected ray


  for (int b = 0; b < bounces && perturbed_ray.intersect_color.w <= 0.0f; b++ )  { // Keep bouncing until bounce count, or object isnt reflective
    perturbed_ray = (perturbed_ray.intersect_color.w == 0.0f) ? reflect_ray(perturbed_ray) : refract_ray(perturbed_ray); // Reflect or refract??
    single_ray_intersections(&perturbed_ray, triangle_vertexes, triangle_normals, triangle_colors, triangle_n); // Check new intersection point
    int intersected = (perturbed_ray.intersect_triangle != -1 && perturbed_ray.intersect_color.w > 0.0f );
    if (intersected)  {
      // Its a boring material, so calculate light and color of material and dull it a bit
      const float3 light = indirect_light + direct_light(perturbed_ray, triangle_vertexes, triangle_colors,  light_pos, triangle_n, perturbed_ray.intersect_normal, global_id);
      return 0.9f*light*perturbed_ray.intersect_color.xyz;
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
{      
  // Copy triangle geometry to local memory
  event_t e = async_work_group_copy(LOC_triangle_vertexes, triangle_vertexes, triangle_n*3, 0); 
  e         = async_work_group_copy(LOC_triangle_normals,  triangle_normals,  triangle_n,   0);
  e         = async_work_group_copy(LOC_triangle_colors,   triangle_colors,   triangle_n,   0);

  const int x = get_global_id(0);  // Pixel to work on
  const int y = get_global_id(1);  
  const int global_id = y*SCREEN_WIDTH+x;
  float3 final_color_total = (float3) (0.0f);


  const float3 base_dir = (float3) (x*rays_x - (SCREEN_WIDTH*rays_x)/2.0f, y*rays_y - (SCREEN_HEIGHT*rays_y)/2.0f, focal_length);



  Ray rays[aa_rays];
  const float3 r0 = rot_matrix[0];
  const float3 r1 = rot_matrix[1];
  const float3 r2 = rot_matrix[2];
  // Initialise all anti-aliasing rays
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
  // Calculate first intersection for all rays

  for (char r = 0; r < aa_rays; r++)  {
    single_ray_intersections(&rays[r], LOC_triangle_vertexes, LOC_triangle_normals, LOC_triangle_colors, triangle_n);
  }



  // For each ray
  for (char r = 0; r < aa_rays; r++)  { 
    if (rays[r].intersect_triangle != -1) { // If intersection

      if(rays[r].intersect_color.w <= 0.0f)  { // Mirror or glass
        final_color_total += secondary_light(rays[r], LOC_triangle_vertexes, LOC_triangle_normals, LOC_triangle_colors, triangle_n, light_pos, global_id);
      } else { // Plain material
        const float3 first_light  = direct_light(rays[r], LOC_triangle_vertexes, LOC_triangle_colors,  light_pos, triangle_n, rays[r].intersect_normal, global_id);
        final_color_total += rays[r].intersect_color.xyz*(indirect_light + first_light);
      }
    }
  }
  // Store in screen buffer
  color_pixel(screen_buffer, (short)x, (short)y, final_color_total/((float)aa_rays));
}