#pragma OPENCL EXTENSION cl_khr_fp64 : enable

constant float focal_length = 500.0;
constant float3 indirect_light = (float3)(0.5f, 0.5f, 0.5f);
#define SCREEN_WIDTH 1024
#define SCREEN_HEIGHT 1024

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

inline float dot_product(float3 a, float3 b) {
	return a.x * b.x + a.y * b.y + a.y * b.y;
}

// void PutPixelSDL(global uint *screen_buffer, int x, int y, float3 colour) {
//   if(x<0 || x>=SCREEN_WIDTH || y<0 || y>=SCREEN_HEIGHT)  {
//     printf("apa\n");
//     return;
//   }
//   uint r = (uint) min(max(255*colour.x, 0.f), 255.f);
//   uint g = (uint) min(max(255*colour.y, 0.f), 255.f);
//   uint b = (uint) min(max(255*colour.z, 0.f), 255.f);
//   // if(x<0 && y==0)  {
//   // screen_buffer[y*SCREEN_WIDTH+x] = (128<<24) + (r<<16) + (g<<8) + b;
// }


bool closest_intersection(float3 start, float3 d, global float3 *triangle_vertexes, private Intersection* closest_intersection, int triangle_n) {
  // Set closest intersection to be the max float value
  float current_t = 99999999999999.0f;
  // Make 4D ray into 3D ray
  for (uint i = 0; i < triangle_n; i++) {
    // Define two corners of triangle relative to the other corner
    float3 v0 = triangle_vertexes[i*3];
    float3 v1 = triangle_vertexes[i*3+1];
    float3 v2 = triangle_vertexes[i*3+2];

    float3 e1 = (float3) (v1.x-v0.x,    v1.y-v0.y,    v1.z-v0.z);
    float3 e2 = (float3) (v2.x-v0.x,    v2.y-v0.y,    v2.z-v0.z);
    float3 b  = (float3) (start.x-v0.x, start.y-v0.y, start.z-v0.z);

    // Cramers, might be det repeated computation..?
    float3 A[3]  = {-d, e1, e2};
    float3 A0[3] = {b, e1, e2};
    float3 A1[3] = {-d, b, e2};
    float3 A2[3] = {-d, e1, b};

    float detA  = det(A);
    float detA0 = det(A0);
    float detA1 = det(A1);
    float detA2 = det(A2); 

    float3 x = (float3) (detA0/detA, detA1/detA, detA2/detA);

    // If ray goes through triangle, and is the closest triangle
    if (x.x >= 0 && x.y >= 0 && x.z >= 0 && (x.y + x.z) <= 1 && x.x < current_t) {
      float3 position = ((float3) (v0.x, v0.y, v0.z)) + (x.y * e1) + (x.z * e2);

      closest_intersection->position = (float3) (position.x, position.y, position.z);
      float3 dist_vec = x.x*d;
      closest_intersection->distance = native_sqrt(dist_vec.x*dist_vec.x + dist_vec.y*dist_vec.y + dist_vec.z*dist_vec.z);
      closest_intersection->triangle_index = i;
      current_t = x.x;
    }
  }
  if (current_t == 99999999999999.0f) return false;
  return true;
}


kernel void draw(global uint  *screen_buffer,    global float3 *triangle_vertexes,   global float3 *triangle_normals,
				 global float3 *triangle_colors, global float3 *rot_matrix,           float3 camera_pos, int triangle_n)
{         /* accumulated magnitudes of velocity for each cell */
  const short x = get_global_id(0);
  const short y = get_global_id(1);

  printf("r 0 %f 1 %f\n", rot_matrix[0].x, rot_matrix[1].x);


  // if(x==0 && y==0)  {
  // 	printf("Triangle n %d\n", triangle_n);
  // // 	printf("Triangle 0: norm  %f col %f\n",       triangle_normals[0].s0, triangle_colors[0].s0);
  // // 	printf("Triangle 1: v0x  %f v1x %f v2x %f\n", triangle_vertexes[3].s0, triangle_vertexes[4].s0, triangle_vertexes[5].s0);
  // // 	printf("Triangle 1: norm  %f col %f\n",       triangle_normals[1].s0, triangle_colors[1].s0);

  // }


  // Declare ray for given position on the screen. Rotate ray by current view angle
  float3 d = (float3) (x - SCREEN_WIDTH/2.0, y - SCREEN_HEIGHT/2.0, focal_length);
  d =        (float3) (dot_product(rot_matrix[0], d), dot_product(rot_matrix[1], d), dot_product(rot_matrix[2], d));

  // Find intersection point with closest geometry. If no intersection, paint the abyss
  Intersection intersection;
  if (closest_intersection(camera_pos, d, triangle_vertexes, &intersection, triangle_n)) {
    float3 p = triangle_colors[intersection.triangle_index];
    // vec3 final_color = p*(direct_light(intersection) + indirect_light);
	if( x<0 || x>=SCREEN_WIDTH || y<0 || y>=SCREEN_HEIGHT )  {
		printf("apa\n");
		return;
	}
	uint r = (uint) min(max(255*p.x, 0.f), 255.f);
	uint g = (uint) min(max(255*p.y, 0.f), 255.f);
	uint b = (uint) min(max(255*p.z, 0.f), 255.f);
  	screen_buffer[y*SCREEN_WIDTH+x] = (128<<24) + (r<<16) + (g<<8) + b;
  } else {
    // Otherwise draw black
  	screen_buffer[y*SCREEN_WIDTH+x] = 0;
  }

}



// vec3 direct_light(const Intersection& intersection) {

//   // Vector from the light to the point of intersection
//   vec4 r = light_position - intersection.position;
//   // Distance of the checked point to the light source
//   float radius = length(r);

//   Intersection obstacle_intersection;
//   float threshold = 0.001f;

//   if (closest_intersection(intersection.position + vec4(r.x * threshold, r.y * threshold, r.z * threshold, 1.0f),
//                            r, triangles, obstacle_intersection)) {
//     if (obstacle_intersection.distance < radius) return vec3(0, 0, 0);
//   }

//   // Get the normal of the triangle that the light has hit
//   vec4 n = triangles.at(intersection.triangle_index).normal;
//   // Intensity of the colour, based on the distance from the light
//   vec3 D = (vec3) (light_color * max(glm::dot(r, n) , 0)) / (float) (4 * M_PI * radius * radius);

//   return D;
// }

