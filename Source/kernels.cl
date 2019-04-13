#pragma OPENCL EXTENSION cl_khr_fp64 : enable

constant float focal_length = 500.0;
constant float3 indirect_light = (float3)(0.5f, 0.5f, 0.5f);
#define SCREEN_WIDTH 1400
#define SCREEN_HEIGHT 1400
/////READ ONLY BUFFERS

kernel void draw(global uint  *screen_buffer,    global float3 *triangle_vertexes,   global float4 *triangle_normals,
				 global float3 *triangle_colors, global float *rot_matrix,           float4 camera_pos)
{         /* accumulated magnitudes of velocity for each cell */
  const short x = get_global_id(0);
  const short y = get_global_id(1);
  if(x==0 && y==0)  {
  	printf("Triangle 0: v0x  %f v1x %f v2x %f\n", triangle_vertexes[0].s0, triangle_vertexes[1].s0, triangle_vertexes[2].s0);
  	printf("Triangle 0: norm  %f col %f\n",       triangle_normals[0].s0, triangle_colors[0].s0);
  	printf("Triangle 1: v0x  %f v1x %f v2x %f\n", triangle_vertexes[3].s0, triangle_vertexes[4].s0, triangle_vertexes[5].s0);
  	printf("Triangle 1: norm  %f col %f\n",       triangle_normals[1].s0, triangle_colors[1].s0);

  }
}

// bool closest_intersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closest_intersection) {
//   // Set closest intersection to be the max float value
//   float current_t = std::numeric_limits<float>::max();
//   // Make 4D ray into 3D ray
//   vec3 d = vec3(dir.x, dir.y, dir.z);
//   for (uint i = 0; i < triangles.size(); i++) {
//     // Define two corners of triangle relative to the other corner
//     vec4 v0 = triangles.at(i).v0;
//     vec4 v1 = triangles.at(i).v1;
//     vec4 v2 = triangles.at(i).v2;
//     vec3 e1 = vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
//     vec3 e2 = vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
//     vec3 b = vec3(start.x-v0.x,start.y-v0.y,start.z-v0.z);

//     // Cramers, might be det repeated computation..?
//     float detA = glm::determinant(mat3(-d, e1, e2));
//     float detA0 = glm::determinant(mat3(b, e1, e2));
//     float detA1 = glm::determinant(mat3(-d, b, e2));
//     float detA2 = glm::determinant(mat3(-d, e1, b));

//     vec3 x(detA0/detA, detA1/detA, detA2/detA);

//     // If ray goes through triangle, and is the closest triangle
//     if (x.x >= 0 && x.y >= 0 && x.z >= 0 && (x.y + x.z) <= 1 && x.x < current_t) {
//       vec3 position = vec3(v0.x, v0.y, v0.z) + (x.y * e1) + (x.z * e2);

//       closest_intersection.position = vec4(position.x, position.y, position.z, 1.0);
//       closest_intersection.distance = length(x.x * d);
//       closest_intersection.triangle_index = i;
//       current_t = x.x;
//     }
//   }
//   if (current_t == std::numeric_limits<float>::max()) return false;
//   return true;
// }

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

// // Place your drawing here
// void draw(screen* screen, t_ocl ocl) {
//   cl_int err;
//   // Clear the buffer
//   // memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
//   mat4 R;

//   float r[16] = {cos(yaw),  sin(pitch)*sin(yaw),   sin(yaw)*cos(pitch),  1.0f,
//                0.0f,      cos(pitch),           -sin(pitch),             1.0f,
//               -sin(yaw),  cos(yaw)*sin(pitch),   cos(pitch)*cos(yaw),    1.0f,
//                1.0f,      1.0f,                  1.0f,                   1.0f};
//   memcpy(glm::value_ptr(R), r, sizeof(r));

//   err = clSetKernelArg(ocl.draw, 2, sizeof(cl_mem), &ocl.rotation_matrix_buffer);
//   checkError(err, "setting draw arg 2", __LINE__);
//   err = clSetKernelArg(ocl.draw, 3, sizeof(cl_float4), &camera_position);
//   checkError(err, "setting draw arg 3", __LINE__);


//   for (int y = 0; y < screen->height; y++) {
//     for (int x = 0; x < screen->width; x++) {
//       // Declare ray for given position on the screen. Rotate ray by current view angle
//       vec4 d = vec4(x - screen->width/2, y - screen->height/2, focal_length, 1.0);
//       d = R * d;

//       // Find intersection point with closest geometry. If no intersection, paint the abyss
//       Intersection intersection;
//       if (closest_intersection(camera_position, d, triangles, intersection)) {
//         // If the ray drawn does intersect with geometry then draw the correct
//         // colour returned by direct_light()
//         // Get colour of the triangle the light has hit
//         vec3 p = triangles.at(intersection.triangle_index).color;
//         vec3 final_color = p*(direct_light(intersection) + indirect_light);
//         PutPixelSDL(screen, x, y, final_color);
//       } else {
//         // Otherwise draw black
//         PutPixelSDL(screen, x, y, vec3(0.0f, 0.0f, 0.0f));
//       }
//     }
//   }
// }