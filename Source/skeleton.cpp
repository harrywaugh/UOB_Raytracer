#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>
#include <glm/gtc/type_ptr.hpp>
#include <math.h>
#include <CL/opencl.h>
#include <chrono> 

using namespace std;
using namespace std::chrono; 
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::cos;
using glm::sin;
using glm::length;

SDL_Event event;


#define OCLFILE "Source/kernels.cl"
#define WORK_SIZE_X 32                 
#define WORK_SIZE_Y 8       


#define SCREEN_WIDTH 3072
#define SCREEN_HEIGHT 3072
#define FULLSCREEN_MODE false

typedef struct
{
  cl_device_id      device;
  cl_context        context;
  cl_command_queue  queue;

  cl_program program;
  //Kernels
  cl_kernel draw;
  cl_kernel average_pixels;

  //Memory Buffers
  cl_mem screen_buffer;
  cl_mem triangles_buffer;
  cl_mem normal_buffer;
  cl_mem color_buffer;
  cl_mem rotation_matrix_buffer;
} t_ocl;
          


struct Intersection {
  vec4 position;
  float distance;
  int triangle_index;
};
float focal_length = 3000.0;
vec4  camera_position(0.0, 0.0, -2.5, 1.0);
cl_float3  cl_camera_position  = {0.0, 0.0, -2.4};

float pitch = 0.0f;
float yaw = 0.0f;
vec4 light_position(0, -0.5, -0.7, 1.0);
vec3 light_color = 14.f * vec3(1, 1, 1);
vec3 indirect_light = 0.5f * vec3(1, 1, 1);
bool quit = false;
vector<Triangle> triangles;


bool update();
void draw(screen* screen, t_ocl ocl);
void checkError(cl_int err, const char *op, const int line);
void die(const char* message, const int line, const char* file);
void opencl_initialise(t_ocl *ocl);
cl_device_id selectOpenCLDevice();
void offload_rendering(screen* screen, t_ocl ocl);

float square(float x) {
  return x * x;
}

float max(float x, float y) {
  if (x >= y) return x;
  else return y;
}

int main(int argc, char* argv[]) {
  t_ocl    ocl; 



  screen *screen = InitializeSDL(SCREEN_WIDTH/2, SCREEN_HEIGHT/2, FULLSCREEN_MODE);

  // Load Cornell Box
  LoadTestModel(triangles);
  printf("Triangles Length size %lu\n",  triangles.size());

  opencl_initialise(&ocl);

  // Draw initial scene
  offload_rendering(screen, ocl);
  SDL_Renderframe(screen);

  // While user hasn't quit
  while (!quit) {
    // If there is an update to the scene, then draw changes. Check if user wants to quit
    if (update())  {
      auto start = high_resolution_clock::now();
      offload_rendering(screen, ocl);
      auto stop = high_resolution_clock::now();
      auto offload_duration = duration_cast<microseconds>(stop - start); 
      cout << "\nOffloaded GPU Rendertime: "<< offload_duration.count() << " micro seconds" <<  endl; 

      // start = high_resolution_clock::now();
      // draw(screen, ocl);
      // stop = high_resolution_clock::now();
      // auto CPU_duration = duration_cast<microseconds>(stop - start); 
      // cout << "CPU Draw Time: "<< CPU_duration.count() << " micro seconds" <<  endl; 

      cout << "GPU is "<< (500000.0f)/((float)offload_duration.count()) << "x faster" <<  endl; 

      SDL_Renderframe(screen);
    }
  }

  SDL_SaveImage(screen, "screenshot.bmp");

  KillSDL(screen);
  return 0;
}

void offload_rendering(screen* screen, t_ocl ocl)  {
  cl_int err;

  float rot_matrix[12] = {cos(yaw),  sin(pitch)*sin(yaw),   sin(yaw)*cos(pitch), 0.0f,
                          0.0f,      cos(pitch),           -sin(pitch),          0.0f,            
                          -sin(yaw), cos(yaw)*sin(pitch),   cos(pitch)*cos(yaw), 0.0f};


  err = clEnqueueWriteBuffer(ocl.queue, ocl.rotation_matrix_buffer, CL_TRUE, 0,
  sizeof(float) * 12, &rot_matrix, 0, NULL, NULL);
  checkError(err, "writing rotation matrix data", __LINE__);

  //Set Camera Position and Rotation Matrix Arguments
  err = clSetKernelArg(ocl.draw, 4, sizeof(cl_mem), &ocl.rotation_matrix_buffer);
  checkError(err, "setting draw arg 4", __LINE__);
  err = clSetKernelArg(ocl.draw, 5, sizeof(cl_float3), &cl_camera_position);
  checkError(err, "setting draw arg 5", __LINE__);
  err = clSetKernelArg(ocl.draw, 6, sizeof(cl_float3), &light_position);
  checkError(err, "setting draw arg 6", __LINE__);
  err = clSetKernelArg(ocl.draw, 8, sizeof(cl_float), &focal_length);
  checkError(err, "setting draw arg 8", __LINE__);

  // Enqueue kernel
  size_t global_size[2] = {SCREEN_WIDTH, SCREEN_HEIGHT};
  size_t work_size[2] = {WORK_SIZE_X, WORK_SIZE_Y};
  err = clEnqueueNDRangeKernel(ocl.queue, ocl.draw , 2, NULL, global_size, work_size, 0, NULL, NULL);
  checkError(err, "enqueueing draw kernel", __LINE__);


  size_t av_global_size[2] = {SCREEN_WIDTH/2, SCREEN_HEIGHT/2};
  err = clEnqueueNDRangeKernel(ocl.queue, ocl.average_pixels , 2, NULL, av_global_size, work_size, 0, NULL, NULL);
  checkError(err, "enqueueing average_pixels kernel", __LINE__);

  // err = clFinish(ocl.queue);
  // checkError(err, "Waiting to finish draw kernel", __LINE__);

  err = clEnqueueReadBuffer(ocl.queue, ocl.screen_buffer, CL_TRUE, 0,
  sizeof(cl_uint) * (SCREEN_WIDTH/2) * (SCREEN_HEIGHT/2), screen->buffer, 0, NULL, NULL);
  checkError(err, "reading screen buffer data", __LINE__);
}

bool closest_intersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closest_intersection) {
  // Set closest intersection to be the max float value
  float current_t = std::numeric_limits<float>::max();
  // Make 4D ray into 3D ray
  vec3 d = vec3(dir.x, dir.y, dir.z);
  for (uint i = 0; i < triangles.size(); i++) {
    // Define two corners of triangle relative to the other corner
    vec4 v0 = triangles.at(i).v0;
    vec4 v1 = triangles.at(i).v1;
    vec4 v2 = triangles.at(i).v2;
    vec3 e1 = vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
    vec3 e2 = vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
    vec3 b = vec3(start.x-v0.x,start.y-v0.y,start.z-v0.z);

    // Cramers, might be det repeated computation..?
    float detA = glm::determinant(mat3(-d, e1, e2));
    float detA0 = glm::determinant(mat3(b, e1, e2));
    float detA1 = glm::determinant(mat3(-d, b, e2));
    float detA2 = glm::determinant(mat3(-d, e1, b));

    vec3 x(detA0/detA, detA1/detA, detA2/detA);

    // If ray goes through triangle, and is the closest triangle
    if (x.x >= 0 && x.y >= 0 && x.z >= 0 && (x.y + x.z) <= 1 && x.x < current_t) {
      vec3 position = vec3(v0.x, v0.y, v0.z) + (x.y * e1) + (x.z * e2);

      closest_intersection.position = vec4(position.x, position.y, position.z, 1.0);
      closest_intersection.distance = length(x.x * d);
      closest_intersection.triangle_index = i;
      current_t = x.x;
    }
  }
  if (current_t == std::numeric_limits<float>::max()) return false;
  return true;
}

vec3 direct_light(const Intersection& intersection) {

  // Vector from the light to the point of intersection
  vec4 r = light_position - intersection.position;
  // Distance of the checked point to the light source
  float radius = length(r);

  Intersection obstacle_intersection;
  float threshold = 0.001f;

  if (closest_intersection(intersection.position + vec4(r.x * threshold, r.y * threshold, r.z * threshold, 1.0f),
                           r, triangles, obstacle_intersection)) {
    if (obstacle_intersection.distance < radius) return vec3(0, 0, 0);
  }

  // Get the normal of the triangle that the light has hit
  vec4 n = triangles.at(intersection.triangle_index).normal;
  // Intensity of the colour, based on the distance from the light
  vec3 D = (vec3) (light_color * max(glm::dot(r, n) , 0)) / (float) (4 * M_PI * radius * radius);

  return D;
}

// Place your drawing here
void draw(screen* screen, t_ocl ocl) {
  // Clear the buffer
  // memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  mat4 R;

  float r[16] = {cos(yaw),  sin(pitch)*sin(yaw),   sin(yaw)*cos(pitch),  1.0f,
               0.0f,      cos(pitch),           -sin(pitch),             1.0f,
              -sin(yaw),  cos(yaw)*sin(pitch),   cos(pitch)*cos(yaw),    1.0f,
               1.0f,      1.0f,                  1.0f,                   1.0f};
  memcpy(glm::value_ptr(R), r, sizeof(r));


  for (int y = 0; y < screen->height; y++) {
    for (int x = 0; x < screen->width; x++) {
      // Declare ray for given position on the screen. Rotate ray by current view angle
      vec4 d = vec4(x - screen->width/2, y - screen->height/2, focal_length, 1.0);
      d = R * d;

      // Find intersection point with closest geometry. If no intersection, paint the abyss
      Intersection intersection;
      if (closest_intersection(camera_position, d, triangles, intersection)) {
        // If the ray drawn does intersect with geometry then draw the correct
        // colour returned by direct_light()
        // Get colour of the triangle the light has hit
        vec3 p = triangles.at(intersection.triangle_index).color;

        vec3 final_color = p*(direct_light(intersection) + indirect_light);

        PutPixelSDL(screen, x, y, final_color);
      } else {
        // Otherwise draw black
        PutPixelSDL(screen, x, y, vec3(0.0f, 0.0f, 0.0f));
      }
    }
  }
}

// Place updates of parameters here
bool update() {
  // static int t = SDL_GetTicks();
  // // Compute frame time
  // int t2 = SDL_GetTicks();
  // float dt = float(t2-t);
  // t = t2;
  // Change scene via key events
  SDL_Event e;
  while(SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
      quit = true;
      return false;
    } else if (e.type == SDL_KEYDOWN) {
      int key_code = e.key.keysym.sym;
      switch(key_code) {
        case SDLK_UP:
          pitch -= 0.1;
          break;
        case SDLK_DOWN:
          pitch += 0.1;
          break;
        case SDLK_LEFT:
          yaw += 0.1;
          break;
        case SDLK_RIGHT:
          yaw -= 0.1;
          break;
        case SDLK_w:
          // camera_position.z += 0.2;
          light_position.z += 0.2;
          break;
        case SDLK_s:
          // camera_position.z -= 0.2;
          light_position.z -= 0.2;
          break;
        case SDLK_a:
          // camera_position.x -= 0.2;
          light_position.x -= 0.2;
          break;
        case SDLK_d:
          // camera_position.x += 0.2;
          light_position.x += 0.2;
          break;
        case SDLK_i:
          focal_length += 10;
          break;
        case SDLK_o:
          focal_length -= 10;
          break;
        case SDLK_ESCAPE:
          quit = true;
          return false;
      }
      return true;
	  }
  }
  return false;
}

#define MAX_DEVICES 32
#define MAX_DEVICE_NAME 1024 

void opencl_initialise(t_ocl *ocl)  {
  cl_int err;
  char*  ocl_src;        /* OpenCL kernel source */
  FILE*   fp;            /* file pointer */
  long   ocl_size;       /* size of OpenCL kernel source */
  char   message[1024];  /* message buffer */


  ocl->device = selectOpenCLDevice();

  // Create OpenCL context
  ocl->context = clCreateContext(NULL, 1, &ocl->device, NULL, NULL, &err);
  checkError(err, "creating context", __LINE__);

    fp = fopen(OCLFILE, "r");
  if (fp == NULL)
  {
    sprintf(message, "could not open OpenCL kernel file: %s", OCLFILE);
    die(message, __LINE__, __FILE__);
  }

  // Create OpenCL command queue
  ocl->queue = clCreateCommandQueue(ocl->context, ocl->device, 0, &err);
  checkError(err, "creating command queue", __LINE__);

  // Load OpenCL kernel source
  fseek(fp, 0, SEEK_END);
  ocl_size = ftell(fp) + 1;
  ocl_src = (char*)malloc(ocl_size);
  memset(ocl_src, 0, ocl_size);
  fseek(fp, 0, SEEK_SET);
  fread(ocl_src, 1, ocl_size, fp);
  fclose(fp);

  // Create OpenCL program
  ocl->program = clCreateProgramWithSource(
    ocl->context, 1, (const char**)&ocl_src, NULL, &err);
  free(ocl_src);
  checkError(err, "creating program", __LINE__);

  // Build OpenCL program
  err = clBuildProgram(ocl->program, 1, &ocl->device, "-cl-fast-relaxed-math -cl-mad-enable", NULL, NULL);
  if (err == CL_BUILD_PROGRAM_FAILURE)
  {
    size_t sz;
    clGetProgramBuildInfo(
      ocl->program, ocl->device,
      CL_PROGRAM_BUILD_LOG, 0, NULL, &sz);
    char *buildlog = (char*)malloc(sz);
    clGetProgramBuildInfo(
      ocl->program, ocl->device,
      CL_PROGRAM_BUILD_LOG, sz, buildlog, NULL);
    fprintf(stderr, "\nOpenCL build log:\n\n%s\n", buildlog);
    free(buildlog);
  }
  checkError(err, "building program", __LINE__);

    // Create OpenCL kernels
  ocl->draw = clCreateKernel(ocl->program, "draw", &err);
  checkError(err, "creating draw kernel", __LINE__);
  ocl->average_pixels = clCreateKernel(ocl->program, "average_pixels", &err);
  checkError(err, "creating average_pixels kernel", __LINE__);

  // Allocate OpenCL buffers
  ocl->screen_buffer          = clCreateBuffer(ocl->context, CL_MEM_READ_WRITE,
                                sizeof(cl_uint)  * SCREEN_WIDTH * SCREEN_HEIGHT, NULL, &err);
  checkError(err, "creating screen buffer", __LINE__);

  ocl->triangles_buffer       = clCreateBuffer(ocl->context, CL_MEM_READ_WRITE,
                                (sizeof(cl_float4) * triangles.size()*3), NULL, &err);
  checkError(err, "creating Triangle buffer", __LINE__);

  ocl->rotation_matrix_buffer = clCreateBuffer(ocl->context, CL_MEM_READ_WRITE,
                                sizeof(cl_float) * 12 , NULL, &err);
  checkError(err, "creating Rot Mat Buffer buffer", __LINE__);
  
  ocl->normal_buffer = clCreateBuffer(ocl->context, CL_MEM_READ_WRITE,
                                sizeof(cl_float4) * triangles.size() , NULL, &err);
  checkError(err, "creating Normal Buffer buffer", __LINE__);
  
  ocl->color_buffer = clCreateBuffer(ocl->context, CL_MEM_READ_WRITE,
                                sizeof(cl_float4) * triangles.size() , NULL, &err);
  checkError(err, "creating Color Buffer buffer", __LINE__);

  err = clSetKernelArg(ocl->average_pixels, 0, sizeof(cl_mem), &ocl->screen_buffer);
  checkError(err, "setting average_pixels arg 0", __LINE__);


  // Set kernel arguments
  err = clSetKernelArg(ocl->draw, 0, sizeof(cl_mem), &ocl->screen_buffer);
  checkError(err, "setting draw arg 0", __LINE__);
  err = clSetKernelArg(ocl->draw, 1, sizeof(cl_mem), &ocl->triangles_buffer);
  checkError(err, "setting draw arg 1", __LINE__);

  err = clSetKernelArg(ocl->draw, 2, sizeof(cl_mem), &ocl->normal_buffer);
  checkError(err, "setting draw arg 2", __LINE__);
  err = clSetKernelArg(ocl->draw, 3, sizeof(cl_mem), &ocl->color_buffer);
  checkError(err, "setting draw arg 3", __LINE__);
  err = clSetKernelArg(ocl->draw, 4, sizeof(cl_mem), &ocl->rotation_matrix_buffer);
  checkError(err, "setting draw arg 4", __LINE__);
  int n = triangles.size();
  err = clSetKernelArg(ocl->draw, 7, sizeof(cl_int), &n);
  checkError(err, "setting draw arg 7", __LINE__);

  err = clSetKernelArg(ocl->draw, 9, sizeof(cl_float4)*triangles.size()*3, NULL);     //Work Item's Local tot_speeds 
  checkError(err, "setting draw arg 9", __LINE__);
  err = clSetKernelArg(ocl->draw, 10, sizeof(cl_float4)*triangles.size(), NULL);     //Work Item's Local tot_speeds 
  checkError(err, "setting draw arg 10", __LINE__);
  err = clSetKernelArg(ocl->draw, 11, sizeof(cl_float4)*triangles.size(), NULL);     //Work Item's Local tot_speeds 
  checkError(err, "setting draw arg 11", __LINE__);

  cl_float4 *triangle_vertexes = (cl_float4 *)malloc(sizeof(cl_float4)*triangles.size()*3);
  cl_float4 *triangle_normals  = (cl_float4 *)malloc(sizeof(cl_float4)*triangles.size());
  cl_float4 *triangle_colors   = (cl_float4 *)malloc(sizeof(cl_float4)*triangles.size());

  for (uint i = 0; i < triangles.size(); i++)  {
    triangle_vertexes[i*3]    = {triangles[i].v0.x,     triangles[i].v0.y,     triangles[i].v0.z,     0.0f};
    triangle_vertexes[i*3+1]  = {triangles[i].v1.x,     triangles[i].v1.y,     triangles[i].v1.z,     0.0f};
    triangle_vertexes[i*3+2]  = {triangles[i].v2.x,     triangles[i].v2.y,     triangles[i].v2.z,     0.0f};
    triangle_normals[i]       = {triangles[i].normal.x, triangles[i].normal.y, triangles[i].normal.z, 0.0f};
    triangle_colors[i]        = {triangles[i].color.x,  triangles[i].color.y,  triangles[i].color.z,  0.0f};
  }

  err = clEnqueueWriteBuffer(ocl->queue, ocl->triangles_buffer, CL_TRUE, 0,
  sizeof(cl_float4) * triangles.size()*3, triangle_vertexes, 0, NULL, NULL);
  checkError(err, "writing triangle buffer data", __LINE__);

  err = clEnqueueWriteBuffer(ocl->queue, ocl->normal_buffer, CL_TRUE, 0,
  sizeof(cl_float4) * triangles.size(), triangle_normals, 0, NULL, NULL);
  checkError(err, "writing triangle buffer data", __LINE__);

  err = clEnqueueWriteBuffer(ocl->queue, ocl->color_buffer, CL_TRUE, 0,
  sizeof(cl_float4) * triangles.size(), triangle_colors, 0, NULL, NULL);
  checkError(err, "writing triangle buffer data", __LINE__);


}
  
void checkError(cl_int err, const char *op, const int line)
{
  if (err != CL_SUCCESS)
  {
    fprintf(stderr, "OpenCL error during '%s' on line %d: %d\n", op, line, err);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
}

void die(const char* message, const int line, const char* file)
{
  fprintf(stderr, "Error at line %d of file %s:\n", line, file);
  fprintf(stderr, "%s\n", message);
  fflush(stderr);
  exit(EXIT_FAILURE);
}
cl_device_id selectOpenCLDevice()
{
  cl_int err;
  cl_uint num_platforms = 0;
  cl_uint total_devices = 0;
  cl_platform_id platforms[8];
  cl_device_id devices[MAX_DEVICES];
  char name[MAX_DEVICE_NAME];

  // Get list of platforms
  err = clGetPlatformIDs(8, platforms, &num_platforms);
  checkError(err, "getting platforms", __LINE__);

  // Get list of devices
  for (cl_uint p = 0; p < num_platforms; p++)
  {
    cl_uint num_devices = 0;
    err = clGetDeviceIDs(platforms[p], CL_DEVICE_TYPE_ALL,
                         MAX_DEVICES-total_devices, devices+total_devices,
                         &num_devices);
    checkError(err, "getting device name", __LINE__);
    total_devices += num_devices;
  }

  // Print list of devices
  printf("\nAvailable OpenCL devices:\n");
  for (cl_uint d = 0; d < total_devices; d++)
  {
    clGetDeviceInfo(devices[d], CL_DEVICE_NAME, MAX_DEVICE_NAME, name, NULL);
    printf("%2d: %s\n", d, name);
  }
  printf("\n");

  // Use first device unless OCL_DEVICE environment variable used
  cl_uint device_index = 0;
  char *dev_env = getenv("OCL_DEVICE");
  if (dev_env)
  {
    char *end;
    device_index = strtol(dev_env, &end, 10);
    if (strlen(end))
      die("invalid OCL_DEVICE variable", __LINE__, __FILE__);
  }

  if (device_index >= total_devices)
  {
    fprintf(stderr, "device index set to %d but only %d devices available\n",
            device_index, total_devices);
    exit(1);
  }

  // Print OpenCL device name
  clGetDeviceInfo(devices[device_index], CL_DEVICE_NAME,
                  MAX_DEVICE_NAME, name, NULL);
  printf("Selected OpenCL device:\n-> %s (index=%d)\n\n", name, device_index);

  return devices[device_index];
}
