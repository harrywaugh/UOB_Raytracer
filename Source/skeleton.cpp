#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>
#include <glm/gtc/type_ptr.hpp>
#include <math.h>
using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::cos;
using glm::sin;
using glm::length;

SDL_Event event;

#define SCREEN_WIDTH 1400
#define SCREEN_HEIGHT 1400
#define FULLSCREEN_MODE false


struct Intersection {
  vec4 position;
  float distance;
  int triangle_index;
};

float focal_length = 500.0;
vec4  camera_position(0.0, 0.0, -3.0, 1.0);
float pitch = 0.0f;
float yaw = 0.0f;

vec4 light_position(0, -0.5, -0.7, 1.0);
vec3 light_color = 14.f * vec3(1, 1, 1);
vec3 indirect_light = 0.5f * vec3(1, 1, 1);

bool quit = false;

vector<Triangle> triangles;

bool update();
void draw(screen* screen);

int main(int argc, char* argv[]) {

  // Initialise screen
  screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);

  // Load Cornell Box
  LoadTestModel(triangles);

  // Draw initial scene
  draw(screen);
  SDL_Renderframe(screen);

  // While user hasn't quit
  while (!quit) {
    // If there is an update to the scene, then draw changes. Check if user wants to quit
    if (update())  {
      draw(screen);
      SDL_Renderframe(screen);
    }
  }

  SDL_SaveImage(screen, "screenshot.bmp");

  KillSDL(screen);
  return 0;
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

    // mat3 A(-d, e1, e2);
    // vec3 x = glm::inverse(A) * b;

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

float square(float x) {
  return x * x;
}

float max(float x, float y) {
  if (x >= y) return x;
  else return y;
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
void draw(screen* screen) {
  // Clear the buffer
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  mat4 R;
  for (int y = 0; y < screen->height; y++) {
    for (int x = 0; x < screen->width; x++) {
      // We only need to implement rotation around y and x axis
      float r[16] = {cos(yaw),  sin(pitch)*sin(yaw),   sin(yaw)*cos(pitch),  1.0f,
                     0.0f,      cos(pitch),           -sin(pitch),           1.0f,
                    -sin(yaw),  cos(yaw)*sin(pitch),   cos(pitch)*cos(yaw),  1.0f,
                     1.0f,      1.0f,                  1.0f,                 1.0f};
      mat4 R;
      memcpy(glm::value_ptr(R), r, sizeof(r));

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
        PutPixelSDL(screen, x, y, vec3(0.0,0.0,0.0));
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
