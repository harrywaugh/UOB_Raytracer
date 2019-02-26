#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>
#include <glm/gtc/type_ptr.hpp>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::cos;
using glm::sin;

SDL_Event event;

#define SCREEN_WIDTH 500
#define SCREEN_HEIGHT 500
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


vector<Triangle> triangles;

bool update();
void draw(screen* screen);

int main(int argc, char* argv[]) {
  screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);

  LoadTestModel(triangles);

  while (update()) {
    draw(screen);
    SDL_Renderframe(screen);
  }

  SDL_SaveImage(screen, "screenshot.bmp");

  KillSDL(screen);
  return 0;
}

bool closest_intersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closest_intersection) {
  float current_t = std::numeric_limits<float>::max();
  vec3 d = vec3(dir.x, dir.y, dir.z);
  for (uint i = 0; i < triangles.size(); i++) {
    vec4 v0 = triangles.at(i).v0;
    vec4 v1 = triangles.at(i).v1;
    vec4 v2 = triangles.at(i).v2;
    vec3 e1 = vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
    vec3 e2 = vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
    vec3 b = vec3(start.x-v0.x,start.y-v0.y,start.z-v0.z);
    mat3 A(-d, e1, e2);
    vec3 x = glm::inverse(A) * b;

    if (x.x >= 0 && x.y >= 0 && x.z >= 0 && (x.y + x.z) <= 1 && x.x < current_t) {
      vec3 position = vec3(v0.x, v0.y, v0.z) + (x.y * e1) + (x.z * e2);
      closest_intersection.position = vec4(position.x, position.y, position.z, 1.0);
      closest_intersection.distance = x.x;
      closest_intersection.triangle_index = i;
      current_t = x.x;
    }
  }
  if (current_t == std::numeric_limits<float>::max()) return false;
  return true;
}

/*Place your drawing here*/
void draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  mat4 R;
  for (int y = 0; y < screen->height; y++) {
    for (int x = 0; x < screen->width; x++) {
      float r[16] = {cos(yaw),                  sin(yaw),            0.0f,       1.0f,
                    -sin(yaw)*glm::cos(pitch),  cos(yaw)*cos(pitch), sin(pitch), 1.0f,
                     sin(yaw)*glm::sin(pitch), -cos(yaw)*sin(pitch), cos(pitch), 1.0f,
                     1.0f,                      1.0f,                1.0f,       1.0f};
      mat4 R;
      memcpy(glm::value_ptr(R), r, sizeof(r));
      vec4 d = vec4(x - screen->width/2, y - screen->height/2, focal_length, 1.0);
      d = R * d;

      Intersection intersection;
      if (closest_intersection(camera_position, d, triangles, intersection)) {
        PutPixelSDL(screen, x, y, triangles.at(intersection.triangle_index).color);
      } else {
        PutPixelSDL(screen, x, y, vec3(0.0,0.0,0.0));
      }
    }
  }
}

/*Place updates of parameters here*/
bool update() {
  // static int t = SDL_GetTicks();
  // /* Compute frame time */
  // int t2 = SDL_GetTicks();
  // float dt = float(t2-t);
  // t = t2;

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
      return false;
    } else if (e.type == SDL_KEYDOWN) {
	    int key_code = e.key.keysym.sym;
	    switch(key_code) {
	      case SDLK_UP:
          pitch += 0.2;
      		break;
	      case SDLK_DOWN:
          pitch -= 0.2;
      		break;
	      case SDLK_LEFT:
      		yaw += 0.2;
      		break;
	      case SDLK_RIGHT:
      		yaw -= 0.2;
      		break;
        case SDLK_w:
          camera_position.z += 0.2;
          break;
        case SDLK_s:
          camera_position.z -= 0.2;
          break;
        case SDLK_a:
          camera_position.x -= 0.2;
          break;
        case SDLK_d:
          camera_position.x += 0.2;
          break;
        case SDLK_i:
          focal_length += 10;
          break;
        case SDLK_o:
          focal_length -= 10;
          break;
        case SDLK_ESCAPE:
      		/* Move camera quit */
      		return false;
      }
	  }
  }
  return true;
}
