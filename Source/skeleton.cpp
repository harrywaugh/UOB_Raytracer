#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

struct Intersection {
  vec4 position;
  float distance;
  int triangleIndex;
};

float focalLength = 500.0;
vec4 cameraPos(0.0, 0.0, -10.0, 1.0);

vector<Triangle> triangles;

bool Update();
void Draw(screen* screen);

int main(int argc, char* argv[]) {

  screen *screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE);

  LoadTestModel(triangles);

  while (Update()) {
    Draw(screen);
    SDL_Renderframe(screen);
  }

  SDL_SaveImage(screen, "screenshot.bmp");

  KillSDL(screen);
  return 0;
}

bool ClosestIntersection(vec4 start, vec4 dir, const vector<Triangle>& triangles, Intersection& closestIntersection) {
  float m = std::numeric_limits<float>::max();
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

    if (x.x >= 0 && x.y > 0 && x.z > 0 && (x.y + x.z) < 1) {
      vec3 v03 = vec3(v0.x, v0.y, v0.z);
      vec3 pos = v03 + (x.y * e1) + (x.z * e2);
      closestIntersection.position = vec4(pos.x, pos.y, pos.z, 1.0);
      closestIntersection.distance = x.x;
      closestIntersection.triangleIndex = i;
      return true;
    }
  }
  return false;
}

/*Place your drawing here*/
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  for (int y = 0; y < screen->height; y++) {
    for (int x = 0; x < screen->width; x++) {
      vec4 d = vec4(x - screen->width/2, y - screen->height/2, focalLength, 1.0);
      Intersection intersection;
      if (ClosestIntersection(cameraPos, d, triangles, intersection)) {
        PutPixelSDL(screen, x, y, triangles.at(intersection.triangleIndex).color);
      } else {
        PutPixelSDL(screen, x, y, vec3(0.0,0.0,0.0));
      }
    }
  }
}

/*Place updates of parameters here*/
bool Update() {
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  printf("pos: %f, %f, %f, %f\n", cameraPos.x, cameraPos.y, cameraPos.z, focalLength);

  SDL_Event e;
  while(SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) {
      return false;
    } else if (e.type == SDL_KEYDOWN) {
	    int key_code = e.key.keysym.sym;
	    switch(key_code) {
	      case SDLK_UP:
      		/* Move camera forward */
      		break;
	      case SDLK_DOWN:
      		/* Move camera backwards */
      		break;
	      case SDLK_LEFT:
      		/* Move camera left */
      		break;
	      case SDLK_RIGHT:
      		/* Move camera right */
      		break;
	      case SDLK_ESCAPE:
      		/* Move camera quit */
      		return false;
      }
	  }
  }
  return true;
}
