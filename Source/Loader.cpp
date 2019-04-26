#include <iostream>
#include <fstream>
#include <sstream>
#include "TestModelH.h"

using namespace std;
using glm::vec3;
using glm::vec4;
using glm::mat4;

vector<Triangle> load_obj(string filename) {

    // Defines colors:
  vec3 red(    0.75f, 0.15f, 0.15f );
  vec3 yellow( 0.75f, 0.75f, 0.15f );
  vec3 green(  0.15f, 0.75f, 0.15f );
  vec3 cyan(   0.15f, 0.75f, 0.75f );
  vec3 blue(   0.15f, 0.15f, 0.75f );
  vec3 mirror( -1.0f, -1.0f, -1.0f);
  vec3 purple( 0.75f, 0.15f, 0.75f );
  vec3 white(  0.75f, 0.75f, 0.75f );
  vec3 grey(  0.3f, 0.3f, 0.3f );

  ifstream source(filename);
  vector<vec4> vertices;
  vector<Triangle> triangles;

  string line;

  while (getline(source, line)) {
    istringstream in(line);

    string s;
    in >> s;

    if (s == "v") {
      float x, y, z;
      in >> x >> y >> z;
      vertices.push_back(vec4(1.5f*x, 1.5f*y, 1.5f*z, 1.f));
    } else if (s == "f") {
      int v1, v2, v3;
      in >> v1 >> v2 >> v3;
      Triangle triangle = Triangle(vertices[v1-1], vertices[v2-1], vertices[v3-1], grey);

      vec4 translate(-0.4f, 1.15f, -0.7f, 1.0f);

      triangle.v1 = (-1.f)*triangle.v1 + translate;
      triangle.v2 = (-1.f)*triangle.v2 + translate;
      triangle.v0 = (-1.f)*triangle.v0 + translate;

      triangles.push_back(triangle);
    }
  }

  return triangles;
}