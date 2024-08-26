/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Kobe Kodachi
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <glm/glm.hpp>
#include <cmath>
#include <iostream>
// #include <omp.h>
// #include <cuda_runtime.h>
#include <vector>
#include <random>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

// The different display modes.
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// While solving the homework, it is useful to make the below values smaller for debugging purposes.
// The still images that you need to submit with the homework should be at the below resolution (640x480).
// However, for your own purposes, after you have solved the homework, you can increase those values to obtain higher-resolution images.
#define WIDTH 640
#define HEIGHT 480

// The field of view of the camera, in degrees.
#define fov 60.0
#define PI 3.14159265
#define PRECISION 1e-3
#define ANTIALIASING true
#define SOFTSHADOW false
#define SHADOWSAMPLE 10.0

// Buffer to store the image when saving it to a JPEG.
unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

bool equals(const double& x,const double& y){
  return abs(x - y) < PRECISION;
}

const double aspect_ratio = (double) WIDTH / (double) HEIGHT;
const double topX = -aspect_ratio * tan(fov * 0.5 * PI / 180.0); // top x coordinate of img plane
const double topY = tan(fov * 0.5 * PI / 180.0); // top y coordinate of img plane
const double x_delta = (2.0 * abs(topX)) / (double) WIDTH; // change in x for each pixel in width
const double y_delta = (2.0 * abs(topY)) / (double) HEIGHT; // change in y for each pixel in width

class Ray {
    public:
    glm::vec3 point; // origin point, camera origin point or intersection for shadow
    glm::vec3 dir; // direction of ray, to image plane or light source

    // create a shadow ray from a point and direction
    Ray(glm::vec3 p,glm::vec3 d){
      point = p;
      dir = d;
      
    }

    // create a ray from the camera to location of pixel(x,y)
    Ray(double x,double y){
      glm::vec3 origin(0.0,0.0,0.0);
      // double dirX = topX + (((double)x + 0.5) * x_delta); // top left of img x + middle of pixel * change per pixel x
      // double dirY = topY + (((double)y + 0.5) * y_delta); // top left of img y + middle of pixel * change per pixel y
      // adjust center pixel to range 0 to 1, expand and shift range to -1 to 1, multiple by distance to edge
      double dirX = ((x / (double) WIDTH) * 2.0 - 1) * aspect_ratio * topY;
      double dirY = ((y / (double) HEIGHT) * 2.0 - 1) * topY;

      glm::vec3 d(dirX,dirY,-1.0);
      point = origin;
      dir = glm::normalize(d);
    }
    // calculate intersection using t0,1 = (-b +/- sqrt(b^2 - 4c)) / 2, if min(t0,t1) > 0 intersection
    bool sphereIntersect(const Sphere& sphere, glm::vec3& i){
      glm::vec3 sphereDist = point - glm::vec3(sphere.position[0],sphere.position[1],sphere.position[2]);
      double b = 2.0 * glm::dot(dir,sphereDist);
      double c = glm::dot(sphereDist,sphereDist) - ((double)sphere.radius * (double)sphere.radius);
      double disc = (b * b) - (4.0 * c);
      if (disc < -PRECISION) return false;
      double t0,t1;
      if (equals(disc,0.0)){
        t0 = -b * 0.5;
        t1 = t0;
      } else {
        t0 = (-b + sqrt(disc)) * 0.5;
        t1 = (-b - sqrt(disc)) * 0.5;
      }
      if ((t0 < -PRECISION) && (t1 < -PRECISION)) return false;
      else if ((t1 < t0-PRECISION) && (t1 > PRECISION)) t0 = t1;
      // intersection point
      i = point + (dir * (float)t0);
      return true;
    }

    // calculate intersection and barycentric coordinates with area 
    bool triIntersect(const Triangle& tri, glm::vec3& i,glm::vec3& bary){
      glm::vec3 A(tri.v[0].position[0],tri.v[0].position[1],tri.v[0].position[2]);
      glm::vec3 B(tri.v[1].position[0],tri.v[1].position[1],tri.v[1].position[2]);
      glm::vec3 C(tri.v[2].position[0],tri.v[2].position[1],tri.v[2].position[2]);
      glm::vec3 norm = glm::normalize(glm::cross(B-A,C-A));

      double denom = glm::dot(norm,dir); // n * d
      if (equals(denom,0.0)) return false;
      double t = -(glm::dot(norm,point) + glm::dot(-A,norm)) / denom; // -(n * p0 + d) / (n * dir)
      if ((t < -PRECISION) || equals(t,0.0)) return false;
      // intersection point
      i = point + (dir * (float)t);

      double C0C1C2 = glm::length(glm::cross((B-A),(C-A)));
      double CC1C2 = glm::length(glm::cross((B-i),(C-i)));
      double C0CC2 = glm::length(glm::cross((i-A),(C-A)));
      double C0C1C = glm::length(glm::cross((B-A),(i-A)));
      double alpha = CC1C2 / C0C1C2;
      double beta = C0CC2 / C0C1C2;
      double gamma = C0C1C / C0C1C2;

      if ((alpha < -PRECISION) || (beta < -PRECISION) || (gamma < -PRECISION) || !equals(alpha+beta+gamma,1.0)) return false;
      // barycentric coordinates
      bary = glm::vec3(alpha,beta,gamma);
      return true;      
    }
};

// check for intersections with spheres for a given Ray and light
bool sphereCheck(Ray& shadow,int i,const glm::vec3& light){
  glm::vec3 si(0.0);
  for (int k=0;k<num_spheres;k++){
    if (i == k) continue;
    bool sIntersect = shadow.sphereIntersect(spheres[k],si);
    double sphereDist = glm::length(si - shadow.point); // distance from og intersection point to new
    double lightDist = glm::length(light - shadow.point); // distance from light to og intersection point
    // check intersection of all spheres except current where blocking sphere is between light and current
    if (sIntersect && (sphereDist < lightDist)){ 
      return false;
    }
  }
  return true;
}

// check for intersections with triangles for a given Ray and light
bool triCheck(Ray& shadow,int i,const glm::vec3& light){
  glm::vec3 si(0.0);
  glm::vec3 garbage(0.0);
  for (int k=0;k<num_triangles;k++){
    if (i == k) continue;
    bool tIntersect = shadow.triIntersect(triangles[k],si,garbage);
    double triDist = glm::length(si - shadow.point); // distance from og intersection point to new
    double lightDist = glm::length(light - shadow.point); // distance from light to og intersection point
    // check intersection of all triangles where blocking triangle is between light and current
    if (tIntersect && (triDist < lightDist)){
      return false;
    }
  }
  return true;
}

// implements soft shadows by randomly generate 10 sublights for each light source
// can vary radius to change effects of lighting, currently x1.0
std::vector<Light> getSublights(const Light l){
  std::vector<Light> sublights;
  // randomly generate numbers 0.0-1.0
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0.0, 1.0);
  for (int i=0;i<SHADOWSAMPLE;i++){
    // randomly generate theta, phi, and radius
    double theta = dis(gen) * 2 * PI;
    double phi = dis(gen) * PI;
    double r = dis(gen) * 1.0;
    Light sublight;
    // spherical to cartesian
    sublight.position[0] = l.position[0] + r * sin(phi) * cos(theta);
    sublight.position[1] = l.position[1] + r * sin(phi) * sin(theta);
    sublight.position[2] = l.position[2] + r * cos(phi);
    // divide original light by number of sublights
    sublight.color[0] = l.color[0] / SHADOWSAMPLE;
    sublight.color[1] = l.color[1] / SHADOWSAMPLE;
    sublight.color[2] = l.color[2] / SHADOWSAMPLE;
    sublights.push_back(sublight);
  }
  return sublights;
}

// I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ sh)
glm::vec3 sphereShading(const Light& l,const Sphere& s,const glm::vec3& i){
  glm::vec3 lightColor(l.color[0],l.color[1],l.color[2]);
  glm::vec3 diffuse(s.color_diffuse[0],s.color_diffuse[1],s.color_diffuse[2]);
  glm::vec3 lightDir = glm::normalize(glm::vec3(l.position[0],l.position[1],l.position[2]) - i);
  glm::vec3 normal = glm::normalize(i - glm::vec3(s.position[0],s.position[1],s.position[2]));
  glm::vec3 specular(s.color_specular[0],s.color_specular[1],s.color_specular[2]);
  glm::vec3 viewDir = glm::normalize(-i);

  double LdotN = glm::dot(lightDir, normal);
  LdotN = glm::clamp(LdotN,0.0,1.0);
  glm::vec3 reflection = 2 * (float)LdotN * normal - lightDir;
  reflection = glm::normalize(reflection);
  double RdotV = glm::dot(reflection,viewDir);
  RdotV = glm::clamp(RdotV,0.0,1.0);
  RdotV = pow(RdotV,s.shininess);

  return lightColor * ((diffuse * (float)LdotN) + (specular * (float)RdotV));
}

// shade triangles using interpolated values with barycentric coordinates
// I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ sh)
glm::vec3 triShading(const Light& l,const Triangle& t,const glm::vec3 i,const glm::vec3 bary){
  glm::vec3 lightColor(l.color[0],l.color[1],l.color[2]);
  // interpolate values using barycentric coordinates
  glm::vec3 diffA(t.v[0].color_diffuse[0],t.v[0].color_diffuse[1],t.v[0].color_diffuse[2]);
  glm::vec3 diffB(t.v[1].color_diffuse[0],t.v[1].color_diffuse[1],t.v[1].color_diffuse[2]);
  glm::vec3 diffC(t.v[2].color_diffuse[0],t.v[2].color_diffuse[1],t.v[2].color_diffuse[2]);
  glm::vec3 interDiff = glm::vec3(bary.x * diffA + bary.y * diffB + bary.z * diffC);

  glm::vec3 lightDir = glm::normalize(glm::vec3(l.position[0],l.position[1],l.position[2]) - i);

  glm::vec3 normA(t.v[0].normal[0],t.v[0].normal[1],t.v[0].normal[2]);
  glm::vec3 normB(t.v[1].normal[0],t.v[1].normal[1],t.v[1].normal[2]);
  glm::vec3 normC(t.v[2].normal[0],t.v[2].normal[1],t.v[2].normal[2]);
  glm::vec3 interNorm = glm::normalize(glm::vec3(bary.x * normA + bary.y * normB + bary.z * normC));

  glm::vec3 specA(t.v[0].color_specular[0],t.v[0].color_specular[1],t.v[0].color_specular[2]);
  glm::vec3 specB(t.v[1].color_specular[0],t.v[1].color_specular[1],t.v[1].color_specular[2]);
  glm::vec3 specC(t.v[2].color_specular[0],t.v[2].color_specular[1],t.v[2].color_specular[2]);
  glm::vec3 interSpec = glm::vec3(bary.x * specA + bary.y * specB + bary.z * specC);

  glm::vec3 viewDir = glm::normalize(-i);

  double shininess = bary.x * t.v[0].shininess + bary.y * t.v[1].shininess + bary.z * t.v[2].shininess;

  double LdotN = glm::dot(lightDir, interNorm);
  LdotN = glm::clamp(LdotN,0.0,1.0);

  glm::vec3 reflection = 2 * (float)LdotN * interNorm - lightDir;
  reflection = glm::normalize(reflection);

  double RdotV = glm::dot(reflection,viewDir);
  RdotV = glm::clamp(RdotV,0.0,1.0);
  RdotV = pow(RdotV,shininess);

  return lightColor * ((interDiff * (float)LdotN) + (interSpec * (float)RdotV));
}

//perform intersections with the scene for a given ray from the camera
glm::vec3 sceneIntersect(Ray& r){
  glm::vec3 col(1.0);
  double curr_near = -1e6; // set to small intial value to be overwritten
  glm::vec3 garbage(0.0);
  // intersections for spheres
  for (int i=0;i<num_spheres;i++){
    glm::vec3 inter(0.0,0.0,-1e6); // set z so don't overwrite unintentionally
    bool intersects = r.sphereIntersect(spheres[i],inter);
    if (intersects && (inter.z > curr_near)){ // check intersection and new nearest
      col = glm::vec3(0.0);
      for (int j=0;j<num_lights;j++){
        if (SOFTSHADOW){
          std::vector<Light> sublights = getSublights(lights[j]);
          glm::vec3 subCol(0.0);
          for (int k=0;k<SHADOWSAMPLE;k++){
            glm::vec3 light(sublights[k].position[0],sublights[k].position[1],sublights[k].position[2]);
            glm::vec3 dir = glm::normalize(light - inter);
            Ray shadow(inter,dir);
            if (sphereCheck(shadow,i,glm::vec3(sublights[k].position[0],sublights[k].position[1],sublights[k].position[2])) 
            && triCheck(shadow,-1,glm::vec3(sublights[k].position[0],sublights[k].position[1],sublights[k].position[2]))){
              col += sphereShading(sublights[k],spheres[i],inter);
            }
          }
        } else{
          glm::vec3 light(lights[j].position[0],lights[j].position[1],lights[j].position[2]);
          glm::vec3 dir = glm::normalize(light - inter);
          Ray shadow(inter,dir);
          if (sphereCheck(shadow,i,light) && triCheck(shadow,-1,light)){
            col += sphereShading(lights[j],spheres[i],inter);
          }
        }
      }
      curr_near = inter.z;
    }
  }
  // intersections for triangles
  for (int i=0;i<num_triangles;i++){ // intersections for triangles
    glm::vec3 inter(0.0,0.0,-1e6); // set z so don't overwrite unintentionally
    glm::vec3 bary(0.0);
    bool intersects = r.triIntersect(triangles[i],inter,bary);
    if (intersects && (inter.z > curr_near)){ // check intersection and new nearest
      col = glm::vec3(0.0);
      for (int j=0;j<num_lights;j++){
        if (SOFTSHADOW){
          std::vector<Light> sublights = getSublights(lights[j]);
          glm::vec3 subCol(0.0);
          for (int k=0;k<SHADOWSAMPLE;k++){
            glm::vec3 light(sublights[k].position[0],sublights[k].position[1],sublights[k].position[2]);
            glm::vec3 dir = glm::normalize(light - inter);
            Ray shadow(inter,dir);
            if (sphereCheck(shadow,-1,glm::vec3(sublights[k].position[0],sublights[k].position[1],sublights[k].position[2])) 
            && triCheck(shadow,i,glm::vec3(sublights[k].position[0],sublights[k].position[1],sublights[k].position[2]))){
              col += triShading(sublights[k],triangles[i],inter,bary);
            }
          }
        } else {
          glm::vec3 light(lights[j].position[0],lights[j].position[1],lights[j].position[2]);
          glm::vec3 dir = glm::normalize(light - inter);
          Ray shadow(inter,dir);
          if (sphereCheck(shadow,-1,light) && triCheck(shadow,i,light)){
            col += triShading(lights[j],triangles[i],inter,bary);
          }
        }
      }
      curr_near = inter.z;
    }
  }
  return col + glm::vec3(ambient_light[0],ambient_light[1],ambient_light[2]);
}

void draw_scene()
{
  // omp_set_num_threads(64);
  // #pragma omp parallel for collapse(2)
  // {
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    // Do not worry about this usage of OpenGL. This is here just so that we can draw the pixels to the screen,
    // after their R,G,B colors were determined by the ray tracer.
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      if (!ANTIALIASING){
        Ray ray(x+0.5,y+0.5); // create ray using center of pixel
        glm::vec3 col = sceneIntersect(ray); // intersect with scene
        col = glm::clamp(col,glm::vec3(0.0),glm::vec3(1.0));
        plot_pixel(x,y,col.x * 255,col.y * 255,col.z * 255);
      } else if (ANTIALIASING){
        glm::vec3 col(0.0);
        glm::vec3 curr(0.0);
        glm::vec3 prev(1.0);
        glm::vec3 holder(0.0);
        float coords[][2] = { // divide pixel into 9 pieces for antialiasing
          {0.25,0.25},{0.25,0.50},{0.25,0.75},
          {0.50,0.25},{0.50,0.50},{0.50,0.75},
          {0.75,0.25},{0.75,0.50},{0.75,0.75},
        };
        bool close = false;
        for (int k=0;k<9;k++){ // trace each ray
          Ray ray(x+coords[k][0],y+coords[k][1]);
          prev = curr;
          curr = sceneIntersect(ray);
          // break if colors close and not on last iteration
          if (equals(prev.x,curr.x) && equals(prev.y,curr.y) && equals(prev.z,curr.z) && k != 8){
            close = true;
            break;
          }
          holder += curr;
        }
        if (close){ // computer normally if colors close
          Ray ray(x+0.5,y+0.5);
          col = sceneIntersect(ray);
        } else{ // average results
          col = holder;
          col /= 9.0;
        }
        col = glm::clamp(col,glm::vec3(0.0),glm::vec3(1.0));
        plot_pixel(x,y,col.x * 255,col.y * 255,col.z * 255);
      }      
    }
    glEnd();
    glFlush();
  // }
  }
  printf("Ray tracing completed.\n"); 
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in saving\n");
  else 
    printf("File saved successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parsing error; abnormal program abortion.\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  if (!file)
  {
    printf("Unable to open input file %s. Program exiting.\n", argv);
    exit(0);
  }

  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  // Hack to make it only draw once.
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

