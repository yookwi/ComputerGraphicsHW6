#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include<algorithm>

using namespace glm;
using namespace std;

// -------------------------------------------------
// Global Variables
// -------------------------------------------------
#define M_PI 3.1415926535f
#define INF 1'000'000'000
int Width = 512;
int Height = 512;
std::vector<float> OutputImage;
// -------------------------------------------------

void create_scene(int*& gIndexBuffer,vec4*& gVertexBuffer, int* vertexNum, int* triangleNum)
{
	int width = 32;
	int height = 16;
	int gNumVertices = (height - 2) * width + 2;
	int gNumTriangles = (height - 2) * (width - 1) * 2;
	*vertexNum = gNumVertices;
	*triangleNum = gNumTriangles;

	float theta, phi;
	int t;

	// TODO: Allocate an array for gNumVertices vertices.
	gIndexBuffer = new int[3 * gNumTriangles];
	gVertexBuffer = new vec4[gNumVertices];

	t = 0;
	for (int j = 1; j < height - 1; ++j)
	{
		for (int i = 0; i < width; ++i)
		{
			theta = (float)j / (height - 1) * M_PI;
			phi = (float)i / (width - 1) * M_PI * 2;

			float   x = sinf(theta) * cosf(phi);
			float   y = cosf(theta);
			float   z = -sinf(theta) * sinf(phi);

			// TODO: Set vertex t in the vertex array to {x, y, z}.
			gVertexBuffer[t] = {x,y,z, 1};
			t++;
		}
	}

	// TODO: Set vertex t in the vertex array to {0, 1, 0}.
	gVertexBuffer[t] = { 0, 1, 0 , 1};
	t++;

	// TODO: Set vertex t in the vertex array to {0, -1, 0}.
	gVertexBuffer[t] = { 0, -1, 0 , 1};
	t++;

	t = 0;
	for (int j = 0; j < height - 3; ++j)
	{
		for (int i = 0; i < width - 1; ++i)
		{
			gIndexBuffer[t++] = j * width + i;
			gIndexBuffer[t++] = (j + 1) * width + (i + 1);
			gIndexBuffer[t++] = j * width + (i + 1);
			gIndexBuffer[t++] = j * width + i;
			gIndexBuffer[t++] = (j + 1) * width + i;
			gIndexBuffer[t++] = (j + 1) * width + (i + 1);
		}
	}
	for (int i = 0; i < width - 1; ++i)
	{
		gIndexBuffer[t++] = (height - 2) * width;
		gIndexBuffer[t++] = i;
		gIndexBuffer[t++] = i + 1;
		gIndexBuffer[t++] = (height - 2) * width + 1;
		gIndexBuffer[t++] = (height - 3) * width + (i + 1);
		gIndexBuffer[t++] = (height - 3) * width + i;
	}

	// The index buffer has now been generated. Here's how to use to determine
	// the vertices of a triangle. Suppose you want to determine the vertices
	// of triangle i, with 0 <= i < gNumTriangles. Define:
	//
	// k0 = gIndexBuffer[3*i + 0]
	// k1 = gIndexBuffer[3*i + 1]
	// k2 = gIndexBuffer[3*i + 2]
	//
	// Now, the vertices of triangle i are at positions k0, k1, and k2 (in that
	// order) in the vertex array (which you should allocate yourself at line
	// 27).
	//
	// Note that this assumes 0-based indexing of arrays (as used in C/C++,
	// Java, etc.) If your language uses 1-based indexing, you will have to
	// add 1 to k0, k1, and k2.
}

void render()
{
	//Create our image. We don't want to do this in 
	//the main loop since this may be too slow and we 
	//want a responsive display of our beautiful image.
	//Instead we draw to another buffer and copy this to the 
	//framebuffer using glDrawPixels(...) every refresh

	int* gIndexBuffer;
	vec4* gVertexBuffer;
	int vertexNum = 0, triangleNum = 0;
	create_scene(gIndexBuffer, gVertexBuffer, &vertexNum, &triangleNum);
	mat4 modelScale(
		2, 0, 0, 0,
		0, 2, 0, 0,
		0, 0, 2, 0,
		0, 0, 0, 1
	);
	mat4 modeltranslate(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, -7, 1
	);
	mat4 modelTransform = modeltranslate * modelScale;
	vec3 e(0, 0, 0);
	vec3 u(1, 0, 0);
	vec3 v(0, 1, 0);
	vec3 w(0, 0, 1);
	
	mat4 cameraTransform(
		u.x, u.y, u.z, 0,
		v.x, v.y, v.z, 0,
		w.x, w.y, w.z, 0,
		e.x, e.y, e.z, 1
	);
	cameraTransform = inverse(cameraTransform);

	float l = -0.1;
	float r = 0.1;
	float b = -0.1;
	float t = 0.1;
	float n = -0.1;
	float f = -1000;
	
	mat4 p1(
		n, 0, 0, 0,
		0, n, 0, 0,
		0, 0, n + f, 1,
		0, 0, -f * n, 0
	);
	mat4 p2(
		2 / (r - l), 0, 0, 0,
		0, 2 / (t - b), 0, 0,
		0, 0, 2 / (n - f), 0,
		-(r + l) / (r - l), -(t + b) / (t - b), -(n + f) / (n - f), 1
	);
	mat4 perspectiveTransform = p2 * p1;


	int nx = 512;
	int ny = 512;

	mat4 viewportTransform(
		nx / 2, 0, 0, 0,
		0, ny / 2, 0, 0,
		0, 0, 1, 0,
		(nx - 1) / 2, (ny - 1) / 2, 0, 1
	);


	vec3 ka = { 0,1,0 };
	vec3 kd = { 0,0.5,0 };
	vec3 ks = { 0.5,0.5,0.5 };
	float p = 32;
	float la = 0.2;
	float intensity = 1;
	vec4 lightPos = { -4,4,-3, 1 };
	vec3 cameraPos = { 0,0,0 };

	//model Transform
	for (int i = 0; i < vertexNum; i++) {
		// Transform each point
		gVertexBuffer[i] = modelTransform * gVertexBuffer[i];
		gVertexBuffer[i] /= gVertexBuffer[i].w;
	}

	//illumination per triangle
	vec3* triangleColor = new vec3[triangleNum];
	for (int i = 0; i < triangleNum; i++) {
		vec4 t1 = gVertexBuffer[gIndexBuffer[3 * i + 0]];
		vec4 t2 = gVertexBuffer[gIndexBuffer[3 * i + 1]];
		vec4 t3 = gVertexBuffer[gIndexBuffer[3 * i + 2]];
		vec3 mid = (vec3)((t1 + t2 + t3) / 3.0f);
		vec3 v1 = vec3(t2 - t1);
		vec3 v2 = vec3(t3 - t1);

	    vec3 normal = normalize(cross(v1, v2));
		vec3 lay = (vec3)lightPos - mid;
		normal = normalize(normal);
		lay = normalize(lay);

		vec3 view = cameraPos - mid;
		view = normalize(view);
		vec3 h = view + lay;
		h = normalize(h);

		vec3 light = kd * la + kd * intensity * std::max(0.0f, dot(normal, lay)) + ks * intensity * pow(std::max(0.0f, dot(normal, h)),p);
		triangleColor[i] = light;
	}

	for (int i = 0; i < vertexNum; i++) {
		// Transform each point
		gVertexBuffer[i] = perspectiveTransform * cameraTransform * gVertexBuffer[i];
		gVertexBuffer[i] /= gVertexBuffer[i].w;
		gVertexBuffer[i] = viewportTransform * gVertexBuffer[i];
	}
	vector<vector<float>> depthBuffer(ny, vector<float>(nx, INF));
	vector<vector<vec3>> colorBuffer(ny, vector<vec3>(nx, {0,0,0}));
	

	for (int i = 0; i < triangleNum; i++) {
		vec4 t1 = gVertexBuffer[gIndexBuffer[3 * i + 0]];
		vec4 t2 = gVertexBuffer[gIndexBuffer[3 * i + 1]];
		vec4 t3 = gVertexBuffer[gIndexBuffer[3 * i + 2]];

		int xmin = floor(std::min({ t1.x, t2.x, t3.x }));
		int xmax = ceil(std::max({ t1.x, t2.x, t3.x }));
		int ymin = floor(std::min({ t1.y, t2.y, t3.y }));
		int ymax = ceil(std::max({ t1.y, t2.y, t3.y }));
		xmin = std::max(0, xmin);
		xmax = std::min(nx - 1, xmax);
		ymin = std::max(0, ymin);
		ymax = std::min(ny - 1, ymax);

		float b = ((t1.y - t3.y) * xmin + (t3.x - t1.x) * ymin + t1.x * t3.y - t3.x * t1.y) /
			((t1.y - t3.y) * t2.x + (t3.x - t1.x) * t2.y + t1.x * t3.y - t3.x * t1.y);

		float r = ((t1.y - t2.y) * xmin + (t2.x - t1.x) * ymin + t1.x * t2.y - t2.x * t1.y) /
			((t1.y - t2.y) * t3.x + (t2.x - t1.x) * t3.y + t1.x * t2.y - t2.x * t1.y);

		float bx = (t1.y - t3.y) / ((t1.y - t3.y) * t2.x + (t3.x - t1.x) * t2.y + t1.x * t3.y - t3.x * t1.y);
		float by = (t3.x - t1.x) / ((t1.y - t3.y) * t2.x + (t3.x - t1.x) * t2.y + t1.x * t3.y - t3.x * t1.y);
		float rx = (t1.y - t2.y) / ((t1.y - t2.y) * t3.x + (t2.x - t1.x) * t3.y + t1.x * t2.y - t2.x * t1.y);
		float ry = (t2.x - t1.x) / ((t1.y - t2.y) * t3.x + (t2.x - t1.x) * t3.y + t1.x * t2.y - t2.x * t1.y);

		for (int y = ymin; y <= ymax; y++) {
			for (int x = xmin; x <= xmax; x++) {
				if (b >= 0 && r >= 0 && b + r <= 1) {
					//mark triangle
					float a = 1 - b - r;
					float depth = -(t1.z * a + t2.z * b + t3.z * r);

					if (depthBuffer[y][x] > depth) {
						depthBuffer[y][x] = depth;
						colorBuffer[y][x] = triangleColor[i];
					}
				}
				b += bx;
				r += rx;
			}
			b -= (xmax - xmin + 1) * bx;
			r -= (xmax - xmin + 1) * rx;
			b += by;
			r += ry;
		}
	}

	float gamma = 2.2;

	OutputImage.clear();
	for (int j = 0; j < ny; ++j)
	{
		for (int i = 0; i < nx; ++i)
		{
			// ---------------------------------------------------
			// --- Implement your code here to generate the image
			// ---------------------------------------------------
			vec3 color = glm::vec3(0.0f, 0.0f, 0.0f);
			if (depthBuffer[j][i] != INF) { /// If it's not the initial value
				color = colorBuffer[j][i];
			}
			color.x = pow(color.x, 1/gamma);
			color.y = pow(color.y, 1/gamma);
			color.z = pow(color.z, 1/gamma);
			
			// set the color
			OutputImage.push_back(color.x); // R
			OutputImage.push_back(color.y); // G
			OutputImage.push_back(color.z); // B
		}
	}
}


void resize_callback(GLFWwindow*, int nw, int nh) 
{
	//This is called in response to the window resizing.
	//The new width and height are passed in so we make 
	//any necessary changes:
	Width = nw;
	Height = nh;
	//Tell the viewport to use all of our screen estate
	glViewport(0, 0, nw, nh);

	//This is not necessary, we're just working in 2d so
	//why not let our spaces reflect it?
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(0.0, static_cast<double>(Width)
		, 0.0, static_cast<double>(Height)
		, 1.0, -1.0);

	//Reserve memory for our render so that we don't do 
	//excessive allocations and render the image
	OutputImage.reserve(Width * Height * 3);
	render();
}


int main(int argc, char* argv[])
{
	// -------------------------------------------------
	// Initialize Window
	// -------------------------------------------------

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(Width, Height, "OpenGL Viewer", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	//We have an opengl context now. Everything from here on out 
	//is just managing our window or opengl directly.

	//Tell the opengl state machine we don't want it to make 
	//any assumptions about how pixels are aligned in memory 
	//during transfers between host and device (like glDrawPixels(...) )
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	//We call our resize function once to set everything up initially
	//after registering it as a callback with glfw
	glfwSetFramebufferSizeCallback(window, resize_callback);
	resize_callback(NULL, Width, Height);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		//Clear the screen
		glClear(GL_COLOR_BUFFER_BIT);

		// -------------------------------------------------------------
		//Rendering begins!
		glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
		//and ends.
		// -------------------------------------------------------------

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();

		//Close when the user hits 'q' or escape
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS
			|| glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		{
			glfwSetWindowShouldClose(window, GL_TRUE);
		}
	}

	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
