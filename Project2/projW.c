
/* 
PLEASE PUT YOUR NAME HERE
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <string>
#include <sstream>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glext.h>
#include <GL/glu.h>
#include <GL/glut.h>

using namespace std;



/*DATA STRUCTRES*/
class vertex{
public:
  double x,y,z;
};

class triangle{
public:
  int a,b,c;
}; 

/* --------------------------------------------- */
/* ----- GLOBAL VARIABLES  --------------------- */
/* --------------------------------------------- */

static const int VPD_DEFAULT = 800;

static const int MENU_SLOWER = 1;
static const int MENU_FASTER = 2;
static const int MENU_STOP_RUN = 3;

static const double TWOPI = (2.0 * M_PI);

int numTriangles,numVertices; 
triangle *triTable;
vertex *vertTable;



/* --------------------------------------------- */

GLint wid;               /* GLUT window id; value asigned in main() and should stay constant */
GLint vpw = VPD_DEFAULT; /* viewport dimensions; change when window is resized (resize callback) */
GLint vph = VPD_DEFAULT;

/* --------------------------------------------- */

GLfloat angle1 = 0;  /* angles used in animation */
GLfloat angle2 = 0;
GLfloat dangle1 = 0.0057;
GLfloat dangle2 = 0.0071;

bool animate = true;   /* animate or not? */

/* --------------------------------------------- */

GLuint ProgramHandle;   // this is an OpenGL name/ID number of the program object
// compiled in SetUpProgram

// handles of uniform variables - values assigned in SetUpProgram
GLint ModelViewMatrixHandle, ProjectionMatrixHandle, NormalMatrixHandle;
GLint LightLocationHandle, LightIntensityHandle, AmbientIntensityHandle;
GLint DiffuseCoefficientHandle;

/* --------------------------------------------- */
/* ----- SHADER RELATED CODE ------------------- */
/* --------------------------------------------- */

// read entire file into a string

string ReadFromFile ( const char *name )
{
  
    ifstream ifs(name);
  if (!ifs)
    {
      cout << "Can't open " << name << endl;
      exit(1);
    }
  stringstream sstr;
  sstr << ifs.rdbuf();
  return sstr.str();
 
}

void readInputFile(const char* name){
  ifstream ifs(name);
  if (!ifs)
    {
      cout << "Can't open " << name << endl;
      exit(1);
    }
   
  ifs >> numTriangles >> numVertices;
  
  triTable = new triangle[numTriangles];
  vertTable = new vertex[numVertices];
  
  for(int i=0;i < numTriangles; ++i) {
    ifs >> triTable[i].a >> triTable[i].b >> triTable[i].c;
  }
  
  for(int i=0;i < numVertices; ++i) {
    ifs >> vertTable[i].x >> vertTable[i].y >> vertTable[i].z; 
  }
  
  ifs.close();
}

/* --------------------------------------------- */

// function for printing compilation error/warning messages etc

void PrintInfoLog ( GLuint obj )
{
  int infologLength = 0;
  int maxLength;
 
  if(glIsShader(obj))
    glGetShaderiv(obj,GL_INFO_LOG_LENGTH,&maxLength);
  else
    glGetProgramiv(obj,GL_INFO_LOG_LENGTH,&maxLength);
 
  char infoLog[maxLength];
  
  if (glIsShader(obj))
    glGetShaderInfoLog(obj, maxLength, &infologLength, infoLog);
  else
    glGetProgramInfoLog(obj, maxLength, &infologLength, infoLog);
  
  if (infologLength > 0)
    cout << infoLog << endl;
}

/* --------------------------------------------- */

// set up program

void SetUpProgram()
{
  GLuint vid = glCreateShader(GL_VERTEX_SHADER);
  GLuint fid = glCreateShader(GL_FRAGMENT_SHADER);

  string VertexSrc = ReadFromFile("vertex.glsl");
  string FragmentSrc = ReadFromFile("fragment.glsl");

  const char *VS = VertexSrc.c_str();
  const char *FS = FragmentSrc.c_str();

  // set source of the shaders
  glShaderSource(vid,1,&VS,NULL);
  glShaderSource(fid,1,&FS,NULL);

  // compile
  glCompileShader(vid);
  glCompileShader(fid);

  // any problems?
  cout << "Vertex Shader Log: " << endl;
  PrintInfoLog(vid);
  cout << "Fragment Shader Log: " << endl;
  PrintInfoLog(fid);

  // create program
  ProgramHandle = glCreateProgram();

  // attach vertex and fragment shadets
  glAttachShader(ProgramHandle,vid);
  glAttachShader(ProgramHandle,fid);

  // link
  glLinkProgram(ProgramHandle);  

  // any problems?
  cout << "Linker Log: " << endl;
  PrintInfoLog(ProgramHandle);

  // get uniform variable locations so that they can be set from C code

  ModelViewMatrixHandle = glGetUniformLocation(ProgramHandle,"ModelViewMatrix");
  assert(ModelViewMatrixHandle!=-1);  // -1 means a problem

  ProjectionMatrixHandle = glGetUniformLocation(ProgramHandle,"ProjectionMatrix");
  assert(ProjectionMatrixHandle!=-1);

  NormalMatrixHandle = glGetUniformLocation(ProgramHandle,"NormalMatrix");
  assert(NormalMatrixHandle!=-1);

  LightLocationHandle = glGetUniformLocation(ProgramHandle,"LightLocation");
  assert(LightLocationHandle!=-1);

  LightIntensityHandle = glGetUniformLocation(ProgramHandle,"LightIntensity");
  assert(LightIntensityHandle!=-1);

  AmbientIntensityHandle = glGetUniformLocation(ProgramHandle,"AmbientIntensity");
  assert(AmbientIntensityHandle!=-1);

  DiffuseCoefficientHandle = glGetUniformLocation(ProgramHandle,"DiffuseCoefficient");
  assert(DiffuseCoefficientHandle!=-1);
}

/* --------------------------------------------- */
/* ---- RENDERING ROUTINES --------------------- */
/* --------------------------------------------- */

// draw a cube extending from -1 to 1 in model coordinate system

GLuint draw_cube ( )
{

  // This array contains normal and vertex data
  // generally, normals precede vertex coordinates
  // note that OpenGL 4 does not support quads any more

  static GLfloat array[] = 
    {
      0.0,0.0,-1.0,
      1.0,1.0,-1.0,
      0.0,0.0,-1.0,
      1.0,-1.0,-1.0,
      0.0,0.0,-1.0,
      -1.0,1.0,-1.0,

      0.0,0.0,-1.0,
      1.0,-1.0,-1.0,
      0.0,0.0,-1.0,
      -1.0,-1.0,-1.0,
      0.0,0.0,-1.0,
      -1.0,1.0,-1.0,

      0.0,0.0,1.0,
      1.0,1.0,1.0,
      0.0,0.0,1.0,
      -1.0,1.0,1.0,
      0.0,0.0,1.0,
      1.0,-1.0,1.0,

      0.0,0.0,1.0,
      -1.0,1.0,1.0,
      0.0,0.0,1.0,
      -1.0,-1.0,1.0,
      0.0,0.0,1.0,
      1.0,-1.0,1.0,

      1.0,0.0,0.0,
      1.0,-1.0,1.0,
      1.0,0.0,0.0,
      1.0,-1.0,-1.0,
      1.0,0.0,0.0,
      1.0,1.0,1.0,

      1.0,0.0,0.0,
      1.0,-1.0,-1.0,
      1.0,0.0,0.0,
      1.0,1.0,-1.0,
      1.0,0.0,0.0,
      1.0,1.0,1.0,

      -1.0,0.0,0.0,
      -1.0,-1.0,1.0,
      -1.0,0.0,0.0,
      -1.0,1.0,1.0,
      -1.0,0.0,0.0,
      -1.0,-1.0,-1.0,

      -1.0,0.0,0.0,
      -1.0,1.0,1.0,
      -1.0,0.0,0.0,
      -1.0,1.0,-1.0,
      -1.0,0.0,0.0,
      -1.0,-1.0,-1.0,

      0.0,1.0,0.0,
      -1.0,1.0,1.0,
      0.0,1.0,0.0,
      1.0,1.0,1.0,
      0.0,1.0,0.0,
      -1.0,1.0,-1.0,

      0.0,1.0,0.0,
      1.0,1.0,1.0,
      0.0,1.0,0.0,
      1.0,1.0,-1.0,
      0.0,1.0,0.0,
      -1.0,1.0,-1.0,

      0.0,-1.0,0.0,
      -1.0,-1.0,-1.0,
      0.0,-1.0,0.0,
      1.0,-1.0,-1.0,
      0.0,-1.0,0.0,
      -1.0,-1.0,1.0,

      0.0,-1.0,0.0,
      1.0,-1.0,-1.0,
      0.0,-1.0,0.0,
      1.0,-1.0,1.0,
      0.0,-1.0,0.0,
      -1.0,-1.0,1.0,    
    };

  // vertex attributes
  // the first argument matches the location qualifier of a vertex shader's input variable
  // here, location 0 means vertex coordinates and 1 - vertex normals
  glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,6*sizeof(GLfloat),&array[3]);
  glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,6*sizeof(GLfloat),array);

  // enable attributes # 0 and 1
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);

  // note that the last argument is the number of VERTICES sent into the pipeline
  // 12*3: 12 triangles, 3 vertices each
  glDrawArrays(GL_TRIANGLES,0,12*3);

  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
}

/* --------------------------------------------- */

// draw all cubes

void draw ( )
{
  if (animate)
    {
      // update angles
      angle1 += dangle1;
      angle2 += dangle2;
    }

  /* ensure we're drawing to the correct GLUT window */
  glutSetWindow(wid);

  // straightforward OpenGL settings
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glEnable(GL_DEPTH_TEST);

  /* clear the color buffers */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // use our program....
  // note that there may be several program handles in a more complex OpenGL code 
  //    - we have just one here though
  glUseProgram(ProgramHandle);
  
  // ----------- this part of code sets the uniform variable values

  // set uniform variable values
  glUniform1f(LightIntensityHandle,0.9);

  glUniform1f(AmbientIntensityHandle,0.1);
  glUniform3f(LightLocationHandle,200.0,200.0,200.0);

  // now the matrices... Arguments: field of view, aspect, front/back clip plane distances
  glm::mat4 PMat = glm::perspective(8.0f,1.0f,15.0f,25.0f);

  // send the matrix to shaders; note that &P[0][0] returns the pointer to matrix entries
  //   in an OpenGL-friendly format
  glUniformMatrix4fv(ProjectionMatrixHandle,1,GL_FALSE,&PMat[0][0]);

  // compute the modelview transformation matrix for the large cube
  glm::mat4 MVMat = glm::translate(glm::mat4(1.0),glm::vec3(0.0,0.0,-20.0)) *
    glm::rotate(glm::mat4(1.0),angle1,glm::vec3(1.0,2.0,3.0)) *
    glm::rotate(glm::mat4(1.0),angle2,glm::vec3(-2.0,-1.0,0.0));
    
  // pass modelview matrix to the shaders...
  glUniformMatrix4fv(ModelViewMatrixHandle,1,GL_FALSE,&MVMat[0][0]);

  // normal matrix should be the rotational component of the modelview matrix
  // note that normal matrix is supposed to be 3x3 (see vertex shader code)

  // select the 3x3 block from MV
  glm::mat3 NMat(MVMat);

  // send normal matrix to shader
  glUniformMatrix3fv(NormalMatrixHandle,1,GL_FALSE,&NMat[0][0]);

  // pass color of the large cube to our shaders....
  glUniform3f(DiffuseCoefficientHandle,1.0,1.0,1.0);

  // ----------- draw the first cube

  draw_cube();

  // change colors...
  glUniform3f(DiffuseCoefficientHandle,1.0,0.0,0.0);

  // now, we need to change the modelview matrix so that the cube is scaled and moved to
  //    the corner of the first one...
  // note the order: first scale, then translate finally apply the previous transformation
  glm::mat4 CurrentMVMat = MVMat * 
    glm::translate(glm::mat4(1.0),glm::vec3(-1.0,-1.0,-1.0)) * glm::scale(glm::mat4(1.0),glm::vec3(.4,.4,.4));

  glUniformMatrix4fv(ModelViewMatrixHandle,1,GL_FALSE,&CurrentMVMat[0][0]);  

  draw_cube();

  // now draw the remaining 4 cubes the same way...

  glUniform3f(DiffuseCoefficientHandle,0.0,1.0,0.0);
  CurrentMVMat = MVMat * 
    glm::translate(glm::mat4(1.0),glm::vec3(-1.0,-1.0,1.0)) * glm::scale(glm::mat4(1.0),glm::vec3(.4,.4,.4));
  glUniformMatrix4fv(ModelViewMatrixHandle,1,GL_FALSE,&CurrentMVMat[0][0]);  
  draw_cube();

  glUniform3f(DiffuseCoefficientHandle,0.0,0.0,1.0);
  CurrentMVMat = MVMat * 
    glm::translate(glm::mat4(1.0),glm::vec3(-1.0,1.0,-1.0)) * glm::scale(glm::mat4(1.0),glm::vec3(.4,.4,.4));
  glUniformMatrix4fv(ModelViewMatrixHandle,1,GL_FALSE,&CurrentMVMat[0][0]);  
  draw_cube();
  
  glUniform3f(DiffuseCoefficientHandle,0.5,0.0,0.5);
  CurrentMVMat = MVMat * 
    glm::translate(glm::mat4(1.0),glm::vec3(1.0,-1.0,-1.0)) * glm::scale(glm::mat4(1.0),glm::vec3(.4,.4,.4));
  glUniformMatrix4fv(ModelViewMatrixHandle,1,GL_FALSE,&CurrentMVMat[0][0]);  
  draw_cube();


  /* flush the pipeline */
  glFlush();

  /* look at our handiwork */
  glutSwapBuffers();

  if (animate)
    glutPostRedisplay();  
}

/* --------------------------------------------- */
/* --- INTERACTION ----------------------------- */
/* --------------------------------------------- */

/* handle mouse events */

void mouse_button(GLint btn, GLint state, GLint mx, GLint my)
{
  switch( btn ) {
    case GLUT_LEFT_BUTTON:
      cout << "Left Button"; // remove this line from your final submission
      switch( state ) {
        case GLUT_DOWN: 
	  cout << " down" << endl; // remove this line from your final submission
	  break;
        case GLUT_UP:   
	  cout << " up" << endl; // remove this line from your final submission
	  break;
      }
      break;
    case GLUT_MIDDLE_BUTTON:
      cout << "Middle Button"; // remove this line from your final submission
      switch( state ) {
        case GLUT_DOWN: 
	  cout << " down" << endl; // remove this line from your final submission
	  break;
        case GLUT_UP:   
	  cout << " up" << endl; // remove this line from your final submission
	  break;
      }
      break;
    case GLUT_RIGHT_BUTTON:
      cout << "Right Button"; // remove this line from your final submission
      switch( state ) {
        case GLUT_DOWN: 
	  cout << " down" << endl; // remove this line from your final submission
	  break;
        case GLUT_UP:   
	  cout << " up" << endl; // remove this line from your final submission
	  break;
      }
      break;
  }
}

/* --------------------------------------------- */

// mouse moves with button down

GLvoid button_motion(GLint mx, GLint my)
{
  cout << "Motion with button down: " << mx << "," << my << endl; // remove this line from your final submission
  return;
}

/* --------------------------------------------- */

// mouse moves with button up

GLvoid passive_motion(GLint mx, GLint my)
{
  cout << "Passive Motion: " << mx << "," << my << endl; // remove this line from your final submission
  return;
}

/* --------------------------------------------- */

/* handle keyboard events; here, just exit if ESC is hit */

void keyboard(GLubyte key, GLint x, GLint y)
{
  switch(key) {
    case 27:  /* ESC */
              exit(0);

    default:  break;
  }
}

/* --------------------------------------------- */

/* menu callback */

void menu ( int value )
{
  switch(value)
    {
    case MENU_SLOWER:
      dangle1 *= .5;
      dangle2 *= .5;
      break;
    case MENU_FASTER:
      dangle1 *= 1.5;
      dangle2 *= 1.5;
      break;
    case MENU_STOP_RUN:
      animate = !animate;
      break;
    }
}

/* --------------------------------------------- */

/* handle resizing the glut window */

GLvoid reshape(GLint sizex, GLint sizey)
{
  glutSetWindow(wid);

  vpw = sizex;
  vph = sizey;

  glViewport(0, 0, vpw, vph);
  glutReshapeWindow(vpw, vph);

  glutPostRedisplay();
}

/* --------------------------------------------- */
/* -------- SET UP GLUT  ----------------------- */
/* --------------------------------------------- */

// initialize glut, callbacks etc.
GLint init_glut(GLint *argc, char **argv)
{
  GLint id;

  glutInit(argc,argv);

  /* size and placement hints to the window system */
  glutInitWindowSize(vpw, vph);
  glutInitWindowPosition(10,10);

  /* double buffered, RGB color mode */
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

  /* create a GLUT window (not drawn until glutMainLoop() is entered) */
  id = glutCreateWindow("MACS441 OpenGL Sample code");    

  /* register callbacks */

  /* window size changes */
  glutReshapeFunc(reshape);

  /* keypress handling when the current window has input focus */
  glutKeyboardFunc(keyboard);

  /* mouse event handling */
  glutMouseFunc(mouse_button);           /* button press/release        */
  glutMotionFunc(button_motion);         /* mouse motion w/ button down */
  glutPassiveMotionFunc(passive_motion); /* mouse motion with button up */

  /* window obscured/revealed event handler */
  glutVisibilityFunc(NULL);

  /* handling of keyboard SHIFT, ALT, CTRL keys */
  glutSpecialFunc(NULL);

  /* what to do when mouse cursor enters/exits the current window */
  glutEntryFunc(NULL);

  /* what to do on each display loop iteration */
  glutDisplayFunc(draw);

  /* create menu */
  GLint menuID = glutCreateMenu(menu);
  glutAddMenuEntry("slower",MENU_SLOWER);
  glutAddMenuEntry("faster",MENU_FASTER);
  glutAddMenuEntry("stop/run",MENU_STOP_RUN);
  glutSetMenu(menuID);
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  return id;
}

/* --------------------------------------------- */
/* --------------------------------------------- */
/* --------------------------------------------- */


GLint main(GLint argc, char **argv)
{
  //CHECK ARGC 
  
  if(argc != 2 ){
	cout << "ERROR: NUMBER OF ARGUMENTS"<<endl << "Correct Use: ./projW input_file.txt "<<endl; 
	return 69;
  }
  
  
  readInputFile(argv[1]);
  /* initialize GLUT: register callbacks, etc */
  
  wid = init_glut(&argc, argv);

  // initialize glew and check for OpenGL 4.0 support
  glewInit();
  if (glewIsSupported("GL_VERSION_4_0"))
    cout << "Ready for OpenGL 4.0" << endl;
  else 
    {
      cout << "OpenGL 4.0 not supported" << endl;;
      return 1;
    }
  
  // set up vertex/fragment programs
  SetUpProgram();
  
  // main loop: keep processing events until exit is requested
  glutMainLoop();
  
  return 0;
}


/* --------------------------------------------- */
