
/*
put your name here, just in case ****
Project 1: Ray Tracing
CSCI 441
*/

#include <math.h>
#include <fstream>
#include <iostream>
#include <assert.h>

using namespace std;

/* ------------ typedefs -------------- */

typedef struct {
	double r,g,b;
}
RGB;

typedef struct {
	double t;
	int pid; 
}
intersectPrim; 

/* ----------------------- */
/*point class*/
class point{
	public:
		double x,y,z;
};
/*vector class*/
class vector{
	public:
		double x,y,z;
		double dot(const vector b);
		void setVec(const point a, const point b);
		void setVec(const vector b);  
		double getMag(); 
		void toVec(const point a); 
};

vector operator+(const vector a, const vector b)
{
	vector result; 
	result.x = a.x + b.x;
	result.y = a.y + b.y; 
	result.z = a.z + b.z; 

	return result; 
}

vector operator-(const vector a, const vector b)
{
	vector result; 
	result.x = a.x - b.x;
	result.y = a.y - b.y; 
	result.z = a.z - b.z; 

	return result; 
}
vector operator*(const double d, const vector a) 
{
	vector result; 
	result.x = d*a.x;
	result.y = d*a.y; 
	result.z = d*a.z; 

	return result; 
}


double vector::dot( const vector b)
{
	return((x*b.x)+(y*b.y)+(z*b.z));
}
/*-------------------------------*/ 



void vector::setVec(const point a, const point b)
{
	x=a.x-b.x; 
	y=a.y-b.y;
	z=a.z-b.z;
}

void vector::setVec(const vector b){
	x=b.x;
	y=b.y;
	z=b.z;
}
double vector::getMag()
{
	return( sqrt((pow(x,2))+(pow(y,2))+(pow(z,2))) );
}

point operator+(const point a, const point b) 
{
	point result; 
	result.x = a.x + b.x;
	result.y = a.y + b.y; 
	result.z = a.z + b.z; 

	return result;
}

void vector::toVec(const point a){
	x=a.x;
	y=a.y;
	z=a.z;
}

point operator-(const point a, const point b) 
{
	point result; 
	result.x = a.x + b.x;
	result.y = a.y + b.y; 
	result.z = a.z + b.z; 

	return result;
} 
/*--------------------------------------------*/
/*material class*/ 
class material{
	public:
		double k_diff_r,k_diff_g,k_diff_b; 
		double k_amb_r,k_amb_g,k_amb_b; 
		double k_spec; 
		double n_spec;
}; 

/*ray class*/ 
class ray{ //charles 
	public:
		point origin; 
		vector direction; 
}; 

/*sphere class*/ 
class sphere{
	public:
		point center;
		double radius;
		material m; 
}; 

/*triangle class*/ 
class triangle{
	public: 
		point a1,a2,a3;
		material m;
};
		

static unsigned char clampnround ( double x )
{
  if (x>255)
    x = 255;
  if (x<0) 
    x = 0;
  return (unsigned char)floor(x+.5);
}

/* ----------------------- */

/*image class*/ 
class image {
	int xsize,ysize; // resolution
  	RGB *rgb;        // pixel intensities
 public:
	image ( int m, int n );       // allocates image of specified size
	RGB &pixel ( int i, int j );  // access to a specific pixel
	void save_to_ppm_file ( char *filename );
};

image::image ( int m, int n ) : xsize(m), ysize(n)
{
	rgb = new RGB[m*n];
}

RGB &image::pixel ( int i, int j )
{
	return rgb[i+xsize*j];
}

void image::save_to_ppm_file ( char *filename )
{
  ofstream ofs(filename,ios::binary);
  assert(ofs);
  ofs << "P6" << endl;
  ofs << xsize << " " << ysize << endl << 255 << endl;
  for ( int i=0; i<xsize*ysize; i++ )
    {
      unsigned char r = clampnround(256*rgb[i].r);
      unsigned char g = clampnround(256*rgb[i].g);
      unsigned char b = clampnround(256*rgb[i].b);
      ofs.write((char*)&r,sizeof(char));
      ofs.write((char*)&g,sizeof(char));
      ofs.write((char*)&b,sizeof(char));
    }
}
/*primative intersection functions*/ 

double intersectionS(ray r,sphere s){
//Define coefficients for Quadtratic eqn
	double A,B,C,Det,t1,t2; 
	vector co; 
	co.setVec(s.center,r.origin); 

	A = pow(co.getMag(),2) - pow(s.radius,2); 
	B = r.direction.dot(2*co);
	C = pow(r.direction.getMag(),2); 

	//Caclulate the Determinant 
	Det = pow(B,2) - (4*A*C); 

	//Take Appropriate Actions (reguarding t-value) 
	if( Det < 0 ){ //Negative Det 
		return -1.0; //ANY Negative val. corresponds to NO Solution 
    }else{ 
		t1=(-B+sqrt(Det))/(2*A);
		t2=(-B-sqrt(Det))/(2*A);
		if(t1 == 0){
			if(t2 > 0){
				return t2;
			}else{
				return -1.0; //line is tangent, but not on ray
			}
		}else if(t2 == 0){
				if(t2 > 0){
					return t2;
				}else{
					return -1.0; //line is tangent, but not on ray
				}
		}else if(t1 < 0 && t2 < 0){
			return -1.0;//no intersection 
		}else if(t1 > 0 && t2<0){ return t1; }
		 else if(t2 > 0 && t1<0){ return t2; } 
         else{
			if(t1 < t2){ return t1; } 
			else if(t2 < t1){ return t2; }
         }
    }	
	
}
//TODO: dse a strucouble intersectionT(ray r,triangle t) 



/* ----------- Reading the input file ---------- */

// global variables
int resolution_x, resolution_y;
triangle *T;
sphere *S;
int n_T,n_S;
point viewpoint,lowerleft; 
vector horz,vert; 


double intersection(ray r, int pid){
	if(pid < n_S){
		return intersectionS(r,S[pid]); 
	}
}

intersectPrim closestIntersect(ray r){
	intersectPrim temp; 
	double mint=-1.0; 
	double t; 
	int pid;
	for(int i=0;i<n_S-1;++i){
		t=intersection(r,i);
		if( t!=-1.0  && ( mint==-1.0 || t<mint )){
			mint=t;		
			pid=i;
 		}
	}
	temp.t=mint;
	temp.pid=pid; 
	return temp;
}

ray eyeRay( int i, int j){
	ray temp;
	vector direct,tempA,tempB,tempC; 

	temp.origin.x = viewpoint.x;
	temp.origin.y = viewpoint.y;
    temp.origin.z = viewpoint.z;
	
    vector l,e; 
	l.toVec(lowerleft);
	e.toVec(viewpoint);
    tempA=(((i*0.5)/resolution_x)*horz);
	tempB=(((i*0.5)/resolution_y)*vert);
	tempC=l+tempA+tempB;
	direct=tempC-e; 
	temp.direction = direct; 
	
	return temp;
}


	

// ... and the input file reading function
void read_input_file()
{
  ifstream ifs("input.txt");
  assert(ifs);
  /*Existing data types
  double viewpoint[3];
  double screen_lower_left_corner[3];
  double screen_horizontal_vector[3];
  double screen_vertical_vector[3];
  */
  double light_source[3];
  double light_intensity;
  double ambient_light_intensity;
  int number_of_primitives;

  ifs >> resolution_x >> resolution_y;
  ifs >> viewpoint.x >> viewpoint.y >> viewpoint.z;
  ifs >> lowerleft.x >> lowerleft.y >> lowerleft.z; 
  ifs >> horz.x >> horz.y >> horz.z;
  ifs >> vert.x >> vert.y >> vert.z; 

  /*Existing code
  ifs >> viewpoint[0] >> viewpoint[1] >> viewpoint[2];
  ifs >> screen_lower_left_corner[0] >> screen_lower_left_corner[1] >> screen_lower_left_corner[2];
  ifs >> screen_horizontal_vector[0] >> screen_horizontal_vector[1] >> screen_horizontal_vector[2];
  ifs >> screen_vertical_vector[0] >> screen_vertical_vector[1] >> screen_vertical_vector[2];
  */
  ifs >> light_source[0] >> light_source[1] >> light_source[2];
  ifs >> light_intensity;
  ifs >> ambient_light_intensity;
  ifs >> number_of_primitives;

  //allocate memory to hold primitives(wasting memory but oh well) 
  n_T = 0; 
  n_S = 0; 
  T = new triangle[number_of_primitives];
  S = new sphere[number_of_primitives];
  // save all this info to your datastructures or global variables here

  for ( int i=0; i<number_of_primitives; i++ )
    {
      char primitive_type;
      ifs >> primitive_type;
      switch(primitive_type)
	{
	case 's':
	case 'S':
	  {
		sphere temp_s; 
        /*EXISTING VARS
	    double center[3];
	    double radius;
	    double k_diffuse[3];
	    double k_ambient[3];
	    double k_specular;
	    double n_specular;
		*/
		ifs >> temp_s.center.x >> temp_s.center.y >> temp_s.center.z; 
		ifs >> temp_s.radius; 
        ifs >> temp_s.m.k_diff_r >> temp_s.m.k_diff_g >> temp_s.m.k_diff_b;
		ifs >> temp_s.m.k_amb_r >> temp_s.m.k_amb_g >> temp_s.m.k_amb_b; 
		ifs >> temp_s.m.k_spec >> temp_s.m.n_spec;
 
        /*EXISTING CODE
	    ifs >> center[0] >> center[1] >> center[2];
	    ifs >> radius;
	    ifs >> k_diffuse[0] >> k_diffuse[1] >> k_diffuse[2];
	    ifs >> k_ambient[0] >> k_ambient[1] >> k_ambient[2];
	    ifs >> k_specular >> n_specular;
        */
		S[i]=temp_s; 
		++n_S; 
	    // add the sphere to your datastructures (primitive list, sphere list or such) here
	  }
	  break;
	case 'T':
	case 't':
	  {
		//TODO: Read in triangle attributes, add triangle object to array 
	    double a1[3];
	    double a2[3];
	    double a3[3];
	    double k_diffuse[3];
	    double k_ambient[3];
	    double k_specular;
	    double n_specular;
	    
	    ifs >> a1[0] >> a1[1] >> a1[2];
	    ifs >> a2[0] >> a2[1] >> a2[2];
	    ifs >> a3[0] >> a3[1] >> a3[2];
	    ifs >> k_diffuse[0] >> k_diffuse[1] >> k_diffuse[2];
	    ifs >> k_ambient[0] >> k_ambient[1] >> k_ambient[2];
	    ifs >> k_specular >> n_specular; 	    

	    // add the triangle to your datastructure (primitive list, sphere list or such) here
		++n_T; 
	  }
	  break;
	default:
	  assert(0);
	}
    }
}

/* ----------- main function ---------- */

int main ( int argc, char *argv[] )
{
  int x,y;
  intersectPrim prim; 

  read_input_file();

  image img(resolution_x,resolution_y);
  for ( x=0; x<resolution_x; x++ )
    for ( y=0; y<resolution_y; y++ )
      {
		RGB &pix = img.pixel(x,y);

		ray r=eyeRay(x,y);
		closestIntersect(r);
		if(prim.t == -1.0){
			pix.r=0.0;
			pix.g=0.0;
			pix.b=0.0;
		}else{
				if( prim.pid < n_T ){
					//return rgb for pid of T 
   				 }else{
					pix.r=S[prim.pid].m.k_diff_r; 
					pix.g=S[prim.pid].m.k_diff_g;
					pix.b=S[prim.pid].m.k_diff_b;
				}
		}
	
	/* 
	   call your raytracer here
	   then assign the rgb values
	   it returns to the pixel 
	*/

	// this is just to produce a fun image...
    /*
	pix.r = 0.5+0.5*sin(sin(x/30.0)+y*y/700.0);
	pix.g = 0.5+0.5*sin(y/71.0);
	pix.b = 0.5+0.5*sin(x*x*x/120000.0+y*y/1700.0);
    */
      }

  /* save the image */
  img.save_to_ppm_file((char*)"output.ppm");

  return 0;
}
