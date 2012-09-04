
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


/*vector class*/
class vector{
	public:
		double x,y,z;
		double dot(const vector b);
		double dot(const vector a, const vector b); 
		vector cross(const vector b); 
		vector cross(const vector a, const vector b); 
		void setVec(const vector b);  
		double getMag(); 		
};

double vector::dot( const vector b)
{
	return((x*b.x)+(y*b.y)+(z*b.z));
}

double vector::dot(const vector a, const vector b)
{
	return((a.x*b.x)+(a.y*b.y)+(a.z*b.z));
}
/*-------------------------------*/ 

vector vector::cross( const vector b){
	vector temp; 
	temp.x = (y*b.z) - (z*b.y);
	temp.y = (z*b.x) - (x*b.z);
	temp.z = (x*b.y) - (y*b.x); 
	return temp;
}

vector vector::cross(const vector a, const vector b){
	vector temp; 
	temp.x = (a.y*b.z) - (a.z*b.y);
	temp.y = (a.z*b.x) - (a.x*b.z);
	temp.z = (a.x*b.y) - (a.y*b.x); 
	return temp;
}

void vector::setVec(const vector b){
	x=b.x;
	y=b.y;
	z=b.z;
}
double vector::getMag()
{
	return( sqrt( (x*x)+(y*y)+(z*z) ) );
}


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
		vector origin; 
		vector direction; 
}; 

/*sphere class*/ 
class sphere{
	public:
		vector center;
		double radius;
		material m; 
}; 

/*triangle class*/ 
class triangle{
	public: 
		vector a1,a2,a3;
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

ofstream outfile("output.txt");

double intersectionS(ray r,sphere s){
//Define coefficients for Quadtratic eqn
	double A,B,C,Det,t1,t2; 
	vector co; 

	co = r.origin - s.center; // o-c GIVES vector co 
 
	//Calculate Coefficients
	C = pow(co.getMag(),2) - pow(s.radius,2); 
	B = 2*r.direction.dot(co);
	A = pow(r.direction.getMag(),2); 
   
	//Caclulate the Determinant 
	Det = pow(B,2) - (4*A*C); 

	//Take Appropriate Actions (reguarding t-value) 
	if( Det < 0 ){ //Negative Det 
		return -1.0; //ANY Negative val. corresponds to NO Solution 
    }else{ 
		t1=(-B+sqrt(Det))/(2*A);
		t2=(-B-sqrt(Det))/(2*A);
		if(Det == 0){
			
			if(t2 > 0){return t2;}
			if(t1 > 0){return t1;}
			else{
				return -1.0; //line is tangent, but not on ray
			}
		}else if(Det > 0){
				 
				if(t2 > 0 && t1 < 0 ){ return t2;}
				if(t1 > 0 && t2 < 0 ){ return t1;}
				if(t2 < t1){ return t2; } 
				if(t1 < t2){ return t1; }				
				else{
					return -1.0; //line is tangent, but not on ray
				}
		}
    }	
	
}
//TODO: dse a strucouble intersectionT(ray r,triangle t) 

double intersectionT(ray r,triangle t){

	vector a1a2,a1a3,n,temp1,temp2; 
	vector v1,v2,v3; 
	vector p;

	//Determine if point p is in triangles plane 
	double t_val,returnValue=-1.0; 
	a1a2 = t.a2 - t.a1; 
	a1a3 = t.a3 - t.a1; 
	
	n = a1a2.cross(a1a3); 
	temp1 = t.a1 - r.origin;
	t_val=(temp1.dot(n))/(r.direction.dot(n)); 
	if(t_val<0){returnValue = -1.0;}
	
		p = r.origin + (t_val*r.direction);
        
		vector pa1 = t.a1 - p;
        vector pa2 = t.a2 - p;
        vector pa3 = t.a3 - p;

        v1 = v1.cross(pa1, pa2);
        v2 = v2.cross(pa2, pa3);
        v3 = v3.cross(pa3, pa1);
		


		if(v1.dot(v2) >= 0 ){
			if(v2.dot(v3) >= 0){
				if(v3.dot(v1) >= 0){
					returnValue = t_val; 
				}
			}
		}

		if(v1.dot(v2) <= 0 ){
			if(v2.dot(v3) <= 0){
				if(v3.dot(v1) <= 0){
					returnValue = t_val; 
				}
			}
		}
	
		
		return returnValue;

}



/* ----------- Reading the input file ---------- */

// global variables
int resolution_x, resolution_y;
triangle *T;
sphere *S;
int n_T,n_S;
vector viewpoint,lowerleft; 
vector horz,vert; 


double intersection(ray r, int pid){
	if(pid < n_T){
		return intersectionT(r,T[pid]); 
	}else{
		return intersectionS(r,S[pid-n_T]);
	}
}

intersectPrim closestIntersect(ray r){
	intersectPrim temp; 
	double mint=-1.0; 
	double t; 
	int pid;
	for(int i=0;i<(n_S+n_T);++i){
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
	l = lowerleft;
	e = viewpoint;
    tempA=(((i+0.5)/resolution_x)*horz);
	tempB=(((j+0.5)/resolution_y)*vert);
	tempC=l+tempA+tempB;
	direct=tempC-e; 
	temp.direction = direct; 
	
	return temp;
}


	

// ... and the input file reading function
void read_input_file()
{
  ifstream ifs("input9r.txt");
  
  
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
	sphere temp_s; 
	triangle temp_t;
  for ( int i=0; i<number_of_primitives; i++ )
    {
      char primitive_type;
      ifs >> primitive_type;
      switch(primitive_type)
	{
	case 's':
	case 'S':
	  {
		
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

		S[i-n_T]=temp_s; 
		//outfile << i << " " << temp_s.center.x <<" "<<temp_s.center.y<<" "<<temp_s.center.z << " "<<temp_s.radius<<endl;
		++n_S; 
	    // add the sphere to your datastructures (primitive list, sphere list or such) here
	  }
	  break;
	case 'T':
	case 't':
	  {
	    ifs >> temp_t.a1.x >> temp_t.a1.y >> temp_t.a1.z;
	    ifs >> temp_t.a2.x >> temp_t.a2.y >> temp_t.a2.z;
	    ifs >> temp_t.a3.x >> temp_t.a3.y >> temp_t.a3.z;
	    ifs >> temp_t.m.k_diff_r >> temp_t.m.k_diff_g >> temp_t.m.k_diff_b;
	    ifs >> temp_t.m.k_amb_r >> temp_t.m.k_amb_g >> temp_t.m.k_amb_b;
	    ifs >> temp_t.m.k_spec >> temp_t.m.n_spec ; 	    

	    // add the triangle to your datastructure (primitive list, sphere list or such) here
		T[i-n_S]=temp_t;
		//outfile << i << " " << temp_t.a1.x <<" "<<temp_t.a1.y <<" "<<temp_t.a1.z  <<endl;
		++n_T; 
		
	  }
	  break;
	default:
	  assert(0);
	}
    }
	//outfile << "NUM_OF_TRI: "<< n_T << "  NUM_OF_SPH: "<< n_S <<endl; 
}

/* ----------- main function ---------- */

int main ( int argc, char *argv[] )
{
  int x,y;
  intersectPrim prim; 
   
  read_input_file();

  //FOR DEBUGGIN PURPOSES ONLY 


  image img(resolution_x,resolution_y);
  for ( x=0; x<resolution_x; x++ )
    for ( y=0; y<resolution_y; y++ )
      {
		RGB &pix = img.pixel(x,y);

		ray r=eyeRay(x,y);
		prim=closestIntersect(r);
		if(prim.t == -1.0){
			pix.r=0.0;
			pix.g=0.0;
			pix.b=0.0;
		}else{	
			if(prim.pid < n_T){		//outfile << S[prim.pid].m.k_diff_r << S[prim.pid].m.k_diff_g << S[prim.pid].m.k_diff_b;
				pix.r=T[prim.pid].m.k_diff_r; 
				pix.g=T[prim.pid].m.k_diff_g;
				pix.b=T[prim.pid].m.k_diff_b;
			}else{
				pix.r=S[prim.pid-n_T].m.k_diff_r; 
				pix.g=S[prim.pid-n_T].m.k_diff_g;
				pix.b=S[prim.pid-n_T].m.k_diff_b;
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
   outfile.close();
  return 0;
}
