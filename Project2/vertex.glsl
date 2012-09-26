#version 400

layout(location=0) in vec3 coord;   // in model coordinate system 
layout(location=1) in vec3 normal;  // also in model coordinate system 
 
flat out float NdotL; // flat is the interpolation qualifier: request flat interpolation 
// other interpolation qualifiers: 
//      nonperspective = linear in screen space (this is the way in which color is interpolated in Gouraud shading)
//      smooth = perspectively correct  
 
uniform mat4 ModelViewMatrix; 
uniform mat4 ProjectionMatrix; 
uniform mat3 NormalMatrix; 
uniform vec3 LightLocation; 
 
void main() {  
  vec4 WorldCoord = ModelViewMatrix * vec4(coord,1.0);   // convert to world coordinates  
  vec3 L = normalize(LightLocation - WorldCoord.xyz);    // L vector for illumination 
  vec3 WorldNormal = NormalMatrix*normal;                // normal in world coordinates 
  vec3 N = normalize(WorldNormal);                       // N vector for illumination 
  NdotL = dot(N,L);                                      // part of diffuse term (multiplied by k_d's etc in the fragment shader) 
  gl_Position = ProjectionMatrix*WorldCoord;             // gl_Position is a predefined variable 
    // a correctly written vertex shader should write screen space coordinates to gl_Position 
    // they are used on the rasterization stage! 
} 
