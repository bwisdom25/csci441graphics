#version 400
 
flat in float NdotL;   // interpolated NdotL values (output of the vertex shader!) 
// note that interpolation qualifiers have to match between vertex and fragment shader

out vec4 fragcolor; 
 
uniform float LightIntensity; 
uniform float AmbientIntensity; 
uniform vec3 DiffuseCoefficient;     // for RGB: this is why it's a 3D vector 
 
void main() { 
  fragcolor = vec4( (LightIntensity*(NdotL > 0.0 ? NdotL : 0.0) + AmbientIntensity) * DiffuseCoefficient, 1);  
      // note that some simplifying assumptions are made in the above formula 
      //  for example, k_a=k_d is used; also, specular term is not used; also, no attenuation here 
} 
