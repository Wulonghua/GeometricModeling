#ifdef GL_ES
// Set default precision to medium
precision mediump int;
precision mediump float;
#endif

attribute vec3 aVertexPosition;
attribute vec3 aVertexNormal;

uniform mat4 uMVPMatrix;
uniform mat4 uMVMatrix;
uniform mat4 uNMatrix;
uniform vec4 uVertexColor;

varying vec4 pos_in_eye;  //vertex position in eye space
varying vec3 v_normal;  // vertex normal
varying vec3 vColor;


void main(void) {

// transform normal from local to eye space: normal matrix is the inverse transpose of the modelview matrix
v_normal =vec3(uNMatrix*vec4(aVertexNormal,0.0));

// transform the vertex position to eye space
pos_in_eye = uMVMatrix*vec4(aVertexPosition, 1.0);

gl_Position = uMVPMatrix*vec4(aVertexPosition, 1.0);

vColor = vec3(uVertexColor);
}
