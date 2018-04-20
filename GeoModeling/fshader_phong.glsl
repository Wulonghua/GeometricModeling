#ifdef GL_ES
// Set default precision to medium
precision mediump int;
precision mediump float;
#endif

uniform mat4 uMVPMatrix;
uniform mat4 uMVMatrix;
uniform mat4 uNMatrix;

uniform vec4 light_pos;
uniform vec4 ambient_coef;
uniform vec4 diffuse_coef;
uniform vec4 specular_coef;

uniform vec4 light_ambient;
uniform vec4 light_diffuse;
uniform vec4 light_specular;

varying vec4 pos_in_eye;
varying vec3 v_normal;
varying vec4 vColor;

void main(void)
{ 
	vec4 light_pos_in_eye = light_pos;
	vec3 vnormal = normalize(v_normal);
    vec3 light_vector = normalize(vec3(light_pos_in_eye - pos_in_eye));

    //float ndotl = max(dot(vnormal, light_vector),0.0);
	float ndotl = dot(vnormal, light_vector);
	vec4 ambient = ambient_coef * light_ambient;
	vec4 diffuse = diffuse_coef * light_diffuse* ndotl;


	vec3 R = normalize(vec3(reflect(-light_vector, v_normal)));
	vec3 eye_vector = normalize(-vec3(pos_in_eye));
	float rdotv = max(dot(R, eye_vector), 0.0);

	vec4 specular;
	if (ndotl>0)
		specular = specular_coef* light_specular*pow(rdotv, 10);
	else{
		//specular = vec4(0, 0, 0, 1);
		specular = specular_coef* light_specular*pow(rdotv, 10);
		//vColor   = vec4(0.25,0.9,0.4,1);
		diffuse  = -diffuse;
	}
	vec4 tmp = clamp(vColor *(ambient + diffuse)+specular,0.0,1.0);
	tmp[3] = 0.3;
	gl_FragColor = tmp;

}

