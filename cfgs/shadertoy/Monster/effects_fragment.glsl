#extension GL_ARB_shader_texture_lod : enable
//#extension GL_ARB_gpu_shader_fp64 : enable

// FXAA shader, GLSL code adapted from:
// http://horde3d.org/wiki/index.php5?title=Shading_Technique_-_FXAA
// Whitepaper describing the technique:
// http://developer.download.nvidia.com/assets/gamedev/files/sdk/11/FXAA_WhitePaper.pdf
// alpha depth mipmap DoF added by marius

uniform sampler2D my_texture;

uniform float dof_scale;  // {min=-29.5. max=1000. step=.1}
uniform float dof_offset; // {min=-5. max=5. step=.01}

uniform float xres, yres;
uniform float speed;

varying vec2 texture_coordinate;

// Fetch texel from mipmap, depth clue in alpha layer 0.
vec4 texture2Dx(sampler2D text, vec2 coord) {
  float z = texture2D(text, coord).w;
  float z_near = abs(speed);
  float z_far = 65535.0 * z_near;
  float dist = -z_far * z_near / (z * (z_far - z_near) - z_far);
  float d = dist / (z_near * (30.0 + dof_scale));
  float lod = abs(log(1.0 + float(d)) - .5 + dof_offset);
  return vec4(texture2DLod(text, coord, lod).xyz, lod);
}

void main() {
  // The inverse of the texture dimensions along X and Y
  vec2 texcoordOffset = vec2(1.0/xres, 1.0/yres);

  // The parameters are hardcoded for now, but could be
  // made into uniforms to control from the program.
  float FXAA_SPAN_MAX = 8.0;
  float FXAA_REDUCE_MUL = 1.0/8.0;
  float FXAA_REDUCE_MIN = (1.0/128.0);

  float z = texture2D(my_texture, texture_coordinate.xy).w;  // save z
  if (z == float(0)) {
    // 0 is used for background hit, i.e. missed fractal. Do not touch.
    gl_FragColor = texture2D(my_texture, texture_coordinate.xy);
    gl_FragDepth = z;
    return;
  }

  vec4 rgbM  = texture2Dx(my_texture, texture_coordinate.xy).xyzw;
  vec3 rgbNW = texture2Dx(my_texture, texture_coordinate.xy + (vec2(-1.0, -1.0) * texcoordOffset)).xyz;
  vec3 rgbNE = texture2Dx(my_texture, texture_coordinate.xy + (vec2(+1.0, -1.0) * texcoordOffset)).xyz;
  vec3 rgbSW = texture2Dx(my_texture, texture_coordinate.xy + (vec2(-1.0, +1.0) * texcoordOffset)).xyz;
  vec3 rgbSE = texture2Dx(my_texture, texture_coordinate.xy + (vec2(+1.0, +1.0) * texcoordOffset)).xyz;
	
  vec3 luma = vec3(0.299, 0.587, 0.114);
  float lumaNW = (dot(rgbNW, luma));
  float lumaNE = (dot(rgbNE, luma));
  float lumaSW = (dot(rgbSW, luma));
  float lumaSE = (dot(rgbSE, luma));
  float lumaM  = (dot( rgbM.xyz, luma));

  float lumaMin = min(lumaM, min(min(lumaNW, lumaNE), min(lumaSW, lumaSE)));
  float lumaMax = max(lumaM, max(max(lumaNW, lumaNE), max(lumaSW, lumaSE)));
	
  vec2 dir;
  dir.x = -((lumaNW + lumaNE) - (lumaSW + lumaSE));
  dir.y =  ((lumaNW + lumaSW) - (lumaNE + lumaSE));
	
  float dirReduce = max((lumaNW + lumaNE + lumaSW + lumaSE) * (0.25 * FXAA_REDUCE_MUL), FXAA_REDUCE_MIN);
	  
  float rcpDirMin = 1.0/(min(abs(dir.x), abs(dir.y)) + dirReduce);
	
  dir = min(vec2(FXAA_SPAN_MAX,  FXAA_SPAN_MAX), 
        max(vec2(-FXAA_SPAN_MAX, -FXAA_SPAN_MAX), dir * rcpDirMin)) * texcoordOffset;
		
  vec3 rgbA = (1.0/2.0) * (
              texture2DLod(my_texture, texture_coordinate.xy + dir * (1.0/3.0 - 0.5), rgbM.w).xyz +
              texture2DLod(my_texture, texture_coordinate.xy + dir * (2.0/3.0 - 0.5), rgbM.w).xyz);
  vec3 rgbB = rgbA * (1.0/2.0) + (1.0/4.0) * (
              texture2DLod(my_texture, texture_coordinate.xy + dir * (0.0/3.0 - 0.5), rgbM.w).xyz +
              texture2DLod(my_texture, texture_coordinate.xy + dir * (3.0/3.0 - 0.5), rgbM.w).xyz);
  float lumaB = (dot(rgbB, luma));

  if((lumaB < lumaMin) || (lumaB > lumaMax)){
    gl_FragColor=vec4(rgbA,z);
  } else {
    gl_FragColor=vec4(rgbB,z);
  }
  gl_FragDepth = z;
}
