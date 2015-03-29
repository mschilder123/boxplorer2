// by @301z; based on http://meatfighter.com/juggler/
//
// This has been kicked out of Shadertoy due to lack of stability.
// Handle with care!

#include "setup.inc"
#line 8

uniform int iters;    // {min=1 max=1000} Number of fractal iterations.
uniform int color_iters;    // {min=1 max=1000} Number of fractal iterations.

//precision mediump float;

//uniform float time;
//uniform vec2 resolution;

struct vec3x3 {
	vec3 x;
	vec3 y;
	vec3 z;
};

vec3 pose(float time) {
	float angle = time * (2.0 * 3.14159);
	float oscillation = cos(angle) * 0.5 + 0.5;
	float armAngle = oscillation * -0.35;
	return vec3(angle, oscillation, armAngle);
}

vec3x3 body(vec3 pose) {
	vec3x3 body;
	body.x = pose * vec3(0.0, 4.0, 0.0) + vec3(151.0, 81.0, -151.0);
	body.y = normalize(sin(pose).xxx * vec3(0.0, 0.0, -4.0) + vec3(0.0, 70.0, 0.0));
	body.z = body.y.xzy * vec3(0.0, 1.0, -1.0);
	return body;
}

vec3x3 appendage(vec3 p, vec3 q, vec3 w, float a, float b) {
	vec3 pq = q - p;
	float dd = dot(pq, pq);
	float rd = inversesqrt(dd);
	float aa = a * a;
	float y = (aa - b * b + dd) * rd * 0.5;
	vec3 v = pq * rd;
	return vec3x3(p, q, p + sqrt(abs(aa - y * y)) * cross(v, w) + y * v);
}

vec3 interpolate(vec3x3 appendage, int i) {
	float a = float((i < 8) ? i : (i - 8)) * (1.0 / 8.0);
	return mix((i < 8) ? appendage.x : appendage.y, appendage.z, a);
}

vec3 head(vec3x3 body) {return vec3(body.x + body.y * 70.0);}
vec3 glabella(vec3x3 body) {return body.x + vec3(-9.0, 0.0, 0.0) + body.y * vec3(0.0, 69.0, 69.0);}
vec3 eyeL(vec3x3 body) {return vec3(glabella(body) + body.z * vec3(0.0, -7.0, -7.0));}
vec3 eyeR(vec3x3 body) {return vec3(glabella(body) + body.z * vec3(0.0, +7.0, +7.0));}
vec3 neck(vec3x3 body) {return vec3(body.x + body.y * 55.0);}
vec3 hair(vec3x3 body) {return vec3(body.x + vec3(1.0, 0.0, 0.0) + body.y * vec3(0.0, 71.0, 71.0));}
vec3 torso(vec3x3 body, int i) {return vec3(body.x + body.y * float(i) * (32.0 / 7.0));}
vec3x3 leg(vec3 p, vec3 q, vec3 w) {return appendage(p, q, w, 42.580, 34.070);}
vec3x3 legL(vec3x3 body) {return leg(vec3(159.0, 2.5, -133.0), body.x + body.y * -9.0 + body.z * -16.0, body.z);}
vec3x3 legR(vec3x3 body) {return leg(vec3(139.0, 2.5, -164.0), body.x + body.y * -9.0 + body.z * +16.0, body.z);}
vec3x3 arm(vec3 p, vec3 q, vec3 w) {return appendage(p, q, w, 44.294, 46.098);}
vec2 arm(vec3 pose) {return vec2(vec2(cos(pose.z), sin(pose.z)) * vec2(41.0, -41.0) + vec2(69.0, 60.0));}
vec3x3 armL(vec3x3 body, vec3 p) {return arm(p, body.x + body.y * 45.0 + body.z * -19.0, normalize(body.y * +0.4 + body.z * -38.9));}
vec3x3 armR(vec3x3 body, vec3 p) {return arm(p, body.x + body.y * 45.0 + body.z * +19.0, normalize(body.y * -0.4 + body.z * -0.9));}
vec3x3 armL(vec3 pose) {return armL(body(pose), vec3(arm(pose), -108.0));}
vec3x3 armR(vec3 pose) {return armR(body(pose), vec3(arm(pose), -182.0));}

vec3 ball0(float time) {return vec3(110.0, time * 96.0 + time * time * -96.0 + 88.0, time * 74.0 - 182.0);}
vec3 ball1(float time) {return vec3(110.0, time * time * -96.0 + 184.0, time * -37.0 - 145.0);}
vec3 ball2(float time) {return vec3(110.0, time * 192.0 + time * time * -96.0 + 88.0, time * -37.0 - 108.0);}
vec3 head(float time) {return head(body(pose(time)));}
vec3 eyeL(float time) {return eyeL(body(pose(time)));}
vec3 eyeR(float time) {return eyeR(body(pose(time)));}
vec3 neck(float time) {return neck(body(pose(time)));}
vec3 hair(float time) {return hair(body(pose(time)));}
vec3 torso(float time, int i) {return torso(body(pose(time)), i);}
vec3 legL(float time, int i) {return interpolate(legL(body(pose(time))), i);}
vec3 legR(float time, int i) {return interpolate(legR(body(pose(time))), i);}
vec3 armL(float time, int i) {return interpolate(armL(pose(time)), i);}
vec3 armR(float time, int i) {return interpolate(armR(pose(time)), i);}

struct Material {
	vec4 ambientDiffuseSpecularReflection;
	vec4 colorShininess;
};

const float kSkyMaterial = 0.0;
const float kGroundMaterial = 1.0;
const float kBallMaterial = 2.0;
const float kSkinMaterial = 3.0;
const float kEyeMaterial = 4.0;
const float kHairMaterial = 5.0;
const float kTorsoMaterial = 6.0;

Material xlat(vec3 position, vec4 incidenceMaterial) {
	const vec4 kMatte = vec4(0.1, 1.5, 0.0, 0.0);
	const vec4 kSky = vec4(1.0, 0.0, 0.0, 0.0);
	const vec4 kMetal = vec4(0.0, 0.0, 1.0, 1.0);
	const vec4 kPlastic = vec4(0.0, 1.0, 1.0, 0.0);
	vec2 stripes = step(vec2(0.5), fract(position.xz * (1.0 / 214.0)));
	Material checkerBoard = Material(kMatte, vec4(mix(stripes.x, 1.0 - stripes.x, stripes.y), 1.0, 0.0, 10.0));
	Material material = Material(kSky, mix(vec4(0.5, 0.5, 1.0, 10.0), vec4(0.0, 0.0, 0.9, 10.0), incidenceMaterial.y));
	material = (incidenceMaterial.w == kGroundMaterial) ? checkerBoard : material;
	material = (incidenceMaterial.w == kBallMaterial) ? Material(kMetal, vec4(1.0, 1.0, 1.0, 20.0)) : material;
	material = (incidenceMaterial.w == kSkinMaterial) ? Material(kPlastic, vec4(1.0, 0.5, 0.5, 10.0)) : material;
	material = (incidenceMaterial.w == kEyeMaterial) ? Material(kPlastic, vec4(0.0, 0.0, 0.6, 10.0)) : material;
	material = (incidenceMaterial.w == kHairMaterial) ? Material(kPlastic, vec4(0.0, 0.0, 0.0, 10.0)) : material;
	material = (incidenceMaterial.w == kTorsoMaterial) ? Material(kPlastic, vec4(0.9, 0.0, 0.0, 10.0)) : material;
	return material;
}

struct Hit {
	float distance;
	vec4 centerMaterial;
};

Hit trace(Hit hit, vec3 rayOrigin, vec3 rayDirection, float material) {
	float distance = -rayOrigin.y / rayDirection.y;
	vec3 center = rayOrigin + rayDirection * distance + vec3(0.0, -1.0, 0.0);
	hit = all(bvec2(distance > 0.0, distance < hit.distance)) ? Hit(distance, vec4(center, material)) : hit;
	return hit;
}

Hit trace(Hit hit, vec3 rayOrigin, vec3 rayDirection, vec3 center, float radius, float material) {
	vec3 ro = center - rayOrigin;
	float b = dot(ro, rayDirection);
	float t = b * b - dot(ro, ro) + radius * radius;
	float distance = b - sqrt(t);
	hit = all(bvec3(t > 0.0, distance > 0.0, distance < hit.distance)) ? Hit(distance, vec4(center, material)) : hit;
	return hit;
}

Hit trace(Hit hit, vec3 rayOrigin, vec3 rayDirection, vec3 center0, vec3 center1, vec2 radii, vec2 materials) {
	vec3 ro0 = center0 - rayOrigin;
	vec3 ro1 = center1 - rayOrigin;
	vec2 b = vec2(dot(ro0, rayDirection), dot(ro1, rayDirection));
	vec2 t = b * b - vec2(dot(ro0, ro0), dot(ro1, ro1)) + radii * radii;
	vec2 distances = b - sqrt(t);
	bvec2 tz = greaterThan(t, vec2(0.0));
	bvec2 dz = greaterThan(distances, vec2(0.0));
	hit = all(bvec3(tz.x, dz.x, distances.x < hit.distance)) ? Hit(distances.x, vec4(center0, materials.x)) : hit;
	hit = all(bvec3(tz.y, dz.y, distances.y < hit.distance)) ? Hit(distances.y, vec4(center1, materials.y)) : hit;
	return hit;
}

Hit trace(vec3 rayOrigin, vec3 rayDirection, float time) {
	const float kFar = 1e6;
	Hit hit = Hit(kFar, vec4(rayOrigin + rayDirection * (kFar + 1.0), kSkyMaterial));
	hit = trace(hit, rayOrigin, rayDirection, kGroundMaterial);
	hit = trace(hit, rayOrigin, rayDirection, ball0(time), ball1(time), vec2(14.0), vec2(kBallMaterial));
	hit = trace(hit, rayOrigin, rayDirection, ball2(time), head(time), vec2(14.0), vec2(kBallMaterial, kSkinMaterial));
	hit = trace(hit, rayOrigin, rayDirection, eyeL(time), eyeR(time), vec2(4.0), vec2(kEyeMaterial));
//	hit = trace(hit, rayOrigin, rayDirection, neck(time), hair(time), vec2(5.0, 14.0), vec2(kSkinMaterial, kHairMaterial));
	for (int i = 0; i < 17; i++) {
		if (i < 8)
			hit = trace(hit, rayOrigin, rayDirection, torso(time, i), float(i) * (4.0 / 7.0) + 16.0, kTorsoMaterial);
		float radius = min(float(i) * (2.5 / 7.0) + 2.5, 5.0);
		hit = trace(hit, rayOrigin, rayDirection, legL(time, i), radius, kSkinMaterial);
		hit = trace(hit, rayOrigin, rayDirection, legR(time, i), radius, kSkinMaterial);
		hit = trace(hit, rayOrigin, rayDirection, armL(time, i), radius, kSkinMaterial);
		hit = trace(hit, rayOrigin, rayDirection, armR(time, i), radius, kSkinMaterial);
	}
	return hit;
}

vec3 lookAt(vec3 eye, vec3 interest, vec2 uv) {
	vec3 forward = normalize(interest - eye);
	vec3 left = forward.zyx * vec3(-1.0, 0.0, 1.0);
	return normalize(forward + left * uv.x + cross(left, forward) * uv.y);
}

void main() {
	vec4 colorReflectivity = vec4(0.0, 0.0, 0.0, 1.0);
#if 0
	vec3 rayOrigin = vec3(2.0, 100.0, -2.0);
	vec3 rayDirection = lookAt(rayOrigin, vec3(1000.0, 77.0, -1000.0), (gl_FragCoord.xy * 2.0 - resolution) / resolution.xx);
#else
  vec3 rayOrigin, rayDirection;
  if (!setup_ray(eye, dir, rayOrigin, rayDirection)) {
    return;
  }
#endif
	vec3 normal;
	vec4 incidenceMaterial;
	for (int i = 0; i < 6; i++) {
		Hit hit = trace(rayOrigin, rayDirection, fract(time));
		if ((i - (i / 2) * 2) < 1) {
			vec3 position = rayOrigin + rayDirection * hit.distance;
			normal = normalize(position - hit.centerMaterial.xyz);
			incidenceMaterial = vec4(rayDirection, hit.centerMaterial.w);
			rayDirection = normalize(vec3(-564.0, 686.0, 147.0) - position);
			rayOrigin = position + rayDirection * 0.001;
		} else {
			float nDotL = dot(normal, rayDirection);
			float nDotH = dot(normal * nDotL * 2.0 - rayDirection, -incidenceMaterial.xyz);
			Material material = xlat(rayOrigin, incidenceMaterial);
			vec3 ambient = material.colorShininess.xyz * material.ambientDiffuseSpecularReflection.x;
			vec3 diffuse = material.colorShininess.xyz * material.ambientDiffuseSpecularReflection.y * max(0.0, nDotL);
			float specular = pow(max(0.0, nDotH), material.colorShininess.w) * material.ambientDiffuseSpecularReflection.z;
	 		colorReflectivity.xyz += (ambient + (diffuse + specular) * step(1000.0, hit.distance)) * colorReflectivity.w;
			colorReflectivity.w *= material.ambientDiffuseSpecularReflection.w;
	 		rayDirection = reflect(incidenceMaterial.xyz, normal);
		}
	}
	//gl_FragColor = vec4(pow(colorReflectivity.xyz, vec3(1.0 / 2.2)), 1.0);
	write_pixel(dir, 0.0, pow(colorReflectivity.xyz, vec3(1.0 / 2.2)));
}
