// "Parthenon" by dr2 - 2015
// License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
#include "setup.inc"
#line 5

float iGlobalTime = time;
vec2 iResolution = vec2(xres, yres);
vec4 iMouse = vec4(0.);

const float pi = 3.14159;
const vec4 cHashA4 = vec4 (0., 1., 57., 58.);
const vec3 cHashA3 = vec3 (1., 57., 113.);
const float cHashM = 43758.54;

vec2 Hashv2f (float p)
{
  return fract (sin (p + cHashA4.xy) * cHashM);
}

vec4 Hashv4f (float p)
{
  return fract (sin (p + cHashA4) * cHashM);
}

vec4 Hashv4v3 (vec3 p)
{
  const vec3 cHashVA3 = vec3 (37.1, 61.7, 12.4);
  const vec3 e = vec3 (1., 0., 0.);
  return fract (sin (vec4 (dot (p + e.yyy, cHashVA3), dot (p + e.xyy, cHashVA3),
     dot (p + e.yxy, cHashVA3), dot (p + e.xxy, cHashVA3))) * cHashM);
}

float Noiseff (float p)
{
  float i, f;
  i = floor (p);  f = fract (p);
  f = f * f * (3. - 2. * f);
  vec2 t = Hashv2f (i);
  return mix (t.x, t.y, f);
}

float Noisefv2 (vec2 p)
{
  vec2 i = floor (p);
  vec2 f = fract (p);
  f = f * f * (3. - 2. * f);
  vec4 t = Hashv4f (dot (i, cHashA3.xy));
  return mix (mix (t.x, t.y, f.x), mix (t.z, t.w, f.x), f.y);
}

float Noisefv3a (vec3 p)
{
  vec3 i, f;
  i = floor (p);  f = fract (p);
  f *= f * (3. - 2. * f);
  vec4 t1 = Hashv4v3 (i);
  vec4 t2 = Hashv4v3 (i + vec3 (0., 0., 1.));
  return mix (mix (mix (t1.x, t1.y, f.x), mix (t1.z, t1.w, f.x), f.y),
              mix (mix (t2.x, t2.y, f.x), mix (t2.z, t2.w, f.x), f.y), f.z);
}

float Fbm3 (vec3 p)
{
  const mat3 mr = mat3 (0., 0.8, 0.6, -0.8, 0.36, -0.48, -0.6, -0.48, 0.64);
  float f, a, am, ap;
  f = 0.;  a = 0.5;
  am = 0.5;  ap = 4.;
  p *= 0.5;
  for (int i = 0; i < 6; i ++) {
    f += a * Noisefv3a (p);
    p *= mr * ap;  a *= am;
  }
  return f;
}

float Fbmn (vec3 p, vec3 n)
{
  vec3 s = vec3 (0.);
  float a = 1.;
  for (int i = 0; i < 5; i ++) {
    s += a * vec3 (Noisefv2 (p.yz), Noisefv2 (p.zx), Noisefv2 (p.xy));
    a *= 0.5;
    p *= 2.;
  }
  return dot (s, abs (n));
}

vec3 VaryNf (vec3 p, vec3 n, float f)
{
  vec3 e = vec3 (0.2, 0., 0.);
  float s = Fbmn (p, n);
  vec3 g = vec3 (Fbmn (p + e.xyy, n) - s,
     Fbmn (p + e.yxy, n) - s, Fbmn (p + e.yyx, n) - s);
  return normalize (n + f * (g - n * dot (n, g)));
}

float PrBoxDf (vec3 p, vec3 b)
{
  vec3 d = abs (p) - b;
  return min (max (d.x, max (d.y, d.z)), 0.) + length (max (d, 0.));
}

float PrSphDf (vec3 p, float s)
{
  return length (p) - s;
}

float PrCapsDf (vec3 p, float r, float h)
{
  return length (p - vec3 (0., 0., h * clamp (p.z / h, -1., 1.))) - r;
}

float PrCylDf (vec3 p, float r, float h)
{
  return max (length (p.xy) - r, abs (p.z) - h);
}

float PrFlatCylDf (vec3 p, float rhi, float rlo, float h)
{
  return max (length (p.xy - vec2 (rhi *
     clamp (p.x / rhi, -1., 1.), 0.)) - rlo, abs (p.z) - h);
}

float SmoothBump (float lo, float hi, float w, float x)
{
  return (1. - smoothstep (hi - w, hi + w, x)) * smoothstep (lo - w, lo + w, x);
}

vec2 Rot2D (vec2 q, float a)
{
  return q * cos (a) * vec2 (1., 1.) + q.yx * sin (a) * vec2 (-1., 1.);
}

int idObj;
vec3 qHit, fCylPos, sunDir;
float tCur, fCylRad, fCylLen, flmFlkr;
const float dstFar = 70.;
bool tryFlm;
const int idLogs = 1, idCoal = 2, idFCyl = 3, idBase = 11, idCol = 12,
   idColEnd = 13, idRoof = 14, idRoofV = 15, idBall = 16, idPool = 17,
   idGal = 18, idAltr = 19, idPost = 20;

vec3 BgCol (vec3 ro, vec3 rd)
{
  vec3 col;
  vec2 p;
  float cloudFac, w, f;
  if (rd.y > 0.) {
    ro.x += 0.1 * tCur;
    p = 0.05 * (rd.xz * (70. - ro.y) / rd.y + ro.xz);
    w = 0.8;  f = 0.;
    for (int j = 0; j < 4; j ++) {
      f += w * Noisefv2 (p);  w *= 0.5;  p *= 2.;
    }
    cloudFac = clamp (3. * f * rd.y - 0.3, 0., 1.);
    f = max (dot (rd, sunDir), 0.);
    col =  mix (vec3 (0.2, 0.3, 0.55) + 0.2 * pow (1. - rd.y, 5.) +
       (0.35 * pow (f, 6.) + 0.65 * min (pow (f, 256.), 0.3)),
       vec3 (0.85), cloudFac);
  } else {
    p = 0.1 * (rd.xz * (10. - ro.y) / rd.y + ro.xz);
    w = 1.;  f = 0.;
    for (int j = 0; j < 4; j ++) {
      f += w * Noisefv2 (p);  w *= 0.7;  p *= 2.5;
    }
    col = mix ((1. + min (f, 1.)) * vec3 (0.15, 0.2, 0.15),
       vec3 (0.2, 0.3, 0.55) + 0.2, pow (1. + rd.y, 5.));
  }
  return col;
}

float BldgDf (vec3 p, float dHit)
{
  vec3 q;
  float d, da, db, wr;
  q = p;
  d = PrBoxDf (q, vec3 (8.6, 0.101, 12.6));
  q.y -= 0.3;
  d = max (min (d, PrBoxDf (q, vec3 (8.2, 0.201, 12.2))),
     - PrBoxDf (q, vec3 (2., 0.25, 6.)));
  if (d < dHit) { dHit = d;  idObj = idBase;  qHit = q; }
  q.y -= 5.52;
  d = max (PrBoxDf (q, vec3 (7.5, 0.05, 11.5)),
     - PrBoxDf (q, vec3 (2.5, 5., 6.5)));
  q.xz = mod (q.xz + vec2 (1.), 2.) - 1.;
  d = max (d, - PrBoxDf (q, vec3 (0.5, 5., 0.5)));
  if (d < dHit) { dHit = d;  idObj = idGal;  qHit = q; }
  q = p;  q.y -= 0.4;
  d = PrBoxDf (q, vec3 (2., 0.01, 6.));
  if (d < dHit) { dHit = d;  idObj = idPool;  qHit = q; }
  q = p;  q.y -= 1.;
  db = max (PrBoxDf (q, vec3 (8., 4.9, 12.)),
     - PrBoxDf (q, vec3 (2., 10., 6.)));
  q = p;  q.xz = mod (q.xz, 2.) - 1.;  q.y -= 3.14;
  wr = q.y / 2.36;
  d = max (PrCylDf (q.xzy, 0.3 * (1.05 - 0.05 * wr * wr), 2.36), db);
  if (d < dHit) { dHit = d;  idObj = idCol;  qHit = q; }
  q = p;  q.xz = mod (q.xz, 2.) - 1.;  q.y = abs (q.y - 3.14) - 2.43;
  d = PrCylDf (q.xzy, 0.4, 0.07);
  q.y -= 0.14;
  d = max (min (d, PrBoxDf (q, vec3 (0.5, 0.07, 0.5))), db);
  if (d < dHit) { dHit = d;  idObj = idColEnd;  qHit = q; }
  q = p;  q.x = abs (q.x) - 3.;  q.y -= 8.2;
  q.xy = Rot2D (q.xy, 0.15 * pi);
  d = PrBoxDf (q, vec3 (6., 0.07, 12.3));
  q.x += 0.4;  q.xz = mod (q.xz, 2.) - 1.;
  d = max (d, - PrBoxDf (q, vec3 (0.5, 5., 0.5)));
  if (d < dHit) { dHit = d;  idObj = idRoof;  qHit = q; }
  q = p;  q.xz = abs (q.xz);  q -= vec3 (4.1, 7.68, 11.6);
  d = PrBoxDf (q, vec3 (4.3, 1.9, 0.1));
  q.xy = Rot2D (q.xy, 0.15 * pi);
  q.xy -= vec2 (-0.4, -2.);
  d = max (d, PrBoxDf (q, vec3 (4.3, 1.9, 0.1)));
  q = p;  q.xz = abs (q.xz);  q -= vec3 (2.89, 7.7, 11.6);
  da = PrBoxDf (q, vec3 (3.2, 1.4, 1.));
  q.xy = Rot2D (q.xy, 0.15 * pi);
  q.xy -= vec2 (-0.25, -1.5);
  d = max (d, - max (da, PrBoxDf (q, vec3 (3.2, 1.4, 1.))));
  if (d < dHit) { dHit = d;  idObj = idRoofV;  qHit = q; }
  q = p;  q.y -= 7.7;  q.z = abs (q.z) - 11.6;
  d = PrCylDf (q.xzy, 0.09, 1.4);
  if (d < dHit) { dHit = d;  idObj = idRoofV;  qHit = q; }
  d = PrSphDf (q, 0.4);
  if (d < dHit) { dHit = d;  idObj = idBall;  qHit = q; }
  q = p;  q.xz = abs (q.xz);  q -= vec3 (8.5, 0.6, 12.5);
  d = PrCylDf (q.xzy, 0.05, 0.5);
  if (d < dHit) { dHit = d;  idObj = idPost;  qHit = q; }
  q.y -= 0.7;
  d = PrSphDf (q, 0.2);
  if (d < dHit) { dHit = d;  idObj = idBall;  qHit = q; }
  q = p;  q.y -= 1.5;  
  d = max (max (PrSphDf (q, 0.78), -0.01 + q.y), -0.3 - q.y);
  q.y -= -0.7;
  d = min (d, PrCylDf (q.xzy, 0.15, 0.42));
  if (d < dHit) { dHit = d;  idObj = idAltr;  qHit = q; }
  return dHit;
}

float FireDf (vec3 p, float dHit)
{
  vec3 q;
  float d;
  q = p;  q.x = abs (q.x) - 0.3;
  q.y -= fCylPos.y - fCylLen + 0.09;
  d = PrCapsDf (q, 0.12 - 0.04 * Noisefv3a (15. * p), 0.5);
  if (d < dHit) { dHit = d;  qHit = p;  idObj = idLogs; }
  q = p;  q.y -= fCylPos.y - fCylLen + 0.25;
  q.z = abs (q.z) - 0.25;
  d = PrCapsDf (q.zyx, 0.12 - 0.04 * Noisefv3a (15. * p), 0.45);
  if (d < dHit) { dHit = d;  qHit = p;  idObj = idLogs; }
  q = p;  q.x = abs (q.x) - 0.2;
  q.y -= fCylPos.y - fCylLen + 0.43;
  d = PrCapsDf (q, 0.12 - 0.04 * Noisefv3a (15. * p), 0.4);
  if (d < dHit) { dHit = d;  qHit = p;  idObj = idLogs; }
  q = p;  q.y -= fCylPos.y - fCylLen - 0.01;
  d = PrCylDf (q.xzy, fCylRad, 0.01);
  if (d < dHit) { dHit = d;  qHit = q;  idObj = idCoal; }
  return dHit;
}

float ObjDf (vec3 p)
{
  vec3 q;
  float dHit, d;
  dHit = dstFar;
  if (tryFlm) {
    q = p;  q -= fCylPos;
    d = PrCylDf (q.xzy, fCylRad, fCylLen);
    if (d < dHit) { dHit = d;  idObj = idFCyl;  qHit = q; }
    q = p;  q.y -= 0.4;
    d = PrBoxDf (q, vec3 (2., 0.01, 6.));
    if (d < dHit) { dHit = d;  idObj = idPool;  qHit = q; }
  } else {
    dHit = 0.9 * BldgDf (p, dHit);
    dHit = FireDf (p, dHit);
  }
  return dHit;
}

float ObjRay (vec3 ro, vec3 rd)
{
  float d;
  float dHit = 0.;
  for (int j = 0; j < 200; j ++) {
    d = ObjDf (ro + dHit * rd);
    dHit += d;
    if (d < 0.001 || dHit > dstFar) break;
  }
  return dHit;
}

vec3 ObjNf (vec3 p)
{
  const vec3 e = vec3 (0.001, -0.001, 0.);
  vec4 v = vec4 (ObjDf (p + e.xxx), ObjDf (p + e.xyy),
     ObjDf (p + e.yxy), ObjDf (p + e.yyx));
  return normalize (vec3 (v.x - v.y - v.z - v.w) + 2. * v.yzw);
}

float ObjSShadow (vec3 ro, vec3 rd, float dLight)
{
  float sh = 1.;
  float d = 0.15;
  for (int i = 0; i < 30; i++) {
    float h = ObjDf (ro + rd * d);
    sh = min (sh, 20. * h / d);
    d += max (0.15, 0.01 * d);
    if (h < 0.01 || d > dLight) break;
  }
  return clamp (sh, 0., 1.);
}

vec3 FireSrcCol (vec3 n)
{ // Inspired by Dave_H's "Campfire", but implementation differs.
  vec3 col, q;
  float f, di, gl, bri;
  if (idObj == idLogs) {
    n = VaryNf (0.4 * qHit, n, 5.);
    q = 45. * qHit;  q.y -= 0.8 * tCur;
    f = Noisefv3a (q);   
    f += Noisefv3a (10. * qHit) + 0.4 * Noisefv3a (123. * qHit);
    f = 0.01 * pow (abs (f), 13.) + 0.7 * max (1. - 0.11 * dot (qHit, qHit), 0.);
    bri = max (dot (normalize (vec3 (0., 0.65 + flmFlkr, 0.) - qHit), n), 0.4);
    col = bri * bri * f * vec3 (0.5, 0.1, 0.);
    di = 0.;
    gl = 0.;
    for (int i = 0; i < 3; i ++) {
      di += 0.02;
      gl += max (0., di - 1.2 * ObjDf (qHit + di * n));
    }
    col += flmFlkr * flmFlkr * clamp (2. * gl, 0., 1.) * vec3 (1., 0.3, 0.05);
  } else if (idObj == idCoal) {
    q.xz = 20. * qHit.xz;  q.y = tCur;
    f = Noisefv3a (11. * q);
    q.y = 0.04 * tCur;
    f += Noisefv3a (3. * q);
    bri = 3. - 0.5 * flmFlkr;
    col = 0.1 * bri * bri * pow (f, 4.) * vec3 (0.7, 0.05, 0.) *
       (1.1 - pow (length (qHit.xz) / fCylRad, 4.));
  }
  return col;
}

float FireLum (vec3 ro, vec3 rd, float dHit)
{
  vec3 p, q, dp;
  float g, s, fh, fr, f;
  p = ro;
  dp = (fCylRad / 40.) * rd;
  g = 0.;
  for (int i = 0; i < 40; i ++) {
    p += dp;
    s = distance (p.xz, fCylPos.xz);
    q = 4. * p;  q.y -= 6. * tCur;
    fh = 0.5 * max (1. - (p.y - fCylPos.y) / fCylLen, 0.);
    fr = max (1. - s / fCylRad, 0.);
    f = Fbm3 (q);
    q = 7. * p;  q.y -= 8.5 * tCur;
    f += Fbm3 (q);
    g += max (0.5 * fr * fr * fh * (f * f - 0.6), 0.);
    q = 23. * p;  q.y -= 11. * tCur;
    g += 1000. * pow (abs (Noisefv3a (q) - 0.11), 64.);
    if (s > fCylRad || p.y < fCylPos.y - 0.99 * fCylLen || g > 1.) break;
  }
  return g;
}

vec3 ObjCol (vec3 n)
{
  vec3 col;
  if (idObj == idBase) col = vec3 (0.8, 0.8, 0.7);
  else if (idObj == idCol || idObj == idColEnd) col = vec3 (0.7, 0.8, 0.6);
  else if (idObj == idRoof) col = vec3 (0.8, 0.1, 0.1) *
      (1. - 0.6 * SmoothBump (0.98, 1.02, 0.01, mod (qHit.z, 2.))) *
      (1. - 0.6 * SmoothBump (0.98, 1.02, 0.01, mod (qHit.x, 2.)));
  else if (idObj == idRoofV) col = vec3 (0.8, 0.1, 0.1);
  else if (idObj == idPost) col = vec3 (0.8, 0.1, 0.1);
  else if (idObj == idBall) col = vec3 (1., 1., 0.1);
  else if (idObj == idGal) col = vec3 (0.1, 0.3, 0.1);
  else if (idObj == idAltr) col = vec3 (0.6, 0.5, 0.2);
  return col;
}

vec4 ShowScene (vec3 ro, vec3 rd)
{
  vec3 roo, rdo, objCol, flmCol, foVec, col, vn;
  float dstHit, dstFlm, dstHitR, dstFlmR, fIntens, fLum, lDist, sh, bk,
     dif, ltExt, reflFac, a, f;
  int idObjT;
  bool flmRefl, objRefl;
  fCylPos = vec3 (0., 3.53, 0.);
  fCylRad = 0.8;
  fCylLen = 2.;
  roo = ro;
  rdo = rd;
  idObj = -1;
  tryFlm = true;
  dstFlm = ObjRay (roo, rdo);
  flmRefl = false;
  if (dstFlm < dstFar && idObj == idPool) {
    dstFlmR = dstFlm;
    roo += dstFlm * rdo;
    vn = ObjNf (roo);
vn += .1 * cross(roo, .1*vec3(sin(time), cos(time), 0.));
vn = normalize(vn);
    rdo = reflect (rdo, vn);
    roo += 0.01 * rdo;
    idObj = -1;
    dstFlm = ObjRay (roo, rdo);
    flmRefl = true;
  }
  fIntens = (dstFlm < dstFar) ? FireLum (roo + dstFlm * rdo, rdo, dstFlm) : 0.;
  flmFlkr = Noiseff (tCur * 64.);
  idObj = -1;
  tryFlm = false;
  dstHit = ObjRay (ro, rd);
  if (idObj < 0) dstHit = dstFar;
  objRefl = false;
  reflFac = 1.;
  if (dstHit < dstFar && idObj == idPool) {
    dstHitR = dstHit;
    ro += dstHit * rd;
    vn = ObjNf (ro);
vn += .1 * cross(ro, .1*vec3(sin(time), cos(time), 0.));
vn = normalize(vn);
    rd = reflect (rd, VaryNf (0.4 * qHit, vn, 0.1));
    ro += 0.01 * rd;
    idObj = -1;
    dstHit = ObjRay (ro, rd);
    if (idObj < 0) dstHit = dstFar;
    reflFac = 0.8;
    objRefl = true;
  }
  idObjT = idObj;
  if (dstHit >= dstFar) col = 0.1 * BgCol (ro, rd);
  else {
    ro += rd * dstHit;
    vn = ObjNf (ro);
    if (idObj == idBase) vn = VaryNf (10. * qHit, vn, 0.4);
    if (idObj == idCol) {
      a = 0.5 - mod (12. * (atan (qHit.x, qHit.z) / (2. * pi) + 0.5), 1.);
      vn.xz = Rot2D (vn.xz, -0.1 * pi * sin (pi * a));
    }
    if (idObj == idCol || idObj == idColEnd) vn = VaryNf (20. * qHit, vn, 0.3);
    else if (idObj == idAltr) vn = VaryNf (10. * qHit, vn, 1.);
    idObj = idObjT;
    if (idObj == idLogs || idObj == idCoal) col = FireSrcCol (vn);
    else {
      objCol = ObjCol (vn);
      foVec = fCylPos - ro;
      lDist = length (foVec);
      foVec /= lDist;
      sh = ObjSShadow (ro, foVec, lDist);
      fLum = 5. * sh * (0.6 + 0.4 * flmFlkr) / pow (lDist, 1.5) *
         max (dot (foVec, vn), 0.);
      dif = max (dot (vn, sunDir), 0.);
      bk = max (dot (vn, - normalize (vec3 (sunDir.x, 0., sunDir.z))), 0.);
      ltExt = 0.1 * (0.4 * (1. + bk) + max (0., dif));
      col = objCol * (ltExt * vec3 (1.) + fLum * vec3 (1., 0.3, 0.2));
    }
  }
  idObj = idObjT;
  if (dstHit < dstFar && idObj == idBall) col *= 3.;
  if (! (dstHit < dstFar && (idObj == idCoal || idObj == idAltr) ||
     flmRefl == objRefl && dstHit < dstFlm ||
     objRefl && ! flmRefl && dstHitR < dstFlm ||
     ! objRefl && flmRefl && dstHit < dstFlmR)) {
    f = clamp (0.7 * fIntens, 0., 1.);
    f *= f;
    f *= f;
    flmCol = 1.5 * mix (vec3 (1., 0.2, 0.2), vec3 (1., 1., 0.5), f);
    fIntens *= fIntens;
    if (idObj == idLogs) fIntens *= 0.03;
    col = mix (col, flmCol, min (fIntens, 1.));
  }
  col *= reflFac;
  col += (1. - reflFac) * vec3(0.,.2,.5);
  return vec4(col, dstHit);
}

vec3 SetVuPt (float t)
{
  vec4 va, vb;
  vec3 wPt[11], ro;
  float tp[11], vel, tVu;
  wPt[0] = vec3 (0., 8., -20.);     wPt[1] = vec3 (0., 2., -12.);
  wPt[2] = vec3 (-0.5, 2.5, -6.2);  wPt[3] = vec3 (-2.6, 3., -6.2);
  wPt[4] = vec3 (-2.6, 4., 8.1);    wPt[5] = vec3 (0., 5.6, 8.9);
  wPt[6] = vec3 (2.6, 4., 8.1);     wPt[7] = vec3 (2.6, 3., -6.2);
  wPt[8] = wPt[2];  wPt[9] = wPt[1];  wPt[10] = wPt[0];
  tp[0]  = 0.;
  tp[1]  = tp[0] + distance (wPt[1],  wPt[0]); 
  tp[2]  = tp[1] + distance (wPt[2],  wPt[1]); 
  tp[3]  = tp[2] + distance (wPt[3],  wPt[2]); 
  tp[4]  = tp[3] + distance (wPt[4],  wPt[3]); 
  tp[5]  = tp[4] + distance (wPt[5],  wPt[4]); 
  tp[6]  = tp[5] + distance (wPt[6],  wPt[5]); 
  tp[7]  = tp[6] + distance (wPt[7],  wPt[6]); 
  tp[8]  = tp[7] + distance (wPt[8],  wPt[7]); 
  tp[9]  = tp[8] + distance (wPt[9],  wPt[8]); 
  tp[10] = tp[9] + distance (wPt[10], wPt[9]);
  vel = 0.5;
  tVu = mod (t, tp[10] / vel) * vel;
  if (tVu < tp[5]) {
    if (tVu < tp[1])      { va = vec4 (wPt[0], tp[0]);  vb = vec4 (wPt[1], tp[1]); }
    else if (tVu < tp[2]) { va = vec4 (wPt[1], tp[1]);  vb = vec4 (wPt[2], tp[2]); }
    else if (tVu < tp[3]) { va = vec4 (wPt[2], tp[2]);  vb = vec4 (wPt[3], tp[3]); }
    else if (tVu < tp[4]) { va = vec4 (wPt[3], tp[3]);  vb = vec4 (wPt[4], tp[4]); }
    else                  { va = vec4 (wPt[4], tp[4]);  vb = vec4 (wPt[5], tp[5]); }
  } else {
    if (tVu < tp[6])      { va = vec4 (wPt[5], tp[5]);  vb = vec4 (wPt[6], tp[6]); }
    else if (tVu < tp[7]) { va = vec4 (wPt[6], tp[6]);  vb = vec4 (wPt[7], tp[7]); }
    else if (tVu < tp[8]) { va = vec4 (wPt[7], tp[7]);  vb = vec4 (wPt[8], tp[8]); }
    else if (tVu < tp[9]) { va = vec4 (wPt[8], tp[8]);  vb = vec4 (wPt[9], tp[9]); }
    else                  { va = vec4 (wPt[9], tp[9]);  vb = vec4 (wPt[10], tp[10]); }
  }
  ro = mix (va.xyz, vb.xyz, (tVu - va.w) / (vb.w - va.w));
  return ro;
}

void main (void)
{
  vec2 uv = 2. * gl_FragCoord.xy / iResolution.xy - 1.;
  uv.x *= iResolution.x / iResolution.y;
  tCur = iGlobalTime;
  vec4 mPtr = iMouse;
  mPtr.xy = mPtr.xy / iResolution.xy - 0.5;
  mat3 vuMat;
  vec3 ro, rd, vd, u;
  float zmFac, hLkAt, f;
  sunDir = normalize (vec3 (1., 1., 1.));
  hLkAt = 3.;
  zmFac = 1.4;
  if (mPtr.z <= 0.) {
    ro = 0.5 * (SetVuPt (tCur + 0.2) + SetVuPt (tCur - 0.2));
  } else {
    ro = SetVuPt (mod (150. * (mPtr.x + 0.45), 1000.));
    hLkAt += 5. * mPtr.y;
  }
  vd = normalize (vec3 (0., hLkAt, 0.) - ro);
  u = - vd.y * vd;
  f = 1. / sqrt (1. - vd.y * vd.y);
  vuMat = mat3 (f * vec3 (vd.z, 0., - vd.x), f * vec3 (u.x, 1. + u.y, u.z), vd);
  rd = vuMat * normalize (vec3 (uv, zmFac));

  if (!setup_ray(eye, dir, ro, rd)) return;

  vec4 col = ShowScene (ro, rd);
  col.rgb = sqrt (clamp (col.rgb, 0., 1.));

//  gl_FragColor = vec4 (col, 1.);
  write_pixel(dir, col.a, col.rgb);
}
