//based on shader from coyote => https://www.shadertoy.com/view/ltfGzS

#include "setup.inc"
#line 6

vec2 iResolution = vec2(xres, yres);
float iGlobalTime = time;
uniform int iters;
uniform vec3 par[10];
#define Scale1 par[0].x  // {min=0 max=10 step=.01}
#define Scale2 par[0].y  // {min=0 max=10 step=.01}
#define Scale3 par[0].z  // {min=0 max=10 step=.01}
#define powerVector par[1]

void main() {
    vec3 p, q, r;
    if (!setup_ray(eye, dir, q, p)) return;
    
    float total = 0.;
    vec3 col = vec3(0.);

    for (float i=1.; i>0.; i-=.01) {
        float d=0.,s=1.;

        for (int j = 0; j < iters; j++) {
            r=max(r=pow(abs(mod(q*s+1.,2.)-1.), powerVector),r.yzx);
            d=max(d,(Scale1-length(r)*Scale3)/s);
            s*=Scale2;
        }

        q+=p*d;
        total+=d;

        col = p-p+i;

        if (d<max(1e-5,1e-5*total)) break;
    }

    write_pixel(dir, total, col);
}

