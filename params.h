#ifndef _F_PARAMS_H_
#define _F_PARAMS_H_

// Set of properties every shader / keyframe likely has.
// Apply some macro FU for meta programming.
// [type, variable name, config name, spline/uniform?]

#define PROCESS_COMMON_PARAMS \
  PROCESS(int, width, "width", false) \
  PROCESS(int, height, "height", false) \
  PROCESS(int, fullscreen, "fullscreen", false) \
  PROCESS(int, multisamples, "multisamples", false) \
  PROCESS(float, keyb_rot_speed, "keyb_rot_speed", false) \
  PROCESS(float, mouse_rot_speed, "mouse_rot_speed", false) \
  PROCESS(float, fov_x, "fov_x", true) \
  PROCESS(float, fov_y, "fov_y", true) \
  PROCESS(double, speed, "speed", true) \
  PROCESS(float, min_dist, "min_dist", true) \
  PROCESS(int, max_steps, "max_steps", true) \
  PROCESS(int, iters, "iters", true) \
  PROCESS(int, color_iters, "color_iters", true) \
  PROCESS(int, loop, "loop", false) \
  PROCESS(float, ao_eps, "ao_eps", true) \
  PROCESS(float, ao_strength, "ao_strength", true) \
  PROCESS(float, glow_strength, "glow_strength", true) \
  PROCESS(float, dist_to_color, "dist_to_color", true) \
  PROCESS(double, delta_time, "delta_time", false) \
  PROCESS(double, time, "time", true) \
  PROCESS(float, fps, "fps", false) \
  PROCESS(int, depth_size, "depth_size", false) \
  PROCESS(float, aperture, "aperture", true) \
  PROCESS(float, dof_offset, "dof_offset", true) \
  PROCESS(int, enable_dof, "enable_dof", false) \
  PROCESS(int, enable_fxaa, "enable_fxaa", false) \
  PROCESS(int, disable_de, "disable_de", false) \
  PROCESS(int, no_spline, "no_spline", false) \
  PROCESS(float, focus, "focus", true) \
  PROCESS(int, nrays, "nrays", true) \
  PROCESS(int, xres, "xres", false) \
  PROCESS(int, yres, "yres", false) \
  PROCESS(int, use_bg_texture, "use_bg_texture", false) \
  PROCESS(int, backbuffer, "backbuffer", false) \
  PROCESS(int, julia, "julia", true) \
  PROCESS(int, iBackbufferCount, "iBackbufferCount", false) \
  PROCESS(int, frameno, "frameno", false) \
  PROCESS(int, fxaa, "fxaa", false) \
  PROCESS(float, ipd, "ipd", false) \
  PROCESS(float, exposure, "exposure", true) \
  PROCESS(float, maxBright, "maxBright", true) \
  PROCESS(float, gamma, "gamma", true)

#define NUMPARS 20

#endif // _F_PARAMS_H_
