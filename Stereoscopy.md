# Stereoscopy #

Boxplorer2 understands the following commandline parameters to tell it how to render 3d output:
  * _--overunder_ : half height, left eye on top. Compatible with 3d tvs.
  * _--sidebyside_ : half width, right eye on right. Compatible with 3d tvs.
  * _--xeyed_ : full width, right eye on left. Good for free viewing.
  * _--interlaced_ : horizontal interlaced. Compatible with Zalman and LG D2342P-PN like devices.
  * _--quadbuffer_ : quad buffer. Press the p key to switch L/R polarity. Tested on Radeon 7970 + HD600X-LV and nVidia Titan.
  * _--oculus_: Oculus Rift demo kit compatible output.
  * _--anaglyph_: Old-school red&blue glasses might work.

# Details #

To make the 3d easy on the eyes, boxplorer2 has a couple of parameters you can tweak.
  * _speed_ : half the distance between the virtual eyes (aka virtual head size).
  * _focus_ : distance of focal plane from eyes. Default value of 0 places the focus at 20x _speed_. Negative focus puts it closer. The focal plane is where the rays cast from each of the left and right eye converge, for their respective pixels. Hence, objects at the focal plane distance appear right at the screen plane.

When navigating a compatible configuration, boxplorer2 computes the distance estimation for the current position and sets _speed_ to that estimate / 10. As the name _speed_ suggests, this value is also used for speed of movement, so getting closer to the fractal surface will both slow the rate of movement down as well as shrink your virtual head.

You can control _speed_ and _focus_ manually by selecting the view controller with the v key and then control the two values with arrows keys.

For _--oculus_ output, there is not asymmetric view fustrum and _focus_ gets used to tweak your inter-pupillary distance instead.

# Try it! #

Not all configuration files have been updated to do 3d the proper way yet, but for a nice selection of 3d views of the mandelbox try:
```
boxplorer2 cfgs/rrrola/ --xeyed
```
Use the tab key to flip forward through the various scenes and the backspace key to flip backwards. Move around in a scene using your mouse and w,a,s,d keys. Adjust the _--xeyed_ to match your preferred 3d output method.

The enter key toggles between windowed and fullscreen view. For fullscreen _--interlaced_ the resolution will be set to 1920x1080. For _--overunder_ or _--sidebyside_ the resolution will be set to 1280x720 (native resolution of Sony's HMZ-T1).

Other configuration files that should work: colonoscopy, mbox-fp64, menger, menger-sphere, knighty/HypTess, knighty/Reflectoids, knighty/SpentFuelRods, syntopia/spherical and most of shadertoys.

# Save it! #

Found an awesome 3d scene you want to record for posterity? Hit the printscreen key and a screen shot and configuration file will be saved in the _cfgs/rrrola/_ directory.
Hit the space key and an additional keyframe will be saved there, and show up as one of the scenes to tab through.