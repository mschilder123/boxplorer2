From http://www.ms.mff.cuni.cz/~kadlj3am/big/boxplorer/screenshots/
Some with tweaks to enable DoF and LoD.
Tweaked shader code.

If you want to make a movie starting with one of these configs,
you need to first copy it to the cfgs directory and copy the shaders to
its .cfg.data directory.

Something like:

cd cfgs\rrrola.cfg.data
copy keyframe-N.cfg ..\mymovie.cfg
mkdir ..\mymovie.cfg.data
copy *.glsl ..\mymove.cfg.data
cd ..\..

Now you can start the scene with 'boxplorer cfgs\mymovie.cfg'
and create keyframes with SPACE. Those keyframe .cfg files will be written
to cfgs\mymovie.cfg.data\keyframe-*.cfg.

Once you have your keyframes done, exit boxplorer and
restart with 'boxplorer cfgs\mymovie.cfg --render'.

It will now go through the keyframes and write screenshot .tga files in cfgs\mymove.cfg.data.
Piece the .tga files together with mencoder or picasa or any other timelapse composer.
