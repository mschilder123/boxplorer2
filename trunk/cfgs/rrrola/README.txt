From http://www.ms.mff.cuni.cz/~kadlj3am/big/boxplorer/screenshots/
Some with tweaks to enable DoF and LoD.
Tweaked shader code.

If you want to make a movie starting with one of these configs,
you need to first copy it and the shaders to its own directory.

Something like:

cd cfgs\rrrola
mkdir ..\mymovie
copy keyframe-N.cfg ..\mymovie\default.cfg
copy *.glsl ..\mymovie\
cd ..\..

Now you can start the scene with 'boxplorer cfgs\mymovie\'
and create keyframes with SPACE. Those keyframe .cfg files will be written
to cfgs\mymovie\keyframe-*.cfg.

Once you have your keyframes done, exit boxplorer and
restart with 'boxplorer cfgs\mymovie\ --render'.

It will now go through the keyframes and write screenshot .tga files in cfgs\mymovie\.
Piece the .tga files together with mencoder or picasa or any other timelapse composer.
