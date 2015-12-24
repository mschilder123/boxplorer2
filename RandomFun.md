# Fun things to try #

Conway's Game of Life.

# Details #

Boxplorer2 keeps a backbuffer containing the previous frame. When making that available as texture, it is straightforward to implement Conway's Game of Life (http://en.wikipedia.org/wiki/Conway's_Game_of_Life) in a pixel shader, as shown by syntopia in a Fragmentarium example shader.

Add loading interesting starting positions and it gets entertaining.

For instance, try
> `./boxplorer2 cfgs/syntopia/life/ --lifeform=fermatprimecalculator_106.lif`
and watch the mayhem unfold when the action wraps around the screen.

Upon resize, the lifeform gets restarted. Larger is better.

If you start that shader without a lifeform, you'll get syntopia's Ring o'Fire generating random action. Try swirling it around.