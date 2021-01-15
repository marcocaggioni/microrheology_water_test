# microrheology_water_test

This repository provide example data and analysis notebook to perform particle tracking microrheology and differential dynamic microscopy analysis.
The required python environment is provided by the docker image at this link https://hub.docker.com/r/mcaggio/microrheologyserver

the example videos can be downloaded form the release page https://github.com/marcocaggioni/microrheology_water_test/releases

These are microscopy videos of 400nm spherical tracer particles diffusing in water at room temperature
The resolution is 0.1 um per pixel for the uncompressed videos
The frame per second is 100 for both videos

* [400nm_100dil_water_01umpix_100fps_short.cin](https://github.com/marcocaggioni/microrheology_water_test/releases/download/1.0/400nm_100dil_water_01umpix_100fps_short.cin)
uncopressed .cin format from Phantom fast camera - low tracers concentration (allows multi particle tracking microrheology approach)

for a video preview you can download the [compressed format](https://github.com/marcocaggioni/microrheology_water_test/releases/download/1.0/400nm_100dil_water_01umpix_100fps_short.mp4)

![image](https://user-images.githubusercontent.com/5448255/104773646-6131a180-5743-11eb-94b8-48fc7f4c2e96.png)

* [400nm_water_01umpix_100fps_short.cine](https://github.com/marcocaggioni/microrheology_water_test/releases/download/1.0/400nm_water_01umpix_100fps_short.cine)
uncopressed .cin format from Phantom fast camera - high tracers concentration (multi particle tracking microrheology not possible but Differential Dynamic Microscopy still feasable)

for a preview you can download the [compressed format](https://github.com/marcocaggioni/microrheology_water_test/releases/download/1.0/400nm_water_01umpix_100fps_short.mp4)
![image](https://user-images.githubusercontent.com/5448255/104773727-7dcdd980-5743-11eb-9edd-0b360d30c636.png)

The mp4 version of the videos is provided for easy visualization
The .cin format can be read directly from python using the pims library


