## createcolormap.m 
[![View createcolormap on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://jp.mathworks.com/matlabcentral/fileexchange/100089-createcolormap)   
This function allows to create colormap Nx3 array (RGB) with an arbitrary combination of colors. 
RGB values between the specified colors will be smoothly connected by linear interpolation.


<p align="center">
<img src="https://github.com/hydrocoast/createcolormap/blob/main/createcolormap_example.png", width="600">
</p>  
<p align="center">
<img src="https://github.com/hydrocoast/createcolormap/blob/main/createcolormap_example_N5.gif", width="720">
</p>  

## Example
  ```matlab
  cmap = createcolormap(C);
  cmap = createcolormap(n,C);
  cmap = createcolormap(colorA, colorB);
  cmap = createcolormap(n, colorA, colorB);
  cmap = createcolormap(colorA, colorB, colorC, colorD, ...);
  cmap = createcolormap(n, colorA, colorB, colorC, colorD, ...);
  ```
  where `n` is the number of segments for the output color scheme, and `C` is the RGB matrix of color junctions.

## Usage
+ blue-white-red (polar)
  ```matlab
  b = [0,0,1];
  w = [1,1,1];
  r = [1,0,0];

  bwr = createcolormap(b,w,r); % 256x3 array
  
  colormap(bwr)
  colorbar
  ```

  If you want to use dark blue and red colors, try below:
  ```matlab
  b = [0.0,0.0,0.5];
  w = [1.0,1.0,1.0];
  r = [0.5,0.0,0.0];

  bwr = createcolormap(b,w,r); % 256x3 array

  colormap(bwr)
  colorbar
  ```

  To create a more discrete color structure, input the number of elements in the first argument as shown below.
  ```matlab
  bwr = createcolormap(16,b,w,r); % 16x3 array 
  ```

+ more complicated combination
  ```matlab
  colorA = [0.0,1.0,0.0];
  colorB = [1.0,0.5,0.5];
  colorC = [0.5,0.5,0.5];
  colorD = [1.0,1.0,0.0];

  cmap = createcolormap(64,colorA,colorB,colorC,colorD); % 64x3 array

  surf(peaks); 
  colormap(cmap);
  colorbar;
  ```

+ RGB matrix
  ```matlab
  cmap = createcolormap(rand(10,3)); % 10 random colors

  surf(peaks);
  colormap(cmap);
  colorbar;
  ```


## License
MIT

## Author
[Takuya Miyashita](https://hydrocoast.jp)   
Disaster Prevention Research Institute, Kyoto University  

## Update
  v0.1  2021/10/01
  v0.2  2021/10/09

