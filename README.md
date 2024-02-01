# Image-Editing
電腦圖學導論 Project 1

In this project, I implemented some features such as  **Transform Quantization, Dithering, Filtering and Resizing**.

- **Transform**: Grayscale
- **Quantization**: Minimize the number of colors
  - Uniform Quantization: Minimize the colors evenly
  - Populosity Quantization: Find critical colors and calculate euclidean distance to do minimizing
- **Dithering**: Add noises to reduce inaccuracy caused by quantization
  - Brightness Dithering: Conserve brightness
  - Random Dithering: Pick threshold randomly
  - Clustered Dithering: Apply matrix mask
  - Floyd-Steinberg Dithering: Consider inaccuracy with neighbor pixels
- **Filtering**: Apply different filters to work with each frequencies
  - Gaussian Filter: Remove high frequency
  - Edge Filter: Remove frequencies except edge of object.
 
_IDE: Visual Studio 2019_
