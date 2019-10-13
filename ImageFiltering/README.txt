Daniel Olson: Image Filtering Commands

Brighten:
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -brightness 2 -output brightness.png
Random noise:
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -noise 20 -output noise.png
Crop:
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -crop 500 200 300 367 -output crop.png
Extract Channel:
	./build/Release/ImageFiltering.exe-input sample_images/tree_dark.jpg -extractChannel 0 -output extract.png
Contrast:
	./build/Release/ImageFiltering.exe-input sample_images/cube.jpg -contrast 4 -output contrast.png
Saturation:
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -saturation -1 -output saturation.png
Sharpen:
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -sharpen 4 -output sharpen.png
Quantize:
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -quantize 1 -output quantize.png
Random Dither:
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -randomDither 1 -output randomDither.png
Blur:
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -blur 4 -output blur.png
Edge Detect:
	./build/Release/ImageFiltering.exe-input sample_images/tree_dark.jpg -edgeDetect -output edgeDetect.png
Floyd-Steinberg Dither:
	./build/Release/ImageFiltering.exe-input sample_images/cube.jpg -FloydSteinbergDither 1 -output floyd.png
Scale:
	./build/Release/ImageFiltering.exe-input sample_images/cube.jpg -sampling 0 -scale 3 3 -output scale_cp.png
	./build/Release/ImageFiltering.exe-input sample_images/cube.jpg -sampling 1 -scale 3 3 -output scale_bi.png
	./build/Release/ImageFiltering.exe-input sample_images/cube.jpg -sampling 2 -scale 3 3 -output scale_ga.png
Rotate:
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -sampling 0 -rotate 1 -output rotate_cp.png
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -sampling 1 -rotate 1 -output rotate_bi.png
	./build/Release/ImageFiltering.exe-input sample_images/Jellyfish.jpg -sampling 2 -rotate 1 -output rotate_ga.png
Fun:
	./build/Release/ImageFiltering.exe-input sample_images/cube.jpg -sampling 0 -fun -output fun_cp.png
	./build/Release/ImageFiltering.exe-input sample_images/cube.jpg -sampling 1 -fun -output fun_bi.png
	./build/Release/ImageFiltering.exe-input sample_images/cube.jpg -sampling 2 -fun -output fun_ga.png
Nonphotorealism:
	./build/Release/ImageFiltering.exe-input sample_images/tree_dark.jpg -nonPhotorealistic -output out.png
	