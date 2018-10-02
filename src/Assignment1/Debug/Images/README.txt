Brighten:
	./Assignment1 -input sample_images/Jellyfish.jpg -brightness 2 -output brightness.png
Random noise:
	./Assignment1 -input sample_images/Jellyfish.jpg -noise 20 -output noise.png
Crop:
	./Assignment1 -input sample_images/Jellyfish.jpg -crop 500 200 300 367 -output crop.png
Extract Channel:
	./Assignment1 -input sample_images/tree_dark.jpg -extractChannel 0 -output extract.png
Contrast:
	./Assignment1 -input sample_images/cube.jpg -contrast 4 -output contrast.png
Saturation:
	./Assignment1 -input sample_images/Jellyfish.jpg -saturation -1 -output saturation.png
Sharpen:
	./Assignment1 -input sample_images/Jellyfish.jpg -sharpen 4 -output sharpen.png
Quantize:
	./Assignment1 -input sample_images/Jellyfish.jpg -quantize 1 -output quantize.png
Random Dither:
	./Assignment1 -input sample_images/Jellyfish.jpg -randomDither 1 -output randomDither.png
Blur:
	./Assignment1 -input sample_images/Jellyfish.jpg -blur 4 -output blur.png
Edge Detect:
	./Assignment1 -input sample_images/tree_dark.jpg -edgeDetect -output edgeDetect.png
Floyd-Steinberg Dither:
	./Assignment1 -input sample_images/cube.jpg -FloydSteinbergDither 1 -output floyd.png
Scale:
	./Assignment1 -input sample_images/cube.jpg -sampling 0 -scale 3 3 -output scale_cp.png
	./Assignment1 -input sample_images/cube.jpg -sampling 1 -scale 3 3 -output scale_bi.png
	./Assignment1 -input sample_images/cube.jpg -sampling 2 -scale 3 3 -output scale_ga.png
Rotate:
	./Assignment1 -input sample_images/Jellyfish.jpg -sampling 0 -rotate 1 -output rotate.png
	./Assignment1 -input sample_images/Jellyfish.jpg -sampling 2 -rotate 1 -output rotate.png
Fun:
	./Assignment1 -input sample_images/millcity.jpg -sampling 0 -fun -output fun.png
	./Assignment1 -input sample_images/millcity.jpg -sampling 1 -fun -output fun.png
	./Assignment1 -input sample_images/millcity.jpg -sampling 2 -fun -output fun.png
Nonphotorealism:
	./Assignment1 -input sample_images/tree_dark.jpg -nonPhotorealistic -output out.png
Ordered Dither:
	TODO
Lossy Compression:
	IMPROVE, use gaussian inverse
	