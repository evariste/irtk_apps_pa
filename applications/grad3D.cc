#include <irtkImage.h>
#include <irtkGaussianBlurring.h>

char *input_name = NULL, *output_name = NULL;
char *sepBasename = NULL;

void usage()
{
	cerr << "Usage: grad3D [in] [out] [sigma] <-options> " << endl << endl;

	cerr << "A blur with S.D. sigma is appled before the gradient is estimated" << endl;
	cerr << "using finite differences." << endl;
	cerr << "Options:" << endl;
	cerr << "-sep basename    Write separate components to basename-x.nii.gz," << endl;
	cerr << "                 basename-y.nii.gz, basename-z.nii.gz." << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	bool ok;
	irtkRealImage input;
	irtkRealImage output;
	irtkRealImage gradX;
	irtkRealImage gradY;
	irtkRealImage gradZ;
	double voxelx, voxely, voxelz;
	double dx, dy, dz;
	int i, j, k, xdim, ydim, zdim;
	irtkRealPixel tmp = 0;
	double ssq, sigma;
	irtkMatrix w2i;

	if (argc < 4){
		usage();
	}

	input_name  = argv[1];
	argc--;
	argv++;
	output_name = argv[1];
	argc--;
	argv++;
	sigma = atof(argv[1]);
	argc--;
	argv++;

	while (argc > 1){
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-sep") == 0)){
			argc--;
			argv++;
			sepBasename = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if (ok == false){
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}


	// Read input
	input.Read(input_name);
	output.Read(input_name);
	gradX.Read(input_name);
	gradY.Read(input_name);
	gradZ.Read(input_name);

	// Blur input image
	if (sigma > 0){
		irtkGaussianBlurring<irtkRealPixel> gaussianBlurring( sigma );
		gaussianBlurring.SetInput (&input);
		gaussianBlurring.SetOutput(&input);
		gaussianBlurring.Run();
	}

	w2i = input.GetWorldToImageMatrix();

	xdim = input.GetX();
	ydim = input.GetY();
	zdim = input.GetZ();

	input.GetPixelSize(&voxelx, &voxely, &voxelz);
	voxelx = abs(voxelx);
	voxely = abs(voxely);
	voxelz = abs(voxelz);

	for (k = 0; k < zdim; ++k) {
		for (j = 0; j < ydim; ++j){
			for (i = 0; i < xdim; ++i){

				if (i == 0 || i == xdim - 1 ||
						j == 0 || j == ydim - 1 ||
						k == 0 || k == zdim - 1) {

					gradX.Put(i, j, k, 0);
					gradY.Put(i, j, k, 0);
					gradZ.Put(i, j, k, 0);
					output.Put(i, j, k, 0);
					continue;

				}

				ssq = 0.0;

				// Gradient in image coordinates.
				dx = (input.Get(i+1, j  , k  ) - input.Get(i-1, j  , k  )) / 2.0;
				dy = (input.Get(i  , j+1, k  ) - input.Get(i  , j-1, k  )) / 2.0;
				dz = (input.Get(i  , j  , k+1) - input.Get(i  , j  , k-1)) / 2.0;

				// Convert to world coordinates.
				// Assume scalar function is f
				// World coords are (w1, w2, w3) and image coordinates are (im1, im2, im3).
				// e.g for w1,
				//
				// d f    d f   d im1   d f   d im2   d f   d im3
				// ---- = ----- ----- + ----- ----- + ----- -----
				// d w1   d im1 d w1    d im2 d w1    d im3 d w1
				//
				// I.e. can obtain world gradient by multiplying image gradient by world to image matrix!
				tmp = w2i(0, 0) * dx + w2i(0, 1) * dy + w2i(0, 2) * dz;
				gradX.Put(i, j, k, tmp);
				ssq += tmp * tmp;

				tmp = w2i(1, 0) * dx + w2i(1, 1) * dy + w2i(1, 2) * dz;
				gradY.Put(i, j, k, tmp);
				ssq += tmp * tmp;

				tmp = w2i(2, 0) * dx + w2i(2, 1) * dy + w2i(2, 2) * dz;
				gradZ.Put(i, j, k, tmp);
				ssq += tmp * tmp;

				output.Put(i, j, k, sqrt(ssq));

			}
		}
	}


	output.Write(output_name);

	if (sepBasename != NULL){
		char buffer[250];

		sprintf(buffer, "%s-x.nii.gz", sepBasename);
		gradX.Write(buffer);

		sprintf(buffer, "%s-y.nii.gz", sepBasename);
		gradY.Write(buffer);

		sprintf(buffer, "%s-z.nii.gz", sepBasename);
		gradZ.Write(buffer);

	}



}

