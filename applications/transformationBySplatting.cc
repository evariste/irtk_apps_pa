#include <irtkImage.h>

#include <irtkTransformation.h>

irtkRealImage *source = NULL;
irtkRealImage *output = NULL;
irtkRealImage *target = NULL;
irtkRealImage *nc_n = NULL;
irtkRealImage *nc_c = NULL;

  int target_padding, source_padding;

// Default filenames
char *input_name = NULL, *output_name = NULL, *dof_name  = NULL;

void fill(int x, int y, int z){

  int a, b;
  int radius = 0;

  double val = 0, temp;
  int count = 0;

  while (count == 0){
    ++radius;

    if (x - radius < 0 || x + radius > output->GetX() - 1 ||
        y - radius < 0 || y + radius > output->GetY() - 1 ||
        z - radius < 0 || z + radius > output->GetZ() - 1  ){
      val = 0.0;
      count = 1;
      break;
    }

    for (a =  -radius; a < radius; ++a){
      for (b = -radius; b < radius; ++b){

        temp = output->Get(x + a, y + b, z + radius);
        if (temp > source_padding){
          ++count;
          val += temp;
        }

        temp = output->Get(x + a, y + b, z - radius);
        if (temp > source_padding){
          ++count;
          val += temp;
        }

        temp = output->Get(x + a, y + radius, z + b);
        if (temp > source_padding){
          ++count;
          val += temp;
        }

        temp = output->Get(x + a, y - radius, z + b);
        if (temp > source_padding){
          ++count;
          val += temp;
        }

        temp = output->Get(x + radius, y + a, z + b);
        if (temp > source_padding){
          ++count;
          val += temp;
        }

        temp = output->Get(x - radius, y + a, z + b);
        if (temp > source_padding){
          ++count;
          val += temp;
        }
      }
    }
  }

  output->Put(x, y, z, val / (double) count);

}


void usage()
{
  cerr << "Usage: transformation2 [target] [output] <options>\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok; 
  int i, j, k, ii, jj, kk;
  int xdim, ydim, zdim;
  double x, y, z, wx, wy, wz;
  irtkRealPixel *ptr2pixN, *ptr2pixC, *ptr2pixS, *ptr2pixOut;
  int voxels;
  double val, temp;



  irtkTransformation *transformation = NULL;


  // Check command line
  if (argc < 3){
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read image
  cout << "Reading image ... "; cout.flush();
  target = new irtkRealImage(input_name);
  cout << "done" << endl;

  // Other options
  source_padding = MIN_GREY;
  target_padding = MIN_GREY;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dof_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-source") == 0)){
      argc--;
      argv++;
      source = new irtkRealImage(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;
      argv++;
      target_padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sp") == 0)){
      argc--;
      argv++;
      source_padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  // If there is no source image use copy of target image.
  if (source == NULL){
    cout << "Initialising source from target." << endl;
    source = new irtkRealImage(*target);
  }

  output = new irtkRealImage(*source);
  nc_n = new irtkRealImage(*source);
  nc_c = new irtkRealImage(*source);

  voxels = source->GetNumberOfVoxels();
  ptr2pixN = nc_n->GetPointerToVoxels();
  ptr2pixC = nc_c->GetPointerToVoxels();
  ptr2pixOut = output->GetPointerToVoxels();
  for (i = 0; i < voxels; ++i){
    *ptr2pixN = *ptr2pixC = *ptr2pixOut = 0;
    ++ptr2pixN;
    ++ptr2pixC;
    ++ptr2pixOut;
  }

  if (dof_name != NULL){
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    // Create identity transformation
    transformation = new irtkRigidTransformation;
  }


  xdim = source->GetX();
  ydim = source->GetY();
  zdim = source->GetZ();


  for (k = 0; k < target->GetZ(); ++k){
    for (j = 0; j < target->GetY(); ++j){
      for (i = 0; i < target->GetX(); ++i){
        x = i;
        y = j;
        z = k;

        target->ImageToWorld(x, y, z);
        transformation->Transform(x, y, z);
        source->WorldToImage(x, y, z);

        if (x < 0 || x >= xdim - 1 || y < 0 || y >= ydim - 1 || z < 0 || z >= zdim - 1)
          continue;

        val = target->Get(i, j, k);

//         if (val <= target_padding)
//           continue;

        ii = (int) floor(x);
        jj = (int) floor(y);
        kk = (int) floor(z);

        wx = x - ii;
        wy = y - jj;
        wz = z - kk;

        temp = (1 - wx) * (1 - wy) * (1 - wz);
        nc_c->Put(ii, jj, kk, nc_c ->Get(ii, jj, kk) + temp);
        nc_n->Put(ii, jj, kk, nc_n ->Get(ii, jj, kk) + temp * val);

        temp = (1 - wx) * (wy) * (1 - wz);
        nc_c->Put(ii, jj + 1, kk, nc_c ->Get(ii, jj + 1, kk) + temp);
        nc_n->Put(ii, jj + 1, kk, nc_n ->Get(ii, jj + 1, kk) + temp * val);

        temp = (wx) * (1 - wy) * (1 - wz);
        nc_c->Put(ii + 1, jj, kk, nc_c ->Get(ii + 1, jj, kk) + temp);
        nc_n->Put(ii + 1, jj, kk, nc_n ->Get(ii + 1, jj, kk) + temp * val);

        temp = (wx) * (wy) * (1 - wz);
        nc_c->Put(ii + 1, jj + 1, kk, nc_c ->Get(ii + 1, jj + 1, kk) + temp);
        nc_n->Put(ii + 1, jj + 1, kk, nc_n ->Get(ii + 1, jj + 1, kk) + temp * val);


        temp = (1 - wx) * (1 - wy) * (wz);
        nc_c->Put(ii, jj, kk + 1, nc_c ->Get(ii, jj, kk + 1) + temp);
        nc_n->Put(ii, jj, kk + 1, nc_n ->Get(ii, jj, kk + 1) + temp * val);

        temp = (1 - wx) * (wy) * (wz);
        nc_c->Put(ii, jj + 1, kk + 1, nc_c ->Get(ii, jj + 1, kk + 1) + temp);
        nc_n->Put(ii, jj + 1, kk + 1, nc_n ->Get(ii, jj + 1, kk + 1) + temp * val);

        temp = (wx) * (1 - wy) * (wz);
        nc_c->Put(ii + 1, jj, kk + 1, nc_c ->Get(ii + 1, jj, kk + 1) + temp);
        nc_n->Put(ii + 1, jj, kk + 1, nc_n ->Get(ii + 1, jj, kk + 1) + temp * val);

        temp = (wx) * (wy) * (wz);
        nc_c->Put(ii + 1, jj + 1, kk + 1, nc_c ->Get(ii + 1, jj + 1, kk + 1) + temp);
        nc_n->Put(ii + 1, jj + 1, kk + 1, nc_n ->Get(ii + 1, jj + 1, kk + 1) + temp * val);

      }
    }
  }

  ptr2pixN = nc_n->GetPointerToVoxels();
  ptr2pixC = nc_c->GetPointerToVoxels();
  ptr2pixOut = output->GetPointerToVoxels();
  ptr2pixS = source->GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){

    if (*ptr2pixC > 0){
      *ptr2pixOut = *ptr2pixN / *ptr2pixC;
    } else if (*ptr2pixS <= source_padding){
       *ptr2pixOut = source_padding;
       *ptr2pixC = 1;
    }

    ++ptr2pixN;
    ++ptr2pixC;
    ++ptr2pixOut;
    ++ptr2pixS;
  }


  // But there still may be holes in the output.  Need to do
  // something about these.
  for (k = 0; k < source->GetZ(); ++k){
    for (j = 0; j < source->GetY(); ++j){
      for (i = 0; i < source->GetX(); ++i){

        if (nc_c->Get(i, j, k) <= 0 && source->Get(i, j, k) > source_padding){
          fill(i, j, k);
        }

      }
    }
  }



  // Write the final transformation estimate
  output->Write(output_name);

//   nc_c->Write("nc_c.nii.gz");

}
