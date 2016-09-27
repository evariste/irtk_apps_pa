set -vx

# Run like this: bash TESTS.sh  &> test-output.txt

# Then can compare test output with previously saved output in
# test-output-archive.txt


tar xf outputFilesArchive.tar  

for f in `ls outputFilesArchive`; do rm -f $f; done


polydatasphere bla.vtk -solid 0

polydatasphere bla2.vtk -solid 1

polydataappend  2 bla.vtk bla2.vtk  bla3.vtk

polydatabumps  bumps.vtk  -l 4  -randomRotation

polydatacurvatureindices  bumps.vtk 


polydatadecimate bumps.vtk bumps-d.vtk

polydatalistscalars bumps.vtk 

polydatadeletescalars bumps.vtk bumps.vtk -name theta

polydatalistscalars bumps.vtk 



binarize tr-aal-to-neo-nr-40-crop.nii.gz  temp.nii.gz 1 1000  1000 0

mcubes temp.nii.gz temp.vtk 500

polydata2maskimage  temp.vtk temp.nii.gz  temp-recovered.nii.gz  -value 10

polydataassignlabels  tr-aal-to-neo-nr-40-crop.nii.gz  temp.vtk    temp-label.vtk -name aal



extract_roi tr-aal-to-neo-nr-40-crop.nii.gz  rois.nii.gz  10 20

rescale rois.nii.gz  rois.nii.gz  0 1000

dilation  rois.nii.gz rois.nii.gz 

mcubes  rois.nii.gz rois.vtk 500 

polydatalcc  rois.vtk  rois-lcc.vtk   


polydatamaths bumps.vtk -array_name phi -mul 5 -add 2 bumps-maths.vtk


polydatascalarstats  bumps.vtk -name phi

# Scalar name   phi
# No of pts     40962
# After masking 40962
# Mean          3.1317
# Mean Sq       13.0942
# S.D.          1.8129
# Mean(abs)     3.1317
# S.D(abs)      1.8129
# Min/Max       0 6.27113
# Area          22.0564
# Area int      69.1771
# sqrt(Area int / 4 pi) 2.34626

polydatascalarstats  bumps-maths.vtk -name phi

# Scalar name   phi
# No of pts     40962
# After masking 40962
# Mean          17.6585
# Mean Sq       393.988
# S.D.          9.06452
# Mean(abs)     17.6585
# S.D(abs)      9.06452
# Min/Max       2 33.3556
# Area          22.0564
# Area int      389.998
# sqrt(Area int / 4 pi) 5.57091


polydatarecalculatenormals bumps.vtk  bumps-recalcN.vtk



polydatareflect bumps.vtk  bumps-r.vtk -x 



polydatadecimate bumps.vtk bumps-dec.vtk  -reduction 0.7

polydataremesh bumps-dec.vtk bumps.vtk bumps-remesh.vtk  

# Two runs of random scalars starting with the same polydata.
polydatarandomscalars  bumps.vtk bumps-rs.vtk

polydatarandomscalars  bumps.vtk bumps-rs-2.vtk

# The stats should be different for each run of the above.
polydatascalarstats  bumps-rs.vtk -name Random

polydatascalarstats  bumps-rs-2.vtk -name Random


polydatascalarsmooth bumps-rs.vtk bumps-ss.vtk 20  2   -name Random



polydatascalarstats  bumps-rs.vtk -name Random

# Scalar name   Random
# No of pts     40962
# After masking 40962
# Mean          0.498624
# Mean Sq       0.33185
# S.D.          0.288486
# Mean(abs)     0.498624
# S.D(abs)      0.288486
# Min/Max       3.4892e-06 0.999979
# Rob Min/Max   0.98999 0.0103589
# Area          22.0564
# Area int      10.9822
# sqrt(Area int / 4 pi) 0.934846


polydatascalarstats  bumps-ss.vtk -name Random
# Scalar name   Random
# No of pts     40962
# After masking 40962
# Mean          0.498257
# Mean Sq       0.248759
# S.D.          0.022353
# Mean(abs)     0.498257
# S.D(abs)      0.022353
# Min/Max       0.40348 0.570199
# Rob Min/Max   0.546683 0.444611
# Area          22.0564
# Area int      10.9763
# sqrt(Area int / 4 pi) 0.934594



polydatasmooth  bumps.vtk  bumps-smooth.vtk   20 0.7


polydatarandomsphericalpoints ranSph-1.vtk
polydatarandomsphericalpoints ranSph-2.vtk

# Following should contain 1000 points at random on a unit sphere
# each. They should be different slightly, check by looking at centre
# of gravity.
polydatacentreofgravity ranSph-1.vtk 
polydatacentreofgravity ranSph-2.vtk 


# File comparisons

for f in `ls outputFilesArchive`
do
    ls -l $f outputFilesArchive/$f
    echo
done

