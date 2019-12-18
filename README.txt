================== GSModelLibraryProjects Build Instructions ==================

-- Created: 12/17/2019
-- Last Updated: 12/18/2019

The GS Model Library project depends on the GSL (GNU Scientific Library), which can be obtained from this GitHub repository: https://github.com/BrianGladman/gsl. For additional information about the GSL, see https://www.gnu.org/software/gsl/. For information about building the GSL on Windows, see https://www.gnu.org/software/gsl/extras/native_win_builds.html.

Before building GSModelLibraryProjects, carry out the following steps to build the GSL on Windows:

1) Clone https://github.com/BrianGladman/gsl into a folder named 'gsl'. The 'gsl' folder should be at the same level as the folder into which the GSModelLibraryProjects was cloned. In other words, the folder structure should be as follows:

  Main Project Folder
    |
     ------ gsl
    |
     ------ GSModelLibraryProjects (contains GSModelLibraryProjects.sln)

2) Open the gsl.lib.sln file inside gsl\build.vc in Visual Studio. At the time of this writing, this is a VS 2015 .sln file, but it was successfully built in VS 2019.

3) As explained in the gsl.readme.txt file (inside gsl\build.vc in the gsl repository), you must first build the gslhdrs project. Go to Build|Batch Build..., check all four boxes for the gslhdrs project, and click Build. This should finish in 10 - 20 seconds.

4) Now, in Build|Batch Build..., uncheck all boxes, check only the boxes for the cblaslib and gsllib projects (8 in all), and click Build.

The build process for GSL will probably take 5 - 10 minutes. Once it has finished, open the GSModelLibraryProjects.sln file (VS 2019 was used at the time of this writing), and build the projects as follows:

1) At the top of the Visual Studio UI, select either 'Release' or 'Debug' and either x64 or x86.

2) Right-click the solution in the Solution Explorer and choose 'Build Solution'

The 32-bit binaries can be found in the 'Debug' and 'Release' folders at the root level (where the .sln file is located). The 64-bit binaries can be found in the 'Debug' and 'Release' folders inside the x64 subfolder. 

The GSSplineInterpolator.dll and GSTimeSeriesAnalysis.dll files can be used with GoldSim. Copy these to the 'GoldSimModels' subfolder and run the model examples therein to confirm that everything works. Open the External Elements in the models and go to the 'Interface' tab to see the required inputs and the outputs of the DLLs.