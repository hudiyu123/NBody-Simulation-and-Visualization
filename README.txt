1. Software Building and Installation

For the purpose of building the program, the CMake tool is used. Since the
program utilize features from C++20 features, the compiler which is used needs
to be compatible with C++20. Note that, only Linux and macOS platform are
supported to build this software.

This software utilizes several libraries other than the C++20 standard library.
Therefore, before building the software, users need to ensure that these
libraries are installed. The versions of the libraries used by this software and
known to work are as below:

Linux

  • OpenMP  version 4.5
  • Boost   version 1.76.0
  • OpenGL  version 2.1
  • GLFW    version 3.2.1


macOS

  • OpenMP  version 4.5
  • Boost   version 1.76.0
  • OpenGL  version 3.3
  • GLFW    version 3.3.4
  • GLEW    version 2.2.0


This program supports three types of precision for floating-point arithmetic.
These types are shown as below:

  • Single   - float in C++
  • Double   - double in C++
  • Extended - long double in C++

This precision is configured by setting Cmake variable SIMULATION_PRECISION. The
default value for this variable is Double. Note that, the actual format depends
on the implementation of hardware.

In what follows, let $SOURCE_DIR denotes the top-level directory of the NBody
Simulation software, $BUILD_DIR denotes a directory where build files are
created, $INSTALL_DIR denotes a directory where the software is to be installed,
and $PRECISION denotes the floating-point precision used for computing by the
software. To build and install the NBody Simulation software, perform the
following steps:

  • Generate the native build files by running the command:

    cmake -H$SOURCE_DIR -B$BUILD_DIR -DSIMULATION_PRECISION=$PRECISION
          -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR


  • Build the software by running the command:

    cmake --build $BUILD_DIR --clean-first


  • Install the software by running the command:

    cmake --build $BUILD_DIR --clean-first --target install



2. Software Running

The nbody program has two input modes and three output modes, these modes are
described as below:

    Input

      • Simulation: Reads the initial data in Initial format for simulation from
                    input.
      • Frame: Reads simulation result in Frame format from input and visualizes
               the result. Only Visualization output mode is supported in Frame
               input mode.

    Output

      • Step: Writes the simulation result in Step format to output.
      • Frame: Writes the simulation result in Frame format to output.
      • Visualization: Visualizes the simulation result.


Options:

For more information about options, please refer to the documentation.

--input-mode $mode
    Set the input mode to $mode. Argument $mode accepts Simulation or Frame.
    The default type for $mode is Simulation.

-I [ --input-file ] $path
    Set file path for input to $path. If this option is not specified, standard
    input is used.

--output-mode $mode
    Set the output mode to $mode. Argument $mode accepts Step, Frame or
    Visualization. The default value for $mode is Visualization.

-O [ --output-file ] $path
    Set file path for output to $path. This option is ignored if the output
    mode is Visualization. If this option is not specified, standard output is
    used.

--output-precision $value
    Set the precision of floating-point number for output to $value. Argument
    $value accepts an integer that is greater than or equal to 0. The default
    value for $value is 6.

--num-threads $num
    Set the number of threads will be used in parallel computing to $num.
    Argument $num accepts an integer that is greater than or equal to 1. The
    default value for $num is 1.

--gravity-solver $type
    Set the type of gravity solver to $type. Argument $type accepts Direct or
    BarnesHut. The default value for $type is Direct.

--time-step $value
    Set the time step of integration to $value. Argument $value accepts a
    floating-point number that is greater than 0. The default value for $value
    is 0.05.

--G $value
    Set the gravitational constant for interaction computation to $value.
    Arguments $value accepts a floating-point number that is greater than 0. The
    default value for $value is 1.

--softening $value
    Set the gravitational softening for interaction computation to $value.
    Arguments $value accepts a floating-point number that is greater than or
    equal to 0. The default value for $value is 0.

--opening-angle $value
    Set the critical opening angle for BarnesHut gravity solver to $value.
    Arguments $value accepts a floating-point number that is greater than or
    equal to 0. The default value for $value is 1.

--integrator $type
    Set the type of integrator to $type. Argument $type only accepts Leapfrog.
    The default value for $type is Leapfrog.

--finish-step $num
    Set the finish step of simulation to $num. Argument $num accepts an integer
    that is greater than or equal to 0. Note that, the simulation will finish
    immediately after start if $num equals 0. The default value for $num is 0.

--collision-detector $type
    Set the type of collision detector to $type. Argument $type accepts Direct
    or Octree. The default value for $type is Direct.

--collision-resolver $type
    Set the type of collision resolver to $type. Argument $type accepts Hard
    or Merge. The default value for $type is Hard.

--COR [ --coefficient-of-restitution ] $value
    Set the coefficient of restitution for Hard collision resolver to $value.
    Argument $value accepts a floating-point number in range [0, 1]. The default
    value for $value is 1. Note that, the Hard collision resolver will perform
    the same as Merge collision resolver if the COR equals 0.

--merge-threshold $value
    Set the merge threshold for Hard collision resolver to $value. Argument
    $value accepts an integer that is greater than or equal to 0. The default
    value for $value 0.

--boundary-condition-handler $type
    Set the type of boundary condition handler to $type. Argument $type accepts
    None, Open, Rebound, and Periodic. The default value for $type is None.

--BoundaryCOR [ --boundary-coefficient-of-restitution ] $value
    Set the boundary coefficient of restitution for Rebound boundary condition
    handler to $value. Argument $value accepts a floating-point number in range
    [1, 0]. The default value for $value is 1. Note that, the Rebound boundary
    condition handler will perform the same as Open boundary condition handler
    if the BoundaryCOR equals 0.

--ratio $value
    Set the ratio of visualizing speed to real time to $value. Argument $value
    accepts a floating-point number that is greater than 0. The default value
    for $value is 1.

--camera-position-hint $value
    Set the initial z-axis position for camera to $value. Argument $value
    accepts a floating-point number that is greater than or equal to 0. The
    default value for $value is 0

--FPS $value
    Set the FPS of visualization to $value. Argument $value accepts an integer
    that is greater than 0. The default value for $value is 60.

--window-width $value
    Set the initial window width to $value. Argument $value accepts an integer
    that is greater than 0. The default value for $value is 1024.

--window-height $value
    Set the initial window height to $value. Argument $value accepts an integer
    that is greater than 0. The default value for $value is 768.

--icosphere-subdivision $level
    Set the level of subdivision for icosphere to $level. Argument $level
    accepts an integer that is greater than or equal to 0. The default value for
    $level is 3.

--colorful
    Enable colorful particles for visualization. If this option is not
    specified, particles will have white color.

--preloading
    Enable preloading of input data for Frame input mode. If this option is not
    specified, input data will be read simultaneously while visualizing.

-h [ --help]
    Prints help information and exit.



3. Graphical User Interface

Key                     Action
Esc                     Stop simulation and close window
Space                   Pause simulation and visualizing
WASD                    Move forward, turn left, move backward, and turn right
↑,↓,→,← (arrow keys)    Pitch up, pitch down, yaw right, yaw left
Shift                   Move or rotate faster
Control                 Move or rotate slower

For example:

    If a camera is at [1, 1, 1] and looks at the origin [0, 0, 0]. The user
    wants to move to [-1, 1, 1] and still looks at the origin, he can achieve
    this by step as below:

      • Pitch up 45 degrees by pressing the ↑ key (now he looks at [0, 1, 0])
      • Yaw left 45 degrees by pressing the ← key (now he looks at [0, 1, 1])
      • Move forward 2 by pressing the W key (now he is at [-1, 1, 1] and looks
        at [-2, 1, 1])
      • Yaw right 135 degrees by pressing the Right key (now he looks at
        [0, 1, 0])
      • Pitch down 45 degrees by pressing the Down key (now he looks at the
        origin)



4. Input Data Format

1) Initial Data Format

BoundaryBox
NumParticles
ParticleIndex ParticleMass ParticleRadius ParticlePosition ParticleVelocity

For example:

    -500 -500 -500 500 500 500
    3
    0 100000 15 -199.324 445.806 483.978 -6.2062 0.299517 -2.03983
    1 100000 15 -455.226 -280.195 291.115 6.19133 1.83838 0.234251
    2 100000 15 136.894 -483.306 148.551 -1.47898 7.98995 3.05997


Note that, all particles must be inside the boundary box and not collided with
others, otherwise the program will fail in input validation process.


2) Frame Format

FPS
Ratio
BoundaryBox
NumParticles
ParticleIndex ParticleRadius ParticlePosition
.
.
.
NumParticles
ParticleIndex ParticleRadius ParticlePosition

For example:

    60
    20
    -1200 -1200 -1200 1200 1200 1200
    3
    0 5.000000e+00 -2.387560e+02 -4.694660e+02 -2.993560e+02
    1 5.000000e+00 -2.205280e+01 -4.235270e+02 -5.363880e+01
    2 5.000000e+00 -4.040610e+02 -3.761490e+02 -3.154960e+02
    3
    0 5.000000e+00 -2.387561e+02 -4.694650e+02 -2.993558e+02
    1 5.000000e+00 -2.205297e+01 -4.235267e+02 -5.363889e+01
    2 5.000000e+00 -4.040602e+02 -3.761487e+02 -3.154956e+02


3) Step Format

TimeStep
StepIndex
NumParticles
ParticleIndex ParticleMass ParticleRadius ParticlePosition ParticleVelocity
.
.
.
StepIndex
NumParticles
ParticleIndex ParticleMass ParticleRadius ParticlePosition ParticleVelocity

For example:

    0.05
    0
    3
    0 1.000000e+02 5.000000e+00 -2.387560e+02 -4.694660e+02 -2.993560e+02 0.000000e+00 0.000000e+00 0.000000e+00
    1 1.000000e+02 5.000000e+00 -2.205280e+01 -4.235270e+02 -5.363880e+01 0.000000e+00 0.000000e+00 0.000000e+00
    2 1.000000e+02 5.000000e+00 -4.040610e+02 -3.761490e+02 -3.154960e+02 0.000000e+00 0.000000e+00 0.000000e+00
    1
    3
    0 1.000000e+02 5.000000e+00 -2.387560e+02 -4.694660e+02 -2.993560e+02 -1.692612e-04 1.039921e-03 1.822790e-04
    1 1.000000e+02 5.000000e+00 -2.205280e+01 -4.235270e+02 -5.363881e+01 -1.984254e-04 3.012267e-04 -1.106167e-04
    2 1.000000e+02 5.000000e+00 -4.040610e+02 -3.761490e+02 -3.154960e+02 8.364275e-04 2.827268e-04 3.373682e-04



5. Input Data File

Three input data files are provided in this repository for test:

1) input/10000_particles_initial_data.txt

This file contains 10000 particles in Initial Data format.

Usage example:

    nbody \
    --input-file input/10000_particles_initial_data.txt \
    --output-mode Visualization \
    --num-threads 16 \
    --finish-step 2000000 \
    --gravity-solver BarnesHut \
    --collision-detector Octree \
    --collision-resolver Merge \
    --ratio 3 \
    --camera-position-hint 20000 \
    --colorful


2) input/demo_initial_data.txt

This file contains 40 particles in Initial Data format.

Usage example:

    nbody \
    --input-file input/demo_initial_data.txt \
    --time-step 0.01 \
    --finish-step 2000000 \
    --collision-detector Direct \
    --ratio 30 \
    --boundary-condition-handler Rebound \
    --camera-position-hint 2000 \
    --colorful


3) input/solar_frame.txt

This file contains the simulation result of solar system in Frame format.

Usage example:

    nbody \
    --input-mode Frame \
    --input-file input/solar_frame.txt \
    --camera-position-hint 25000 \
    --colorful
