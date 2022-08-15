#ifndef CPP_PROJECT_INCLUDE_NBODY_CONFIGURATION_CONFIGURATION_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_CONFIGURATION_CONFIGURATION_HPP_

#include <NBody/Configuration/ColorGenerator.hpp>
#include <NBody/Simulation/Configuration.hpp>
#include <NBody/Simulation/InitialData.hpp>
#include <NBody/Visualization/Configuration.hpp>
#include <boost/program_options.hpp>
#include <concepts>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>

namespace NBody::Configuration {

/// \brief Input mode of the nbody program.
enum class InputMode {
  /// \brief The nbody program will read simulation result from input in Frame
  /// format and then visualize it.
  Frame,
  /// \brief The nbody program will read the initial data from input and then
  /// output the result according to the output mode.
  /// \details The format of initial data:\n
  /// \n
  /// BoundaryBox\n
  /// NumParticles\n
  /// ParticleIndex ParticleMass ParticleRadius ParticlePosition
  /// ParticleVelocity\n
  /// \n
  /// For example:
  /// \n
  /// -500 -500 -500 500 500 500\n
  /// 3\n
  /// 0 100000 15 -199.324 445.806 483.978 -6.2062 0.299517 -2.03983\n
  /// 1 100000 15 -455.226 -280.195 291.115 6.19133 1.83838 0.234251\n
  /// 2 100000 15 136.894 -483.306 148.551 -1.47898 7.98995 3.05997\n
  /// \note Although the boundary will be ignored if the boundary condition
  /// handler type is None, the camera speed in visualizer depends on the value
  /// of boundary box. Therefore, the position of camera will be fixed if the
  /// value of boundary box is [[0, 0, 0], [0, 0, 0]].
  /// \note All particles should locate inside the boundary box and do not
  /// collided with others.
  Simulation
};

/// \brief The output mode of the nbody program.
enum class OutputMode {
  /// \brief The nbody program will output the simulation result in Step format.
  /// \details The Step format:\n
  /// \n
  /// TimeStep\n
  /// StepIndex\n
  /// NumParticles\n
  /// ParticleIndex ParticleMass ParticleRadius ParticlePosition
  /// ParticleVelocity\n
  /// ...\n
  /// StepIndex\n
  /// NumParticles\n
  /// ParticleIndex ParticleMass ParticleRadius ParticlePosition
  /// ParticleVelocity\n
  /// \n
  /// Example:\n
  /// \n
  /// 0.05\n
  /// 0\n
  /// 3\n
  /// 0 1.000000e+02 5.000000e+00 -2.387560e+02 -4.694660e+02 -2.993560e+02
  /// 0.000000e+00 0.000000e+00 0.000000e+00\n
  /// 1 1.000000e+02 5.000000e+00 -2.205280e+01 -4.235270e+02 -5.363880e+01
  /// 0.000000e+00 0.000000e+00 0.000000e+00\n
  /// 2 1.000000e+02 5.000000e+00 -4.040610e+02 -3.761490e+02 -3.154960e+02
  /// 0.000000e+00 0.000000e+00 0.000000e+00\n
  /// 1\n
  /// 3\n
  /// 0 1.000000e+02 5.000000e+00 -2.387560e+02 -4.694660e+02 -2.993560e+02
  /// -1.692612e-04 1.039921e-03 1.822790e-04\n
  /// 1 1.000000e+02 5.000000e+00 -2.205280e+01 -4.235270e+02 -5.363881e+01
  /// -1.984254e-04 3.012267e-04 -1.106167e-04\n
  /// 2 1.000000e+02 5.000000e+00 -4.040610e+02 -3.761490e+02 -3.154960e+02
  /// 8.364275e-04 2.827268e-04 3.373682e-04
  Step,
  /// \brief The nbody program will output the simulation result in Frame
  /// format.
  /// \details The Frame format:\n
  /// \n
  /// FPS\n
  /// Ratio\n
  /// BoundaryBox\n
  /// NumParticles\n
  /// ParticleIndex ParticleRadius ParticlePosition\n
  /// ...\n
  /// NumParticles\n
  /// ParticleIndex ParticleRadius ParticlePosition\n
  /// \n
  /// Example:\n
  /// \n
  /// 60\n
  /// 20\n
  /// -1200 -1200 -1200 1200 1200 1200\n
  /// 3\n
  /// 0 5.000000000000e+00 -2.387560000000e+02 -4.694660000000e+02
  /// -2.993560000000e+02\n
  /// 1 5.000000000000e+00 -2.205280000000e+01 -4.235270000000e+02
  /// -5.363880000000e+01\n
  /// 2 5.000000000000e+00 -4.040610000000e+02 -3.761490000000e+02
  /// -3.154960000000e+02\n
  /// 3\n
  /// 0 5.000000000000e+00 -2.387561523355e+02 -4.694650640701e+02
  /// -2.993558359486e+02\n
  /// 1 5.000000000000e+00 -2.205297858293e+01 -4.235267288959e+02
  /// -5.363889955513e+01\n
  /// 2 5.000000000000e+00 -4.040602472142e+02 -3.761487455457e+02
  /// -3.154956963683e+02
  Frame,
  /// \brief The nbody program will visualize the simulation result.
  Visualization
};

/// \brief Parses the command line options.
template<std::floating_point T>
class Parser {
 public:
  using ValueType = T;
  using SimulationConfiguration = NBody::Simulation::Configuration<ValueType>;
  using VisualizationConfiguration = NBody::Visualization::Configuration;
  using Particle = NBody::Simulation::Particle<ValueType>;

 private:
  boost::program_options::variables_map vm_;

 public:
  /// \brief Constructs a parser with command line options.
  Parser(int argc, char* argv[]) : vm_{} {
    namespace bpo = boost::program_options;
    namespace NS = NBody::Simulation;

    auto desc = bpo::options_description{"Options"};
    desc.add_options()
        // Input Mode
        ("input-mode",
         bpo::value<InputMode>()->default_value(InputMode::Simulation))

        // Input File
        ("input-file,I", bpo::value<std::string>())

        // Output Mode
        ("output-mode",
         bpo::value<OutputMode>()->default_value(OutputMode::Visualization))

        // Output Data Path
        ("output-file,O", bpo::value<std::string>())

        // Output Precision
        ("output-precision", bpo::value<int>()
            ->default_value(6)
            ->notifier(makeRangeValidator(
                0,
                std::numeric_limits<ValueType>::max_digits10,
                "output-precision"
            )))

        // Number of Threads
        ("num-threads", bpo::value<int>()
            ->default_value(1)
            ->notifier(makeRangeValidator(
                1,
                std::numeric_limits<int>::max(),
                "num-threads"
            )))

        // Gravity Solver
        ("gravity-solver", bpo::value<NS::GravitySolverType>()
            ->default_value(NS::GravitySolverType::Direct))

        // Time Step
        ("time-step", bpo::value<ValueType>()
            ->default_value(ValueType{0.05})
            ->notifier(makeRangeValidator(
                std::nextafter(ValueType{0}, ValueType{1}),
                std::numeric_limits<ValueType>::max(),
                "time-step"
            )))

        // G
        ("G", bpo::value<ValueType>()
            ->default_value(ValueType{1})
            ->notifier(makeRangeValidator(
                std::nextafter(ValueType{0}, ValueType{1}),
                std::numeric_limits<ValueType>::max(),
                "G"
            )))

        // Softening
        ("softening", bpo::value<ValueType>()
            ->default_value(ValueType{0})
            ->notifier(makeRangeValidator(
                ValueType{0},
                std::numeric_limits<ValueType>::max(),
                "softening"
            )))

        // Opening Angle
        ("opening-angle", bpo::value<ValueType>()
            ->default_value(ValueType{1})
            ->notifier(makeRangeValidator(
                ValueType{0},
                std::numeric_limits<ValueType>::max(),
                "opening-angle"
            )))

        // Integrator
        ("integrator", bpo::value<NS::IntegratorType>()
            ->default_value(NS::IntegratorType::Leapfrog))

        // Finish Step
        ("finish-step", bpo::value<std::uint64_t>()->default_value(0))

        // Collision Detector
        ("collision-detector", bpo::value<NS::CollisionDetectorType>()
            ->default_value(NS::CollisionDetectorType::None))

        // Collision Resolver
        ("collision-resolver", bpo::value<NS::CollisionResolverType>()
            ->default_value(NS::CollisionResolverType::Hard))

        // Coefficient of Restitution
        ("coefficient-of-restitution,COR", bpo::value<ValueType>()
            ->default_value(ValueType{1})
            ->notifier(makeRangeValidator(
                ValueType{0},
                ValueType{1},
                "coefficient-of-restitution"
            )))

        // Merge Threshold
        ("merge-threshold", bpo::value<std::uint64_t>()
            ->default_value(0))

        // Boundary Condition Handler
        ("boundary-condition-handler",
         bpo::value<NS::BoundaryConditionHandlerType>()
             ->default_value(NS::BoundaryConditionHandlerType::None))

        // Coefficient of Restitution for Boundary
        ("boundary-coefficient-of-restitution,BoundaryCOR",
         bpo::value<ValueType>()
             ->default_value(ValueType{1})
             ->notifier(makeRangeValidator(
                 ValueType{0},
                 ValueType{1},
                 "boundary-coefficient-of-restitution"
             )))

        // Ratio
        ("ratio", bpo::value<ValueType>()
            ->default_value(ValueType{1})
            ->notifier(makeRangeValidator(
                std::nextafter(ValueType{0}, ValueType{1}),
                std::numeric_limits<ValueType>::max(),
                "ratio"
            )))

        // Camera Position Hint
        ("camera-position-hint", bpo::value<ValueType>()
            ->default_value(0)
            ->notifier(makeRangeValidator(
                ValueType{0},
                std::numeric_limits<ValueType>::max(),
                "camera-position-hint"
            )))

        // FPS
        ("FPS", bpo::value<std::uint64_t>()
            ->default_value(60)
            ->notifier(makeRangeValidator(
                1,
                std::numeric_limits<std::uint64_t>::max(),
                "FPS"
            )))

        // Window Width
        ("window-width", bpo::value<int>()
            ->default_value(1024)
            ->notifier(makeRangeValidator(
                1,
                std::numeric_limits<int>::max(),
                "window-width"
            )))

        // Window Height
        ("window-height", bpo::value<int>()
            ->default_value(768)
            ->notifier(makeRangeValidator(
                1,
                std::numeric_limits<int>::max(),
                "window-height"
            )))

        // Icosphere Subdivision
        ("icosphere-subdivision", bpo::value<std::uint64_t>()
            ->default_value(3))

        // Colorful
        ("colorful", bpo::bool_switch()->default_value(false))

        // Preloading
        ("preloading", bpo::bool_switch()->default_value(false))

        // Help
        ("help,h", "Print usage information and exit");

    bpo::store(bpo::parse_command_line(argc, argv, desc), vm_);
    bpo::notify(vm_);

    // Prints help information and then exit.
    if (vm_.contains("help")) {
      printHelpInfo();
      std::exit(0);
    }
  }

  /// \brief Returns the input mode of the nbody program.
  /// \return Input mode
  [[nodiscard]] InputMode getInputMode() const {
    return getOptionValue<InputMode>("input-mode");
  }

  /// \brief Returns the optional input data path.
  /// \details The function will return an empty optional value if input data
  /// path is not specified in command lint options.
  /// \return The optional input data path
  [[nodiscard]] std::optional<std::string> getInputDataPath() const {
    if (vm_.contains("input-file")) {
      return std::make_optional(getOptionValue<std::string>("input-file"));
    } else {
      return std::optional<std::string>{};
    }
  }

  /// \brief Returns the output mode of the nbody program.
  /// \return Output mode
  [[nodiscard]] OutputMode getOutputMode() const {
    return getOptionValue<OutputMode>("output-mode");
  }

  /// \brief Returns the optional output data path.
  /// \details The function will return an empty optional value if output data
  /// path is not specified in command lint options.
  /// \return The optional output data path
  [[nodiscard]] std::optional<std::string> getOutputDataPath() const {
    if (vm_.contains("output-file")) {
      return std::make_optional(getOptionValue<std::string>("output-file"));
    } else {
      return {};
    }
  }

  /// \brief Returns the precision of floating-point number in output.
  /// \return Precision of floating-point number in output
  [[nodiscard]] int getOutputPrecision() const {
    return getOptionValue<int>("output-precision");
  }

  /// \brief Returns the number of threads will be used in this program.
  /// \return The number of threads
  [[nodiscard]] int getNumThreads() const {
    return getOptionValue<int>("num-threads");
  }

  /// \brief Returns the ratio of visualization time to real time.
  /// \return The ratio of visualization time to real time
  [[nodiscard]] ValueType getRatio() const {
    return getOptionValue<ValueType>("ratio");
  }

  /// \brief Returns the FPS of visualization.
  /// \return FPS
  [[nodiscard]] std::uint64_t getFps() const {
    return getOptionValue<std::uint64_t>("FPS");
  }

  /// \brief Returns whether particles are colorful in visualization.
  /// \return true if particles are colorful, false otherwise
  [[nodiscard]] bool getColorful() const {
    return getOptionValue<bool>("colorful");
  }

  /// \brief Returns whether preloading is enabled.
  /// \return true if preloading is enabled, false otherwise
  [[nodiscard]] bool getPreloading() const {
    return getOptionValue<bool>("preloading");
  }

  /// \brief Returns the simulation configuration.
  /// \return Simulation configuration.
  [[nodiscard]] SimulationConfiguration getSimulationConfiguration() const {
    namespace NS = NBody::Simulation;
    auto config = SimulationConfiguration{};
    config.timeStep = getOptionValue<ValueType>("time-step");
    config.finishStep = getOptionValue<std::uint64_t>("finish-step");
    config.G = getOptionValue<ValueType>("G");
    config.softening = getOptionValue<ValueType>("softening");
    config.openingAngle = getOptionValue<ValueType>("opening-angle");
    config.COR = getOptionValue<ValueType>("coefficient-of-restitution");
    config.boundaryCOR =
        getOptionValue<ValueType>("boundary-coefficient-of-restitution");
    config.mergeThreshold = getOptionValue<std::uint64_t>("merge-threshold");
    config.gravitySolverType =
        getOptionValue<NS::GravitySolverType>("gravity-solver");
    config.integratorType = getOptionValue<NS::IntegratorType>("integrator");
    config.boundaryConditionHandlerType =
        getOptionValue<NS::BoundaryConditionHandlerType>(
            "boundary-condition-handler"
        );
    config.collisionDetectorType = getOptionValue<NS::CollisionDetectorType>(
        "collision-detector"
    );
    config.collisionResolverType = getOptionValue<NS::CollisionResolverType>(
        "collision-resolver"
    );
    return config;
  }

  /// \brief Returns the visualization configuration.
  /// \return Visualization configuration.
  [[nodiscard]] VisualizationConfiguration
  getVisualizationConfiguration() const {
    auto config = VisualizationConfiguration{};
    config.cameraPositionHint
        = getOptionValue<ValueType>("camera-position-hint");
    config.fps = getOptionValue<std::uint64_t>("FPS");
    config.ratio = getOptionValue<ValueType>("ratio");
    config.windowWidth = getOptionValue<int>("window-width");
    config.windowHeight = getOptionValue<int>("window-height");
    config.icosphereSubdivision =
        getOptionValue<std::uint64_t>("icosphere-subdivision");
    return config;
  }

 private:
  static auto makeRangeValidator(
      const auto& min,
      const auto& max,
      const std::string& option
  ) {
    return [min, max, option](const auto& value) {
      if (value < min || value > max) {
        throw boost::program_options::validation_error{
            boost::program_options::validation_error::invalid_option_value,
            option,
            std::to_string(value)
        };
      }
    };
  }

  template<typename OptionValueType>
  [[nodiscard]] OptionValueType getOptionValue(const std::string& key) const {
    return vm_[key].template as<OptionValueType>();
  }

  void printHelpInfo() const {
    std::cout << "Usage:\n"
              << "\n"
              << "  nbody [options]\n";

    std::cout << "\n"
              << "Simulation options:\n"
              << "\n";

    // Input mode
    std::cout << std::left << std::setw(40)
              << "  --input-mode <mode>"
              << "Specify input mode\n";

    // Input file
    std::cout << std::left << std::setw(40)
              << "  -I [ --input-file ] <file>"
              << "Read input from <file>\n";

    // Output mode
    std::cout << std::left << std::setw(40)
              << "  --output-mode <mode>"
              << "Specify output mode\n";

    // Output file
    std::cout << std::left << std::setw(40)
              << "  -O [ --output-file ] <file>"
              << "Write output to <file>\n";

    // Output Precision
    std::cout << std::left << std::setw(40)
              << "  --output-precision <value>"
              << "Specify output precision\n";

    // Number of Threads
    std::cout << std::left << std::setw(40)
              << "  --num-threads <value>"
              << "Specify the number of threads will be used in parallel "
                 "computing\n";

    // Gravity solver
    std::cout << std::left << std::setw(40)
              << "  --gravity-solver <type>"
              << "Specify the type of gravity solver\n";

    // Time Step
    std::cout << std::left << std::setw(40)
              << "  --time-step <value>"
              << "Specify time step of integration\n";

    // G
    std::cout << std::left << std::setw(40)
              << "  --G <value>"
              << "Specify gravitational constant for interaction computation\n";

    // Softening
    std::cout << std::left << std::setw(40)
              << "  --softening <value>"
              << "Specify gravitational softening for interaction computation"
                 "\n";

    // Opening Angle
    std::cout << std::left << std::setw(40)
              << "  --opening-angle <value>"
              << "Specify critical opening angle for Barnes-Hut gravity solver"
                 "\n";

    // Integrator
    std::cout << std::left << std::setw(40)
              << "  --integrator <type>"
              << "Specify the type of integrator\n";

    // Finish Step
    std::cout << std::left << std::setw(40)
              << "  --finish-step <value>"
              << "Specify the finish step of simulation\n";

    // Collision Detector
    std::cout << std::left << std::setw(40)
              << "  --collision-detector <type>"
              << "Specify the type of collision detector\n";

    // Collision Resolver
    std::cout << std::left << std::setw(40)
              << "  --collision-resolver <type>"
              << "Specify the type of collision resolver\n";

    // Coefficient of Restitution
    std::cout << std::left << std::setw(40)
              << "  --COR [ --coefficient-of-restitution ] <value>"
              << "Specify the coefficient of restitution for hard collision "
                 "resolver\n";

    // Merge Threshold
    std::cout << std::left << std::setw(40)
              << "  --merge-threshold <value>"
              << "Specify the merge threshold for Hard collision resolver\n";

    // Boundary Condition Handler
    std::cout << std::left << std::setw(40)
              << "  --boundary-condition-handler <type>"
              << "Specify the type of boundary condition handler\n";

    // Coefficient of Restitution for Boundary
    std::cout << std::left << std::setw(40)
              << "  --BoundaryCOR [ --boundary-coefficient-of-restitution ]"
                 " <value>\n"
              << std::left << std::setw(40)
              << " "
              << "Specify the coefficient of restitution for Open boundary "
                 "condition handler \n";

    std::cout << "\n"
              << "Visualization options:\n"
              << "\n";

    // Ratio
    std::cout << std::left << std::setw(40)
              << "  --ratio <value>"
              << "Specify visualization speedup ratio\n";

    // Camera Position Hint
    std::cout << std::left << std::setw(40)
              << "  --camera-position-hint <value>"
              << "Specify initial z-axis position of camera\n";

    // FPS
    std::cout << std::left << std::setw(40)
              << "  --FPS <value>"
              << "Specify maximum FPS\n";

    // Window Width
    std::cout << std::left << std::setw(40)
              << "  --window-width <value>"
              << "Specify initial window width\n";

    // Window Height
    std::cout << std::left << std::setw(40)
              << "  --window-height <value>"
              << "Specify initial window height\n";

    // Icosphere Subdivision
    std::cout << std::left << std::setw(40)
              << "  --icosphere-subdivision <level>"
              << "Specify the subdivision level of icosphere\n";

    // Colorful
    std::cout << std::left << std::setw(40)
              << "  --colorful"
              << "Enable colorful particle\n";

    // Preloading
    std::cout << std::left << std::setw(40)
              << "  --preloading"
              << "Enable input data preloading for Frame input mode\n";

    std::cout << "\n"
              << "Help options:\n"
              << "\n";

    // Help
    std::cout << std::left << std::setw(40)
              << "  -h [ --help ]"
              << "Print help information and exit\n"
              << "\n";

    std::cout << "For more information about options (e.g., range and default "
                 "value), please refer to documentation.\n"
              << "\n";
  }
};

/// \brief Reads the input mode from a character input stream.
/// \param istream The character input stream
/// \param mode The input mode to be inserted
/// \return The character output stream
std::istream& operator>>(std::istream& istream, InputMode& mode) {
  std::string typeStr;
  if (!(istream >> typeStr)) { return istream; }
  if (typeStr == "Frame") {
    mode = InputMode::Frame;
  } else if (typeStr == "Visualization") {
    mode = InputMode::Simulation;
  } else {
    istream.setstate(std::ios_base::failbit);
  }
  return istream;
}

/// \brief Writes the input mode to a character output stream.
/// \param ostream The character output stream
/// \param mode The input mode  to be extracted
/// \return The character output stream
std::ostream& operator<<(std::ostream& ostream, const InputMode& mode) {
  switch (mode) {
    case InputMode::Frame:
      ostream << "Frame";
      break;
    case InputMode::Simulation:
      ostream << "Simulation";
      break;
  }
  return ostream;
}

/// \brief Reads the output mode from a character input stream.
/// \param istream The character input stream
/// \param mode The output mode to be inserted
/// \return The character output stream
std::istream& operator>>(std::istream& istream, OutputMode& mode) {
  std::string typeStr;
  if (!(istream >> typeStr)) { return istream; }
  if (typeStr == "Step") {
    mode = OutputMode::Step;
  } else if (typeStr == "Frame") {
    mode = OutputMode::Frame;
  } else if (typeStr == "Visualization") {
    mode = OutputMode::Visualization;
  } else {
    istream.setstate(std::ios_base::failbit);
  }
  return istream;
}

/// \brief Writes the output mode to a character output stream.
/// \param ostream The character output stream
/// \param mode The output mode to be extracted
/// \return The character output stream
std::ostream& operator<<(std::ostream& ostream, const OutputMode& mode) {
  switch (mode) {
    case OutputMode::Step:
      ostream << "Step";
      break;
    case OutputMode::Frame:
      ostream << "Frame";
      break;
    case OutputMode::Visualization:
      ostream << "Visualization";
      break;
  }
  return ostream;
}

} // namespace NBody::Configuration

#endif //CPP_PROJECT_INCLUDE_NBODY_CONFIGURATION_CONFIGURATION_HPP_
