#ifndef CPP_PROJECT_INCLUDE_NBODY_NBODYPROGRAM_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_NBODYPROGRAM_HPP_

#include <NBody/Concurrency/Queue.hpp>
#include <NBody/Configuration/Parser.hpp>
#include <NBody/Configuration/ColorGenerator.hpp>
#include <NBody/Simulation/Simulator.hpp>
#include <NBody/Visualization/Visualizer.hpp>
#include <boost/io/ios_state.hpp>
#include <concepts>
#include <future>
#include <iomanip>
#include <unordered_map>

namespace NBody {

/// \brief The nbody program.
template<std::floating_point T>
class NBodyProgram {
 public:
  using ValueType = T;
  using Vector = Math::Vector<ValueType, 3>;
  using Particle = Simulation::Particle<ValueType>;
  using Particles = std::vector<Particle>;
  using ColorGenerator = Configuration::ColorGenerator<ValueType>;
  using Parser = Configuration::Parser<ValueType>;
  using Simulator = Simulation::Simulator<ValueType>;
  using InitialData = Simulation::InitialData<ValueType>;
  using Visualizer = Visualization::Visualizer<ValueType>;
  using InputMode = Configuration::InputMode;
  using OutputMode = Configuration::OutputMode;
  using Data = std::optional<Particles>;
  using DataQueue = Concurrency::Queue<Data>;

 private:
  const Parser parser_;

 public:
  /// \brief Constructs the nbody program with command line options.
  /// \param argc The count of arguments
  /// \param argv The vector of arguments
  NBodyProgram(int argc, char* argv[])
      : parser_{argc, argv} {}

  /// \brief Launches the program.
  void launch() const {
    // Set the number of threads will be used in OpenMP
    omp_set_num_threads(parser_.getNumThreads());
    switch (parser_.getInputMode()) {
      case InputMode::Frame:
        handleFrameInput();
        break;
      case InputMode::Simulation:
        handleSimulationInput();
        break;
    }
  }

 private:
  void handleFrameInput() const {
    auto inputFileStream = getInputFileStream();
    // Standard input stream will be used if input file stream is not provided.
    std::istream& istream =
        inputFileStream ? inputFileStream.value() : std::cin;

    // A cache used to store the color of each particle.
    // Since the raw data of particles do not contain color information,
    // NBody program will create a random color for each particle and store
    // the value in this cache for reuse.
    auto colorCache = [
        colorGenerator = makeColorGenerator(),
        colorMap = std::unordered_map<std::size_t, Vector>{}
    ](std::size_t particleIndex) mutable {
      if (!colorMap.contains(particleIndex)) {
        colorMap[particleIndex] = colorGenerator();
      }
      return colorMap[particleIndex];
    };

    auto readParticles = [&istream, &colorCache]() {
      std::size_t numParticles;
      istream >> numParticles;
      Particles particles;
      particles.reserve(numParticles);
      Particle particle;
      for (std::size_t i = 0; i < numParticles; ++i) {
        istream >> particle.index >> particle.radius >> particle.position;
        particle.color = colorCache(particle.index);
        particles.emplace_back(std::move(particle));
      }
      return particles;
    };

    auto saver = boost::io::basic_ios_exception_saver{istream};

    try {
      istream.exceptions(std::istream::failbit | std::istream::badbit);
    } catch (...) {
      throw std::runtime_error{"Error occurred in parsing input data."};
    }


    auto visConfig = parser_.getVisualizationConfiguration();
    try {
      istream >> visConfig.fps
              >> visConfig.ratio
              >> visConfig.cameraSpeedHint;
    } catch (...) {
      throw std::runtime_error{"Error occurred in parsing input data."};
    }

    // A queue used to buffer the data read from input.
    auto dataQueue = DataQueue{};

    auto fetchData = [&istream, &dataQueue, &readParticles]() {
      try {
        while (!istream.eof()) {
          dataQueue.push(readParticles());
          istream >> std::ws;
        }
        dataQueue.close();
      } catch (...) {
        dataQueue.close();
        throw std::runtime_error{"Error occurred in parsing input data."};
      }
    };

    if (parser_.getPreloading()) {
      fetchData();
      auto vis = Visualizer{visConfig, [&dataQueue](double) {
        auto data = Data{};
        dataQueue.pop(data);
        return data;
      }};
      vis.start();
    } else {
      auto inputThread = std::async(fetchData);
      auto vis = Visualizer{visConfig, [&](double) {
        auto data = Data{};
        if (dataQueue.pop(data) == DataQueue::Status::Closed) {
          inputThread.get();
        }
        return data;
      }};
      vis.start();
    }
  }

  void handleSimulationInput() const {
    const auto simConfig = parser_.getSimulationConfiguration();
    const auto initialData = getInitialData();
    auto sim = Simulator{simConfig, initialData};
    switch (parser_.getOutputMode()) {
      case OutputMode::Step:
        outputStepResult(sim);
        break;
      case OutputMode::Frame:
        outputFrameResult(sim);
        break;
      case OutputMode::Visualization: {
        visualizeResult(sim, initialData.boundaryBox.maxLength());
        break;
      }
    }
  }

  void outputStepResult(Simulator& sim) const {
    const auto outputFunction = [
        &sim,
        precision = parser_.getOutputPrecision()
    ](std::ostream& ostream) {
      const auto outputSingleStep = [=, &ostream](
          const std::uint64_t step,
          const Particles& particles
      ) {
        const auto saver = boost::io::ios_flags_saver{ostream};
        ostream << step << "\n" << particles.size() << "\n";
        ostream << std::scientific << std::setprecision(precision);
        for (auto it = particles.cbegin(); it != particles.cend(); ++it) {
          ostream << *it;
          if (std::next(it) != particles.cend()) {
            ostream << "\n";
          }
        }
      };

      ostream << sim.getTimeStep() << "\n";
      outputSingleStep(0, sim.getParticles());
      while (!sim.isFinished()) {
        ostream << "\n";
        outputSingleStep(sim.getCurrentStep(), sim.forwardStep(1));
      }
    };

    outputResult(outputFunction);
  }

  void outputFrameResult(Simulator& sim) const {
    const auto outputFunction = [
        &sim,
        precision = parser_.getOutputPrecision(),
        fps = parser_.getFps(),
        ratio = parser_.getRatio()
    ](std::ostream& ostream) {
      const auto outputSingleFrame = [=, &ostream](const Particles& particles) {
        const auto saver = boost::io::ios_flags_saver{ostream};
        ostream << particles.size() << "\n";
        ostream << std::scientific << std::setprecision(precision);
        for (auto it = particles.cbegin(); it != particles.cend(); ++it) {
          ostream << it->index << " " << it->radius << " " << it->position;
          if (std::next(it) != particles.cend()) {
            ostream << "\n";
          }
        }
      };

      ostream << fps << "\n";
      ostream << ratio << "\n";
      ostream << sim.getBoundaryBox().maxLength() << "\n";
      outputSingleFrame(sim.getParticles());
      while (!sim.isFinished()) {
        ostream << "\n";
        outputSingleFrame(sim.forwardTime(ratio / static_cast<ValueType>(fps)));
      }
    };

    outputResult(outputFunction);
  }

  void visualizeResult(
      Simulator& sim,
      ValueType cameraSpeedHint
  ) const {
    auto visConfig = parser_.getVisualizationConfiguration();
    visConfig.cameraSpeedHint = cameraSpeedHint;
    auto vis = Visualizer{visConfig, [&sim](double dt) {
      if (sim.isFinished()) {
        return Data{};
      } else {
        return Data{sim.forwardTime(static_cast<ValueType>(dt))};
      }
    }};
    vis.start();
  }

  [[nodiscard]] std::optional<std::ifstream> getInputFileStream() const {
    const auto inputDataPath = parser_.getInputDataPath();
    if (inputDataPath) {
      const auto path = inputDataPath.value();
      auto ifs = std::ifstream{path};
      if (!ifs) {
        throw std::runtime_error{"Cannot open file '" + path + "'."};
      }
      return std::make_optional(std::move(ifs));
    } else {
      return std::optional<std::ifstream>{};
    }
  }

  [[nodiscard]] std::optional<std::ofstream> getOutputFileStream() const {
    const auto outputPath = parser_.getOutputDataPath();
    if (outputPath) {
      const auto path = outputPath.value();
      auto ifs = std::ofstream{path};
      if (!ifs) {
        throw std::runtime_error{"Cannot open file '" + path + "'."};
      }
      return std::make_optional(std::move(ifs));
    } else {
      return std::optional<std::ofstream>{};
    }
  }

  [[nodiscard]] InitialData getInitialData() const {
    auto inputFileStream = getInputFileStream();
    // Standard input stream will be used if input file stream is not provided.
    std::istream& istream =
        inputFileStream ? inputFileStream.value() : std::cin;

    auto saver = boost::io::basic_ios_exception_saver{istream};

    try {
      istream.exceptions(std::istream::failbit | std::istream::badbit);
    } catch (...) {
      throw std::runtime_error{"Error occurred in parsing input data."};
    }

    auto initialData = InitialData{};
    try {
      istream >> initialData;
    } catch (...) {
      throw std::runtime_error{"Error occurred in parsing input data."};
    }

    const auto& boundaryBox = initialData.boundaryBox;
    const auto& particles = initialData.particles;
    auto particleIndices = std::vector<std::size_t>{};
    for (std::size_t i = 0; i < particles.size(); ++i) {
      particleIndices.emplace_back(i);
    }

    // Checks whether the boundaryBox is valid.
    if (!boundaryBox.isValid()) {
      throw std::runtime_error("Found invalid boundary box input.");
    }

    // Check whether all particles are within the boundary box.
    using OpenBoundaryConditionHandler =
    Simulation::OpenBoundaryConditionHandler<ValueType>;
    const auto passedParticleIndices =
        OpenBoundaryConditionHandler::computePassedParticleIndices(
            particles,
            particleIndices,
            boundaryBox
        );
    if (!passedParticleIndices.empty()) {
      auto oss = std::ostringstream{};
      for (const auto& index: passedParticleIndices) {
        oss << particles[index].index << "\n";
      }
      throw std::runtime_error("Found outside particles:\n" + oss.str());
    }

    // Checks whether all particles are not collided with others.
    using DirectCollisionDetector =
    Simulation::DirectCollisionDetector<ValueType>;
    const auto collidedParticleIndexPairs =
        DirectCollisionDetector::computeCollidedParticleIndexPairs(
            particles,
            particleIndices
        );
    if (!collidedParticleIndexPairs.empty()) {
      auto oss = std::ostringstream{};
      for (const auto& pair: collidedParticleIndexPairs) {
        oss << "[" << particles[pair.first].index << ", "
            << particles[pair.second].index << "]\n";
      }
      throw std::runtime_error("Found collided particle pairs:\n" + oss.str());
    }

    // The raw data do not contain color information, so the color generator
    // will create the color for each particle.
    auto colorGenerator = makeColorGenerator();
    for (auto& particle: initialData.particles) {
      particle.color = colorGenerator();
    }
    return initialData;
  }

  void outputResult(
      const std::function<void(std::ostream&)>& outputFunction
  ) const {
    auto outputFileStream = getOutputFileStream();
    std::ostream& ostream =
        outputFileStream ? outputFileStream.value() : std::cout;

    auto saver = boost::io::basic_ios_exception_saver{ostream};
    try {
      outputFunction(ostream);
    } catch (...) {
      throw std::runtime_error{"Error occurred in outputting data."};
    }
  }

  [[nodiscard]] ColorGenerator makeColorGenerator() const {
    return ColorGenerator{parser_.getColorful() ? 0.1f : 0.9f, 0.9f};
  }
};

} // namespace NBody

#endif //CPP_PROJECT_INCLUDE_NBODY_NBODYPROGRAM_HPP_
