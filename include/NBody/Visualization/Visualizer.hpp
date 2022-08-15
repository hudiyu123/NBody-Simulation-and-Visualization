#ifndef CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VISUALIZER_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VISUALIZER_HPP_

#include <NBody/Visualization/GLLoader.hpp>
#include <NBody/Visualization/Configuration.hpp>
#include <NBody/Visualization/Camera.hpp>
#include <NBody/Visualization/Icosphere.hpp>
#include <NBody/Visualization/VertexArray.hpp>
#include <NBody/Visualization/VertexBuffer.hpp>
#include <NBody/Visualization/VertexBufferLayout.hpp>
#include <NBody/Visualization/IndexBuffer.hpp>
#include <NBody/Visualization/Shader.hpp>
#include <NBody/Visualization/Renderer.hpp>
#include <NBody/Simulation/Particle.hpp>
#include <GLFW/glfw3.h>
#include <chrono>
#include <concepts>
#include <functional>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace NBody::Visualization {

#ifdef __APPLE__
/// \brief The vertex shader source used by visualizer.
constexpr std::basic_string_view<GLchar> vShaderSource = R"(
#version 330

in vec3 aPosition;

uniform float radius;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

mat4 scale(float x, float y, float z) {
  return mat4(
    vec4(x, 0.0f, 0.0f, 0.0f),
    vec4(0.0f, y, 0.0f, 0.0f),
    vec4(0.0f, 0.0f, z, 0.0f),
    vec4(0.0f, 0.0f, 0.0f, 1.0f)
  );
}

void main() {
  mat4 scaleMatrix = scale(radius, radius, radius);
  gl_Position = projection * view * model * scaleMatrix * vec4(aPosition, 1.0f);
}
)";

/// \brief The fragment shader source used by visualizer.
constexpr std::basic_string_view<GLchar> fShaderSource = R"(
#version 330
out vec4 fColor;
uniform vec4 customColor;
void main() {
  fColor = customColor;
}
)";
#else
/// \brief The vertex shader source used by visualizer.
constexpr std::basic_string_view<GLchar> vShaderSource = R"(
#version 120

attribute vec3 aPosition;

uniform float radius;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

mat4 scale(float x, float y, float z) {
  return mat4(
    vec4(x, 0.0f, 0.0f, 0.0f),
    vec4(0.0f, y, 0.0f, 0.0f),
    vec4(0.0f, 0.0f, z, 0.0f),
    vec4(0.0f, 0.0f, 0.0f, 1.0f)
  );
}

void main() {
  mat4 scaleMatrix = scale(radius, radius, radius);
  gl_Position = projection * view * model * scaleMatrix * vec4(aPosition, 1.0f);
}
)";

/// \brief The fragment shader source used by visualizer.
constexpr std::basic_string_view<GLchar> fShaderSource = R"(
#version 120
uniform vec4 customColor;
void main() {
  gl_FragColor = customColor;
}
)";
#endif

/// \brief This class is used to visualize the nbody simulation results.
template<std::floating_point T>
class Visualizer {
 public:
  using ValueType = T;
  using Particle = Simulation::Particle<ValueType>;
  using Particles = std::vector<Particle>;
  /// \brief The type if the return value from the simulation callback.
  /// \details This optional will be empty if the simulation is finished or
  /// stopped.
  using CallbackResult = std::optional<Particles>;
  /// \brief The type of simulation callback.
  /// \details The parameter deltaTime indicates how long time has passed since
  /// the last call.
  using SimulationCallback = std::function<CallbackResult(double deltaTime)>;

 private:
  const Icosphere<GLfloat, GLuint, GLsizeiptr, GLsizei> sphere_;
  Camera camera_;
  const std::uint64_t fps_;
  const double ratio_;
  int windowWidth_;
  int windowHeight_;
  const SimulationCallback simulationCallback_;
  bool isPaused_;
  bool isFinished_;
  int lastSpaceKayAction_;
  const std::basic_string<GLchar> vShaderSource_;
  const std::basic_string<GLchar> fShaderSource_;
  std::vector<Particle> currentParticles_;

 public:
  /// \brief Constructs a visualizer with configuration and the data callback.
  /// \param config The configuration
  /// \param callback The data callback
  Visualizer(const Configuration& config, SimulationCallback&& callback)
      : sphere_{static_cast<GLfloat>(1.0f), config.icosphereSubdivision},
        camera_{
            static_cast<GLfloat>(config.cameraSpeedHint / 10.0f),
            Math::Vector<GLfloat, 3>{
                static_cast<GLfloat>(0.0f),
                static_cast<GLfloat>(0.0f),
                static_cast<GLfloat>(config.cameraPositionHint)
            }
        },
        fps_{config.fps},
        ratio_{config.ratio},
        windowWidth_{config.windowWidth},
        windowHeight_{config.windowHeight},
        simulationCallback_{callback},
        isPaused_{false},
        isFinished_{false},
        lastSpaceKayAction_{GLFW_RELEASE},
        vShaderSource_{vShaderSource},
        fShaderSource_{fShaderSource},
        currentParticles_{} {}

  /// \brief Starts the visualization.
  void start() {
    glfwSetErrorCallback([](int errorCode, const char* description) {
      auto oss = std::ostringstream{};
      oss << "Error code " << errorCode << ": " << description << std::endl;
      throw std::runtime_error(oss.str());
    });
    if (!glfwInit()) {
      throw std::runtime_error("GLFW initialize failed.");
    }
    GLFWwindow* window = makeWindow(
        windowWidth_,
        windowHeight_,
        "N-Body Simulation"
    );
    if (!window) {
      throw std::runtime_error("Create window failed.");
    }
    // Set the user pointer in order to retrieve the visualizer in callback.
    glfwSetWindowUserPointer(window, this);

    auto framebufferSizeCallback = [](
        GLFWwindow* window,
        int width,
        int height
    ) {
      auto vis = static_cast<Visualizer*>(glfwGetWindowUserPointer(window));
      vis->framebufferSizeCallback(width, height);
    };
    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);

    auto keyCallback = [](
        GLFWwindow* window,
        int key,
        int scancode,
        int action,
        int mods
    ) {
      auto vis = static_cast<Visualizer*>(glfwGetWindowUserPointer(window));
      vis->keyCallback(window, key, action);
    };
    glfwSetKeyCallback(window, keyCallback);

    glfwMakeContextCurrent(window);
#ifdef __APPLE__
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
      throw std::runtime_error("GLEW initialize failed.");
    }
#else
    if (!glLoaderInit()) {
      throw std::runtime_error("GLLoader initialize failed.");
    }
#endif

    glEnable(GL_DEPTH_TEST);

    auto shader = Shader{vShaderSource_, fShaderSource_, {{0, "aPosition"}}};
    auto va = VertexArray{};
    auto vb = VertexBuffer{
        sphere_.getRawVertices().data(),
        sphere_.getRawVerticesSize(),
        GL_STATIC_DRAW
    };
    auto layout = VertexBufferLayout{};
    layout.pushGLfloat(3);
    va.addVertexBuffer(vb, layout);
    auto ib = IndexBuffer{
        sphere_.getVertexIndices().data(),
        sphere_.getVertexIndicesSize(),
        sphere_.getVertexIndicesCount(),
        GL_UNSIGNED_INT,
        GL_STATIC_DRAW
    };

    double frameTime = 1.0 / static_cast<double>(fps_);
    double lastTime = glfwGetTime();
    while (!glfwWindowShouldClose(window)) {
      const double currentTime = glfwGetTime();
      const double deltaTime = currentTime - lastTime;
      lastTime = currentTime;

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      camera_.process(getCameraOperations(window), deltaTime);
      const auto view = camera_.getViewMatrix();
      const auto width = static_cast<GLfloat>(windowWidth_);
      const auto height = static_cast<GLfloat>(windowHeight_);
      const auto projection = Math::infinitePerspective(
          Math::radians(static_cast<GLfloat>(45.0f)),
          width / height,
          static_cast<GLfloat>(0.1f)
      );
      shader.bind();
      // Sets the projection uniform.
      shader.setUniformMatrix4fv("projection", &projection[0][0], true);
      // Sets the view uniform
      shader.setUniformMatrix4fv("view", &view[0][0], true);

      if (!isPaused_ && !isFinished_) {
        auto result = simulationCallback_(frameTime * ratio_);
        if (result) {
          currentParticles_ = std::move(result.value());
        } else {
          isFinished_ = true;
        }
      }

      // Change the window title to show current simulation status.
      if (isFinished_) {
        glfwSetWindowTitle(window, "NBody Simulation (Finished)");
      } else if (isPaused_) {
        glfwSetWindowTitle(window, "NBody Simulation (Paused)");
      } else {
        glfwSetWindowTitle(window, "NBody Simulation");
      }

      for (const auto& particle: currentParticles_) {
        // Sets the color uniform for each particle.
        shader.setUniform4f(
            "customColor",
            static_cast<GLfloat>(particle.color[0]),
            static_cast<GLfloat>(particle.color[1]),
            static_cast<GLfloat>(particle.color[2]),
            static_cast<GLfloat>(1)
        );
        // Sets the radius uniform for each particle.
        shader.setUniform1f("radius", static_cast<GLfloat>(particle.radius));
        // Computes the model matrix with the position.
        auto model = Math::Matrix<GLfloat, 4, 4>{static_cast<GLfloat>(1.0f)};
        model[0][3] = static_cast<GLfloat>(particle.position[0]);
        model[1][3] = static_cast<GLfloat>(particle.position[1]);
        model[2][3] = static_cast<GLfloat>(particle.position[2]);
        // Sets the model uniform for each particle.
        shader.setUniformMatrix4fv("model", &model[0][0], true);
        // Draw each particle separately.
        Renderer::draw(va, ib, shader);
      }

      glfwSwapBuffers(window);
      glfwPollEvents();

      double elapsedTime = glfwGetTime() - currentTime;
      double gap = frameTime - elapsedTime;
      // Sleep until current frame time is over.
      if (gap > 0.0) {
        const auto microGap = static_cast<long>(gap * 1000000);
        std::this_thread::sleep_for(std::chrono::microseconds{microGap});
      }
    }
    glfwTerminate();
  }

 private:
  static std::set<CameraOperationType> getCameraOperations(GLFWwindow* window) {
    auto operations = std::set<CameraOperationType>{};
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
      operations.insert(CameraOperationType::MoveForward);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
      operations.insert(CameraOperationType::MoveBackward);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
      operations.insert(CameraOperationType::TurnLeft);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
      operations.insert(CameraOperationType::TurnRight);
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS) {
      operations.insert(CameraOperationType::YawLeft);
    }
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS) {
      operations.insert(CameraOperationType::YawRight);
    }
    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {
      operations.insert(CameraOperationType::PitchUp);
    }
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {
      operations.insert(CameraOperationType::PitchDown);
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) {
      operations.insert(CameraOperationType::Shift);
    }
    if (glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS) {
      operations.insert(CameraOperationType::Shift);
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS) {
      operations.insert(CameraOperationType::Control);
    }
    if (glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS) {
      operations.insert(CameraOperationType::Control);
    }
    return operations;
  }

  static GLFWwindow* makeWindow(
      int width,
      int height,
      const std::string& title
  ) {
#ifdef __APPLE__
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#else
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
#endif
    GLFWwindow* window = glfwCreateWindow(
        width,
        height,
        title.c_str(),
        nullptr,
        nullptr
    );
    return window;
  }

  void framebufferSizeCallback(int width, int height) {
    windowWidth_ = width;
    windowHeight_ = height;
    glViewport(0, 0, width, height);
  }

  void keyCallback(GLFWwindow* window, int key, int action) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
      glfwSetWindowShouldClose(window, true);
    }
    if (key == GLFW_KEY_SPACE) {
      if (lastSpaceKayAction_ == GLFW_RELEASE && action == GLFW_PRESS) {
        isPaused_ = !isPaused_;
      }
      lastSpaceKayAction_ = action;
    }
  }
};

} // namespace NBody::Visualization

#endif // CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_VISUALIZER_HPP_
