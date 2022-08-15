#ifndef CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_CAMERA_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_CAMERA_HPP_

#include <NBody/Math/Vector.hpp>
#include <NBody/Math/Matrix.hpp>
#include <NBody/Visualization/GLLoader.hpp>
#include <array>
#include <set>

namespace NBody::Visualization {

/// \brief The type of camera operation.
enum class CameraOperationType {
  /// \brief Move forward.
  MoveForward,
  /// \brief Move backward.
  MoveBackward,
  /// \brief Turn left.
  TurnLeft,
  /// \brief Turn right.
  TurnRight,
  /// \brief Yaw left.
  YawLeft,
  /// \brief Yaw right.
  YawRight,
  /// \brief Pitch up.
  PitchUp,
  /// \brief Pitch down
  PitchDown,
  /// \brief Shift key is pressed.
  Shift,
  /// \brief Control key is pressed.
  Control
};

class Camera {
 public:
  using ValueType = GLfloat;
  using Vector3 = Math::Vector<ValueType, 3>;
  using Matrix4 = Math::Matrix<ValueType, 4, 4>;

 private:
  Vector3 position_;
  Vector3 front_;
  Vector3 worldUp_;
  Vector3 up_;
  Vector3 right_;
  ValueType yaw_;
  ValueType pitch_;
  ValueType movementSpeed_;

 public:
  /// \brief Constructs a camera with parameters.
  /// \param movementSpeed The movement speed of the camera
  /// \param position The initial position of the camera
  /// \param up The up vector
  /// \param yaw The initial yaw angle in degrees
  /// \param pitch The initial pitch angle in degrees
  explicit Camera(
      ValueType movementSpeed,
      Vector3 position = {
          static_cast<ValueType>(0.0f),
          static_cast<ValueType>(0.0f),
          static_cast<ValueType>(0.0f)
      },
      Vector3 up = {
          static_cast<ValueType>(0.0f),
          static_cast<ValueType>(1.0f),
          static_cast<ValueType>(0.0f)
      },
      ValueType yaw = static_cast<ValueType>(-90.0f),
      ValueType pitch = static_cast<ValueType>(0.0f)
  );

  /// \brief Returns the view matrix of the camera.
  /// \return The view matrix
  [[nodiscard]] Matrix4 getViewMatrix() const;

  /// \brief Processes the give operations with delta time.
  /// \param operations The operations
  /// \param deltaTime The delta time
  void process(
      const std::set<CameraOperationType>& operations,
      double deltaTime
  );

 private:
  void updateCameraVectors();
};

} // namespace NBody::Visualization

#endif // CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_CAMERA_HPP_
