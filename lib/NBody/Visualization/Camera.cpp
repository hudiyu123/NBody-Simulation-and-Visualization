#include <NBody/Visualization/Camera.hpp>
#include <NBody/Math/Common.hpp>
#include <cmath>

namespace NBody::Visualization {

Camera::Camera(
    ValueType movementSpeed,
    Vector3 position,
    Vector3 up,
    ValueType yaw,
    ValueType pitch
) : movementSpeed_{movementSpeed},
    position_{position},
    worldUp_{up},
    up_{},
    right_{},
    yaw_{yaw},
    pitch_{pitch},
    front_{
        static_cast<ValueType>(0.0f),
        static_cast<ValueType>(0.0f),
        static_cast<ValueType>(-1.0f)
    } {
  updateCameraVectors();
}

Camera::Matrix4 Camera::getViewMatrix() const {
  return Math::lookAt(position_, position_ + front_, up_);
}

void Camera::process(
    const std::set<CameraOperationType>& operations,
    double deltaTime
) {
  bool shift = operations.contains(CameraOperationType::Shift);
  bool control = operations.contains(CameraOperationType::Control);

  auto cameraMoveVelocity = movementSpeed_ * static_cast<ValueType>(deltaTime);
  auto cameraSpinVelocity = static_cast<ValueType>(0.3f);
  if (shift) {
    cameraMoveVelocity *= static_cast<ValueType>(5.0f);
    cameraSpinVelocity *= static_cast<ValueType>(3.0f);
  }
  if (control) {
    cameraMoveVelocity /= static_cast<ValueType>(5.0f);
    cameraSpinVelocity /= static_cast<ValueType>(3.0f);
  }

  for (const auto& operation: operations) {
    switch (operation) {
      case CameraOperationType::MoveForward:
        position_ += front_ * cameraMoveVelocity;
        break;
      case CameraOperationType::MoveBackward:
        position_ -= front_ * cameraMoveVelocity;
        break;
      case CameraOperationType::TurnLeft:
        position_ -= right_ * cameraMoveVelocity;
        break;
      case CameraOperationType::TurnRight:
        position_ += right_ * cameraMoveVelocity;
        break;
      case CameraOperationType::YawLeft:
        yaw_ -= cameraSpinVelocity;
        break;
      case CameraOperationType::YawRight:
        yaw_ += cameraSpinVelocity;
        break;
      case CameraOperationType::PitchUp:
        pitch_ += cameraSpinVelocity;
        break;
      case CameraOperationType::PitchDown:
        pitch_ -= cameraSpinVelocity;
        break;
      default:
        break;
    }

    if (yaw_ < static_cast<ValueType>(-180)) {
      yaw_ += static_cast<ValueType>(360);
    } else if (yaw_ > static_cast<ValueType>(180)) {
      yaw_ -= static_cast<ValueType>(360);
    }

    if (pitch_ > static_cast<ValueType>(89.9)) {
      pitch_ = static_cast<ValueType>(89.9);
    } else if (pitch_ < static_cast<ValueType>(-89.9)) {
      pitch_ = static_cast<ValueType>(-89.9);
    }
  }
  updateCameraVectors();
}

void Camera::updateCameraVectors() {
  Vector3 front{
      std::cos(Math::radians(yaw_)) * std::cos(Math::radians(pitch_)),
      std::sin(Math::radians(pitch_)),
      std::sin(Math::radians(yaw_)) * std::cos(Math::radians(pitch_))
  };
  front_ = Math::normalize(front);
  right_ = Math::normalize(Math::cross(front_, worldUp_));
  up_ = Math::normalize(Math::cross(right_, front_));
}

} // namespace NBody::Visualization