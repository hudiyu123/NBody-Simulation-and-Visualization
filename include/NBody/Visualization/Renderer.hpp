#ifndef CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_RENDERER_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_RENDERER_HPP_

#include <NBody/Visualization/VertexArray.hpp>
#include <NBody/Visualization/IndexBuffer.hpp>
#include <NBody/Visualization/Shader.hpp>

namespace NBody::Visualization {

/// \brief The abstraction of draw process of OpenGL.
class Renderer {
 public:
  /// \brief Draws to the current OpenGL context with vertex array, index buffer
  /// and shader.
  /// \param va The vertex array
  /// \param ib The index buffer
  /// \param shader The shader
  static void draw(
      const VertexArray& va,
      const IndexBuffer& ib,
      const Shader& shader
  );
};

} // NBody::Visualization

#endif // CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_RENDERER_HPP_
